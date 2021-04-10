//
//  GLLPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  GLL point for preloop processing
//  generator of Point and boundary conditions in core

#ifndef GLLPoint_hpp
#define GLLPoint_hpp

#include "eigen_sem.hpp"
#include "window_tools.hpp"
#include <map>
#include <vector>
#include <memory>

// release
class ABC;
class TimeScheme;
class Domain;
class Point;
class SolidPointWindow;
class FluidPointWindow;

typedef std::tuple<int, int, double> winmap;

class GLLPoint {
public:
    // setup by Quad
    void setup(int globalTag, const eigen::DCol2 &coords, 
               double distTol) {
        if (mGlobalTag == -1) {
            // set by the first element
            mGlobalTag = globalTag;
            mCoords = coords;
            mAngleTol = distTol / mCoords(0);
        } else {
            // check subsequent elements
            if (mGlobalTag != globalTag ||
                (mCoords - coords).norm() > distTol) {
                throw std::runtime_error("GLLPoint::setup || "
                                         "Conflict in GLL-point setup.");
            }  
        }
        // reference count
        mElementCount++;
    }
    
    // add window
    eigen::IColX addWindows(const std::vector<eigen::DMatX2> &wins) {      
        eigen::IColX tags = eigen::IColX::Constant(wins.size(), -1);
        
        // first set of windows added
        if (mWindows.size() == 0) {
            for (int mn = 0; mn < wins.size(); mn++) {
                tags(mn) = mWindows.size();
                mWindows.push_back(wins[mn]);
            }
            mNumNewWins = mWindows.size();
            allocateNewWindows();
            mNumNewWins = 0;
            return tags;
        }
        
        // pre-existing windows
        std::vector<bool> found(mWindows.size(), false);
        for (int mn = 0; mn < wins.size(); mn++) {
            for (int mo = 0; mo < mWindows.size(); mo++) {
                if (std::abs(mWindows[mo](0, 0) - wins[mn](0, 0)) < mAngleTol &&
                    std::abs(mWindows[mo](mWindows[mo].rows() - 1, 0) - wins[mn](mWindows[mo].rows() - 1, 0))  < mAngleTol) {
                    if (wins[mn].rows() != mWindows[mo].rows()) {
                        throw std::runtime_error("GLLPoint::addWindows || "
                                         "Conflict in nr window setup.");
                    }
                    found[mo] = true;
                    tags(mn) = mo;
                    break;
                }
            }
        }
        
        // all windows equal
        if (tags.minCoeff() >= 0) return tags;
        mNumNewWins = (tags.array() < 0).count();
        
        // differing windows -> create map with intepolants for exchanging 3d material/geometry info between window sets
        
        // collect azimuthal samples in existing windows
        std::map<double, winmap> angle_map_old;
        for (int mo = 0; mo < mWindows.size(); mo++) {
            if (found[mo]) continue;
            for (int alpha = 0; alpha < mWindows[mo].rows(); alpha++) {
                angle_map_old.insert(std::pair<double, winmap>(mWindows[mo](alpha, 0), {mo, alpha, 0.}));
            }
        }
        
        // create new windows and
        // collect azimuthal samples in new windows
        std::map<double, winmap> angle_map_new;
        for (int mn = 0; mn < wins.size(); mn++) {
            if (tags(mn) >= 0) continue;
            mWindows.push_back(wins[mn]);
            tags(mn) = mWindows.size();
            found.push_back(false);
            for (int alpha = 0; alpha < wins[mn].rows(); alpha++) {
                angle_map_new.insert(std::pair<double, winmap>(wins[mn](alpha, 0), {tags(mn), alpha, 0.}));
            }
        }
        
        // create map and linear interpolants
        mMapWindows.clear();
        mMapWindows.resize(mWindows.size());
        for (int m = 0; m < mWindows.size(); m++) {
            if (found[m]) continue;
            mMapWindows[m].resize(mWindows[m].rows());
            for (int alpha = 0; alpha < mWindows[m].rows(); alpha++) {
                if (m < mWindows.size() - mNumNewWins) { // from new wins to previous wins
                    mMapWindows[m][alpha] = linearInterpWinMap(angle_map_new, 
                                            mWindows[m](alpha, 0));
                } else { // from previous wins to new wins
                    mMapWindows[m][alpha] = linearInterpWinMap(angle_map_old, 
                                            mWindows[m](alpha, 0));
                }
            }
        }
        allocateNewWindows();
        return tags;
    }
    
    std::pair<winmap, winmap> linearInterpWinMap(const std::map<double, winmap> &map, double value) {
        if (map.size() == 1) {
            winmap onlymap = map.begin()->second;
            std::get<2>(onlymap) = 1.0;
            return std::pair<winmap, winmap>(onlymap, map.begin()->second);
        }
        
        std::pair<winmap, winmap> contributors;
        auto it = map.upper_bound(value);
        double val1, val2;
        if (it == map.begin()) { 
            contributors = std::pair(map.rbegin()->second, map.begin()->second);
            val1 = map.rbegin()->first - 2 * numerical::dPi;
            val2 = map.begin()->first;
        } else if (it == map.end()) {
            contributors = std::pair(map.rbegin()->second, map.begin()->second);
            val1 = map.rbegin()->first;
            val2 = map.begin()->first + 2 * numerical::dPi;
        } else {
            auto it_prev = it;
            it_prev--;
            contributors = std::pair(it_prev->second, it->second);
            val1 = it_prev->first;
            val2 = it->first;
        }

        double factor = (value - val1) / (val2 - val1);
        std::get<2>(contributors.first) = 1. - factor;
        std::get<2>(contributors.second) = factor;
        
        return contributors;
    }
    
    void allocateNewWindows() {
        if (mNumNewWins == 0) return;

        mMassFluid.insert(mMassFluid.end(), mNumNewWins, eigen::DColX::Zero(0));
        mMassSolid.insert(mMassSolid.end(), mNumNewWins, eigen::DColX::Zero(0));
        mNormalSFU.insert(mNormalSFU.end(), mNumNewWins, eigen::DMatX3::Zero(0, 3));
        mNormalSFA.insert(mNormalSFA.end(), mNumNewWins, eigen::DMatX3::Zero(0, 3));
        mClaytonABC.resize(mWindows.size());
        mGamma.insert(mGamma.end(), mNumNewWins, eigen::DColX::Zero(0));
        mCountGammasAdded.insert(mCountGammasAdded.end(), mNumNewWins, 0);
        mNormalTop.insert(mNormalTop.end(), mNumNewWins, eigen::DMatX3::Zero(0, 3));
        mSumRhoDepth.insert(mSumRhoDepth.end(), mNumNewWins, eigen::DColX::Zero(0));
    }
    
    void exchange3DInfoBetweenWindows() {
        if (mNumNewWins == 0) return;
        
        // interpolate between window sets using map
        std::vector<eigen::DColX> unchanged_property1;
        std::vector<eigen::DMatX3> unchanged_property3;
        unchanged_property1 = mMassFluid;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mMapWindows[m].size() == 0) continue;
            addInfo(mMassFluid, unchanged_property1, m, mMapWindows[m]);
        }
        
        unchanged_property1 = mMassSolid;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mMapWindows[m].size() == 0) continue;
            addInfo(mMassSolid, unchanged_property1, m, mMapWindows[m]);
        }
        
        unchanged_property3 = mNormalSFA;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mMapWindows[m].size() == 0) continue;
            bool set = addInfo(mNormalSFA, unchanged_property3, m, mMapWindows[m]);
            if (set) mNormalSFU[m] = mNormalSFA[m];
        }
        
        unchanged_property1 = mGamma;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mMapWindows[m].size() == 0) continue;
            bool set = addInfo(mGamma, unchanged_property1, m, mMapWindows[m]);
            if (set) mCountGammasAdded[m]++;
        }

        unchanged_property3 = mNormalTop;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mMapWindows[m].size() == 0) continue;
            addInfo(mNormalTop, unchanged_property3, m, mMapWindows[m]);
        }
        
        // finished window processing
        mNumNewWins = 0;
    }
    
    void add3DInfoFromWindowsOnOtherRanks() {
        if (mNumNewWins == 0) return;
        
        // add info to existing windows
        std::vector<eigen::DColX> empty_vector1;
        std::vector<eigen::DMatX3> empty_vector3;
        for (int m = 0; m < mWindows.size() - mNumNewWins; m++) {
            if (mMapWindows[m].size() == 0) continue;
            addInfo(mMassFluid, empty_vector1, m, mMapWindows[m]);
            addInfo(mMassSolid, empty_vector1, m, mMapWindows[m]);
            addInfo(mNormalSFA, empty_vector3, m, mMapWindows[m]);
            bool set = addInfo(mGamma, empty_vector1, m, mMapWindows[m]);
            if (set) mCountGammasAdded[m] += mCountGammasAdded[std::get<0>(mMapWindows[m][0].first)];
            addInfo(mNormalTop, empty_vector3, m, mMapWindows[m]);
        }
        
        // discard windows on other ranks
        for (int m = mWindows.size() - mNumNewWins; m < mWindows.size(); m++) {
            eigen::DCol2 startAndSpacing;
            startAndSpacing << mWindows[m](0, 0), mWindows[m](1, 0) - mWindows[m](0, 0);
            mWindowsFromOtherRanks.push_back(startAndSpacing);
        }
        mWindows.resize(mWindows.size() - mNumNewWins);
        mMassFluid.resize(mWindows.size());
        mMassSolid.resize(mWindows.size());
        mNormalSFA.resize(mWindows.size());
        mGamma.resize(mWindows.size());
        mCountGammasAdded.resize(mWindows.size());
        mNormalTop.resize(mWindows.size());
        mNumNewWins = 0;
    }
    
    template <class Mat>
    bool addInfo(std::vector<Mat> &property, 
                 std::vector<Mat> &unchanged_property, 
                 const int &m, const std::vector<std::pair<winmap, winmap>> &map) {
                   
        std::vector<Mat> *prop;
        if (unchanged_property.size() == 0) {
            prop = &property;
        } else {
            prop = &unchanged_property;
        }
        eigen::DMatXX property_add = eigen::DMatXX::Zero(map.size(), (*prop)[m].cols());
        bool property_set = false;
        for (int alpha = 0; alpha < map.size(); alpha++) {
            int m1 = std::get<0>(map[alpha].first);
            int m2 = std::get<0>(map[alpha].second);
            if ((*prop)[m1].rows() > 0) {
                property_set = true;
                if ((*prop)[m1].rows() == 1) {
                    property_add.row(alpha) += (*prop)[m1].row(0) * std::get<2>(map[alpha].first);
                } else {
                    property_add.row(alpha) += (*prop)[m1].row(std::get<1>(map[alpha].first)) * std::get<2>(map[alpha].first);
                }
            }
            
            if ((*prop)[m2].rows() > 0) {
                property_set = true;
                if ((*prop)[m2].rows() == 1) {
                    property_add.row(alpha) += (*prop)[m2].row(0) * std::get<2>(map[alpha].second);
                } else {
                    property_add.row(alpha) += (*prop)[m2].row(std::get<1>(map[alpha].second)) * std::get<2>(map[alpha].second);
                }
            }
        }
        if (!property_set) return false;
        op1D_3D::tryReduceTo1D(property_add);
        op1D_3D::addTo(property_add, property[m]);
        
        return true;
    }
    
    // add mass
    void addMass(int tag, const eigen::DColX &mass, bool fluid) {
        eigen::DColX &myMass = fluid ? mMassFluid[tag] : mMassSolid[tag];
        op1D_3D::addTo(mass, myMass);
    }
    
    // add solid-fluid
    void addNormalSF(int tag, const eigen::DMatX3 &nSF) {
        op1D_3D::addTo(nSF, mNormalSFU[tag]);
        mNormalSFA[tag] = mNormalSFU[tag];
    }
    
    // add Clayton ABC
    // NOTE: the normal is NOT assembled because velocity and density
    //       can be discontinuous; at a boundary point shared by M
    //       elements, M Clayton instances will be created
    void addClaytonABC(int tag, const std::string &key, bool fluid,
                       const eigen::DMatX3 &nABC, const eigen::DColX &rho,
                       const eigen::DColX &vp, const eigen::DColX &vs) {
        eigen::DColX rhoVp, rhoVs;
        op1D_3D::times(rho, vp, rhoVp);
        op1D_3D::times(rho, vs, rhoVs);
        // init (key, vector<tuple>)
        mClaytonABC[tag].insert({key, std::vector<
            std::tuple<bool, eigen::DMatX3, eigen::DColX, eigen::DColX>>()});
        mClaytonABC[tag].at(key).push_back({fluid, nABC, rhoVp, rhoVs});
        // try reduce to 1D
        op1D_3D::tryReduceTo1D(std::get<1>(mClaytonABC[tag].at(key).back()));
        op1D_3D::tryReduceTo1D(std::get<2>(mClaytonABC[tag].at(key).back()));
        op1D_3D::tryReduceTo1D(std::get<3>(mClaytonABC[tag].at(key).back()));
    }
    
    // add gamma
    void addGamma(int tag, const eigen::DColX &Gamma) {
        op1D_3D::addTo(Gamma, mGamma[tag]);
        mCountGammasAdded[tag]++;
    }
    
    // add ocean load
    void addOceanLoad(int tag, const eigen::DMatX3 &nTop,
                      const eigen::DColX &sumRhoDepth) {
        op1D_3D::addTo(nTop, mNormalTop[tag]);
        // this is done repeatedly at shared points
        mSumRhoDepth[tag] = sumRhoDepth;
    }
    
    // set axial
    void setAxial() {
        mAxial = true;
    }
    
    // set surface
    void setSurface() {
        mSurface = true;
    }
    
    // comm size
    int getCommSize() const {
        // must use this full size because the compact size
        // of recv buffer is unknown
        int sumNr = 0;
        for (auto &win: mWindows) sumNr += win.size();
        return 1 + 1 + mWindows.size() * (4 + 1) + mWindows.size() * 6 + sumNr + sumNr + sumNr * 3 + sumNr * 3 + sumNr;
    }
    
    // feed comm
    void feedComm(eigen::DColX &buffer, int &row) const {
        // reference element count
        buffer(row, 0) = mElementCount;
        
        // window info
        buffer(row + 1, 0) = mWindows.size();
        row += 2;
        for (auto &win: mWindows) {
            buffer(row + 0, 0) = win(0, 0);
            buffer(row + 1, 0) = win(win.rows() - 1, 0);
            buffer(row + 2, 0) = win.rows();
            row += 3;
        }
        
        for (int m = 0; m < mWindows.size(); m++) {
            // size info
            // sizes must be sent as they can have different sizes on ranks
            buffer(row + 0, 0) = mMassFluid[m].size();
            buffer(row + 1, 0) = mMassSolid[m].size();
            buffer(row + 2, 0) = mNormalSFA[m].size();
            buffer(row + 3, 0) = mNormalTop[m].size();
            buffer(row + 4, 0) = mGamma[m].size();
            buffer(row + 5, 0) = mCountGammasAdded[m];
            row += 6;
            // mass
            buffer.block(row, 0, mMassFluid[m].size(), 1) = mMassFluid[m];
            row += mMassFluid[m].size();
            buffer.block(row, 0, mMassSolid[m].size(), 1) = mMassSolid[m];
            row += mMassSolid[m].size();
            // solid-fluid
            buffer.block(row, 0, mNormalSFA[m].size(), 1) =
            Eigen::Map<const eigen::DColX>(mNormalSFA[m].data(), mNormalSFA[m].size(), 1);
            row += mNormalSFA[m].size();
            // ocean load
            buffer.block(row, 0, mNormalTop[m].size(), 1) =
            Eigen::Map<const eigen::DColX>(mNormalTop[m].data(), mNormalTop[m].size(), 1);
            row += mNormalTop[m].size();
            // sponge boundary
            buffer.block(row, 0, mGamma[m].size(), 1) = mGamma[m];
            row += mGamma[m].size();
        }
    }
    
    // extract comm
    void extractComm(const eigen::DColX &buffer, int &row) {  
        // reference element count
        mElementCount += (int)round(buffer(row, 0));
        
        // number of windows
        int M = (int)round(buffer(row + 1, 0));
        row += 2;
        // gather window information and add new ones if needed
        std::vector<eigen::DMatX2> wins(M);
        for (int m = 0 ; m < M; m++) {
            wins[m] = eigen::DMatX2::Zero(buffer(row + 2, 0), 2);
            wins[m].col(0) = eigen::DColX::LinSpaced(buffer(row + 2, 0), buffer(row + 0, 0), buffer(row + 1, 0));
            row += 3;
        }
        eigen::IColX winTags = addWindows(wins);
        for (int m = 0 ; m < M; m++) {
            // size info
            int sizeMassFluid = (int)round(buffer(row, 0));
            int sizeMassSolid = (int)round(buffer(row + 1, 0));
            int sizeNormalSFA = (int)round(buffer(row + 2, 0));
            int sizeNormalTop = (int)round(buffer(row + 3, 0));
            int sizeGamma = (int)round(buffer(row + 4, 0));
            int countGammasAdded = (int)round(buffer(row + 5, 0));
            row += 6;
            
            // mass
            op1D_3D::addTo(buffer.block(row, 0, sizeMassFluid, 1), mMassFluid[winTags(m)]);
            row += sizeMassFluid;
            op1D_3D::addTo(buffer.block(row, 0, sizeMassSolid, 1), mMassSolid[winTags(m)]);
            row += sizeMassSolid;
            // solid-fluid
            op1D_3D::addTo(Eigen::Map<const eigen::DMatX3>
                           (buffer.block(row, 0, sizeNormalSFA, 1).data(),
                            sizeNormalSFA / 3, 3), mNormalSFA[winTags(m)]);
            row += sizeNormalSFA;
            // ocean load
            op1D_3D::addTo(Eigen::Map<const eigen::DMatX3>
                           (buffer.block(row, 0, sizeNormalTop, 1).data(),
                            sizeNormalTop / 3, 3), mNormalTop[winTags(m)]);
            row += sizeNormalTop;
            // sponge boundary
            op1D_3D::addTo(buffer.block(row, 0, sizeGamma, 1), mGamma[winTags(m)]);
            mCountGammasAdded[winTags(m)] += countGammasAdded;
            row += sizeGamma;
        }
        add3DInfoFromWindowsOnOtherRanks();
    }
    
    
    
    // release to domain
    void release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain);
    
    // get Point after release
    const std::shared_ptr<Point> &getPoint() const {
        return mPoint;
    }
    const std::shared_ptr<FluidPointWindow> &getFluidPointWindow(int m) const {return mFluidPointWindows[m];};
    const std::shared_ptr<SolidPointWindow> &getSolidPointWindow(int m) const {return mSolidPointWindows[m];};
    
    // get global tag
    int getGlobalTag() const {
        return mGlobalTag;
    }
    
    // get reference element count
    int getElementCount() const {
        return mElementCount;
    }
    
private:
    void computeWindowSumSampling(eigen::DColX &windowSumPhi, std::vector<int> &aligned) const {
        if (mWindows.size() == 1 && mWindowsFromOtherRanks.size() == 0) {
            aligned.push_back(true);
            windowSumPhi = mWindows[0].col(0);
            return;
        } 
      
        double min_dphi = std::numeric_limits<double>::max();
        double min_phi1 = std::numeric_limits<double>::max();
        for (auto &win: mWindows) {
            double dphi = window_tools::setBounds2Pi(win(1, 0) - win(0, 0));
            if (dphi < min_dphi) {
                min_dphi = dphi;
                min_phi1 = win(0, 0);
            } else if (dphi == min_dphi && win(0, 0) < min_phi1) {
                min_phi1 = win(0, 0);
            }
        }
        
        bool not_aligned = false;
        for (auto &win: mWindowsFromOtherRanks) {
            if (win(1) < min_dphi) {
                min_dphi = win(1);
                min_phi1 = win(0);
                not_aligned = true;
            } else if (win(1) == min_dphi && win(0, 0) < min_phi1) {
                min_phi1 = win(0);
                not_aligned = true;
            }
        }
        
        int nr = (int)ceil((2 * numerical::dPi - mAngleTol) / min_dphi);
        double dphi_sum = 2 * numerical::dPi / nr;
        
        if (std::abs(dphi_sum - min_dphi) > mAngleTol) not_aligned = true;
        
        windowSumPhi = eigen::DColX::LinSpaced(nr, min_phi1, min_phi1 + 2 * numerical::dPi - dphi_sum);
        window_tools::wrapPhi(windowSumPhi);
        
        aligned.resize(mWindows.size(), -1);
        if (not_aligned) return;
        
        for (int m; m < mWindows.size(); m++) {
            double dphi = window_tools::setBounds2Pi(mWindows[m](1, 0) - mWindows[m](0, 0));
            if (std::abs(dphi - dphi_sum) > mAngleTol) continue;
            for (int i = 0; i < windowSumPhi.rows(); i++) {
                if (std::abs(windowSumPhi(i) - mWindows[m](0, 0)) < mAngleTol) {
                    aligned[m] = i;
                    break;
                }
            }
        }
    };
    
    ///////////////////////////// data /////////////////////////////
private:
    // tag
    int mGlobalTag = -1;
    
    // nr
    std::vector<eigen::DMatX2> mWindows; // phi + winFunc
    std::vector<eigen::DCol2> mWindowsFromOtherRanks; // phi1 + dphi (only used to determine finest azimuthal sampling for window sum)
    double mAngleTol;
    
    // coords
    eigen::DCol2 mCoords = eigen::DCol2::Zero(2, 1);
    
    // reference element count
    int mElementCount = 0;
    
    std::vector<std::vector<std::pair<winmap, winmap>>> mMapWindows; // for interpolation between different window sets
    int mNumNewWins = 0;
    
    // mass
    std::vector<eigen::DColX> mMassFluid;
    std::vector<eigen::DColX> mMassSolid;
    
    // solid-fuild
    // unassembled
    std::vector<eigen::DMatX3> mNormalSFU;
    // assembled
    std::vector<eigen::DMatX3> mNormalSFA;
    
    // Clayton ABC {key, {fluid/solid, normal, rho * vp, rho * vs}}
    std::vector<std::map<std::string, std::vector<
    std::tuple<bool, eigen::DMatX3, eigen::DColX, eigen::DColX>>>> mClaytonABC;
    
    // Kosloff_Kosloff sponge boundary
    std::vector<eigen::DColX> mGamma;
    std::vector<int> mCountGammasAdded;
    
    // ocean load
    std::vector<eigen::DMatX3> mNormalTop;
    std::vector<eigen::DColX> mSumRhoDepth;
    
    // axial
    bool mAxial = false;
    
    // surface
    bool mSurface = false;
    
    // pointers after release
    std::shared_ptr<Point> mPoint = nullptr;
    std::vector<std::shared_ptr<SolidPointWindow>> mSolidPointWindows;
    std::vector<std::shared_ptr<FluidPointWindow>> mFluidPointWindows;
};

#endif /* GLLPoint_hpp */
