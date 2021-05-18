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
#include "SFRCPointWindow.hpp"
#include "WindowInterpolator.hpp"
#include <map>
#include <vector>
#include <memory>
#include <iostream>
// release
class ABC;
class TimeScheme;
class Domain;
class Point;

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
            mAngleTol = (mCoords(0) > numerical::dEpsilon) ? distTol / mCoords(0) : numerical::dPi / 5;
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
    eigen::IColX addWindows(const std::vector<eigen::DMatX2> &wins, const bool globalWin) {      
        eigen::IColX tags = eigen::IColX::Constant(wins.size(), -1);
        
        // first set of windows added
        if (mWindows.size() == 0) {
            for (int mn = 0; mn < wins.size(); mn++) {
                tags(mn) = mWindows.size();
                mWindows.push_back(wins[mn]);
                mGlobalWin.push_back(globalWin);
            }
            mNumNewWins = mWindows.size();
            allocateNewWindows();
            mNumNewWins = 0;
            return tags;
        }

        // pre-existing windows
        std::vector<bool> found(mWindows.size(), false);
        for (int mo = 0; mo < mWindows.size(); mo++) {
            for (int mn = 0; mn < wins.size(); mn++) {
                if (std::abs(mWindows[mo](0, 0) - wins[mn](0, 0)) < mAngleTol) {
                    if ((mGlobalWin[mo] && globalWin)          // case 1 : both global windows (no need for equal nr)
                        || (wins[mn].rows() == mWindows[mo].rows() // case 2 : local windows with equal nr and dphi
                            && std::abs(mWindows[mo](1, 0) - wins[mn](1, 0)) < mAngleTol)) {
                        found[mo] = true;
                        tags(mn) = mo;
                        break;
                    }
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
            tags(mn) = mWindows.size();
            mWindows.push_back(wins[mn]);
            mGlobalWin.push_back(globalWin);
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
        mGlobalWin.resize(mWindows.size() - mNumNewWins);
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
        for (int m = 0; m < mWindows.size(); m++) {
            buffer(row + 0, 0) = mWindows[m](0, 0);
            buffer(row + 1, 0) = mWindows[m](mWindows[m].rows() - 1, 0);
            buffer(row + 2, 0) = mWindows[m].rows();
            buffer(row + 3, 0) = mGlobalWin[m];
            row += 4;
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
        bool globalWin;
        for (int m = 0 ; m < M; m++) {
            wins[m] = eigen::DMatX2::Zero(buffer(row + 2, 0), 2);
            wins[m].col(0) = eigen::DColX::LinSpaced(buffer(row + 2, 0), buffer(row + 0, 0), buffer(row + 1, 0));
            globalWin = (bool)round(buffer(row + 3, 0));
            row += 4;
        }
        eigen::IColX winTags = addWindows(wins, globalWin);
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
    void release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain, 
         const std::shared_ptr<WindowInterpolator<numerical::Real>> &interpolator);
    
    // get Point after release
    const std::shared_ptr<Point> &getPoint() const {
        return mPoint;
    }
    
    const std::shared_ptr<FluidRCPointWindow<1, numerical::ComplexR>> getFluidPointWindowC(int m) const {
        return std::dynamic_pointer_cast<FluidRCPointWindow<1, numerical::ComplexR>>(mFluidPointWindows[m]);
    };
    const std::shared_ptr<SolidRCPointWindow<3, numerical::ComplexR>> getSolidPointWindowC(int m) const {
        return std::dynamic_pointer_cast<SolidRCPointWindow<3, numerical::ComplexR>>(mSolidPointWindows[m]);
    };
    const std::shared_ptr<FluidRCPointWindow<1, numerical::Real>> getFluidPointWindowR(int m) const {
        return std::dynamic_pointer_cast<FluidRCPointWindow<1, numerical::Real>>(mFluidPointWindows[m]);
    };
    const std::shared_ptr<SolidRCPointWindow<3, numerical::Real>> getSolidPointWindowR(int m) const {
        return std::dynamic_pointer_cast<SolidRCPointWindow<3, numerical::Real>>(mSolidPointWindows[m]);
    };
    
    // get global tag
    int getGlobalTag() const {
        return mGlobalTag;
    }
    
    // get reference element count
    int getElementCount() const {
        return mElementCount;
    }
    
private:
    void computeWindowSumSampling(eigen::DColX &knotsWhole, std::vector<eigen::DColX> &relPhis, 
                                  std::vector<eigen::IColX> &posIndices, int &nr) {
        // no windows sum required except for messaging
        if (mWindows.size() == 1 && mWindowsFromOtherRanks.size() == 0) {
            nr = mWindows[0].rows();
            return;
        } 
        
        // find minimal sampling
        double min_dphi = window_tools::setBounds2Pi(mWindows[0](1, 0) - mWindows[0](0, 0));
        double min_phi1 = mWindows[0](0, 0);
        for (int m = 1; m < mWindows.size(); m++) {
            double dphi = window_tools::setBounds2Pi(mWindows[m](1, 0) - mWindows[m](0, 0));
            if (dphi < min_dphi - mAngleTol) {
                min_dphi = dphi;
                min_phi1 = mWindows[m](0, 0);
            } else if (abs(dphi - min_dphi) < mAngleTol && mWindows[m](0, 0) < min_phi1) {
                min_phi1 = mWindows[m](0, 0);
            }
        }
        
        // check if there is smaller sampling on other ranks
        for (auto &win: mWindowsFromOtherRanks) {
            if (win(1) < min_dphi) {
                min_dphi = win(1);
                min_phi1 = win(0);
            } else if (abs(win(1) - min_dphi) < mAngleTol && win(0) < min_phi1) {
                min_phi1 = win(0);
            }
        }
        
        // set physical sampling for window sum
        // (this must be consistent between ranks
        // for messaging)
        nr = (int)ceil((2 * numerical::dPi - mAngleTol) / min_dphi);
        double dphi_sum = 2 * numerical::dPi / nr;
        eigen::DColX windowSumPhi = eigen::DColX::LinSpaced(nr, min_phi1, min_phi1 + 2 * numerical::dPi - dphi_sum);
        window_tools::wrapPhi(windowSumPhi);
        
        // process windows which are aligned (i.e. azimuthal sampling points 
        // are co-located with azimuthal sampling points from window sum)
        std::vector<bool> aligned(mWindows.size());
        for (int m = 0; m < mWindows.size(); m++) {
            double dphi = window_tools::setBounds2Pi(mWindows[m](1, 0) - mWindows[m](0, 0));
            aligned[m] = false;
            if (std::abs(dphi - dphi_sum) <= mAngleTol) {
                for (int i = 0; i < windowSumPhi.rows(); i++) {
                    if (std::abs(windowSumPhi(i) - mWindows[m](0, 0)) < mAngleTol) {
                        posIndices[m] = eigen::IColX::LinSpaced(mWindows[m].rows(), i, i + mWindows[m].rows() - 1);
                        posIndices[m] = posIndices[m].unaryExpr([windowSumPhi](const int idx) { return idx%((int)windowSumPhi.rows()); });
                        aligned[m] = true;
                        break;
                    }
                }
            }
        }
        
        // no interpolation required
        if (std::all_of(aligned.begin(), aligned.end(), [](bool a){ return a; })) return;
        
        // calculate knots for spline interpolation
        // (extend on both sides to create periodic BC)
        int nr_ext = nr + 2 * interpolation::NrExtend;
        knotsWhole = eigen::DColX::LinSpaced(nr_ext, 0, 1);
        
        // process windows which are not aligned 
        window_tools::unwrapPhi(windowSumPhi);
        for (int m = 0; m < mWindows.size(); m++) {
            if (!aligned[m]) {
                double phi_start = mWindows[m](0, 0);
                double phi_end = mWindows[m](mWindows[m].rows() - 1, 0);
                if (phi_start < windowSumPhi(0)) phi_start += 2 * numerical::dPi;
                if (phi_end < windowSumPhi(0)) phi_end += 2 * numerical::dPi;
                int i1 = 0, i2 = 2;
                while (i1 < windowSumPhi.rows() && windowSumPhi(i1) < phi_start) i1++;
                while (i2 < windowSumPhi.rows() && windowSumPhi(i2) < phi_end) i2++;
                i2 -= 1;
                if (i2 < i1) i2 += windowSumPhi.rows();
                posIndices[m] = eigen::IColX::LinSpaced(i2 - i1 + 1, i1, i2);
                posIndices[m] = posIndices[m].unaryExpr([windowSumPhi](const int idx) { return idx%((int)windowSumPhi.rows()); });
                relPhis[m] = windowSumPhi(i1%windowSumPhi.rows())
                           + eigen::DColX::LinSpaced(posIndices[m].rows(), 0, (posIndices[m].rows() - 1)).array() * dphi_sum;
                double windowSize = phi_end - phi_start;
                if (windowSize < 0) windowSize += 2 * numerical::dPi;
                // relative phis needed for interpolation from window to window sum
                relPhis[m] = (relPhis[m].array() - phi_start) / windowSize;
                
                // relative phis needed for interpolation from window sum to window
                mWindows[m].col(0) = (mWindows[m].col(0).array() + interpolation::NrExtend * dphi_sum) / ((nr_ext - 1) * dphi_sum);
            }
        }
    };
    
    ///////////////////////////// data /////////////////////////////
private:
    // tag
    int mGlobalTag = -1;
    
    // nr
    std::vector<eigen::DMatX2> mWindows; // phi + winFunc
    std::vector<bool> mGlobalWin;
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
    std::vector<std::shared_ptr<PointWindow>> mSolidPointWindows;
    std::vector<std::shared_ptr<PointWindow>> mFluidPointWindows;
};

#endif /* GLLPoint_hpp */
