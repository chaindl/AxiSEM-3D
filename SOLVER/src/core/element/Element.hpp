//
//  FluidElement.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#ifndef Element_hpp
#define Element_hpp

// point
#include <array>
class Point;
class CoordTransform;

// output
#include "channel.hpp"
#include "vector_tools.hpp"
#include "ElementWindow.hpp"
#include "FluidElementWindow.hpp"
#include "SolidElementWindow.hpp"

class Element {
public:
    // constructor
    Element(const int &quadTag, std::vector<std::unique_ptr<ElementWindow>> &windows, 
            const std::array<std::shared_ptr<Point>, spectral::nPEM> &points,
            const std::shared_ptr<WinInt> &interpolator,
            const eigen::IMatX2 &interpTags);
    
    // copy constructor
    Element(const Element &other);
    
private:

    void expandWorkspace(const int nr) {
        if (sWin1ToWin2.rows() < nr) {
            sWin1ToWin2.resize(nr, spectral::nPEM * 6);
            sWin2ToWin1.resize(nr, spectral::nPEM * 6);
        }
        
        if (sStrainWindows.size() < mWindows.size()) sStrainWindows.resize(mWindows.size());
        for (int m = 0; m < mWindows.size(); m++) {
            if (sStrainWindows[m].rows() < mWindows[m]->getNr()) sStrainWindows[m].resize(mWindows[m]->getNr(), spectral::nPEM * 6);
        }
    };
    
    std::vector<std::unique_ptr<ElementWindow>> copyWindows() const {
        typedef FluidElementWindow<FluidRCPointWindow<1, numerical::Real>> rfw;
        typedef FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>> cfw;
        typedef SolidElementWindow<SolidRCPointWindow<3, numerical::Real>> rsw;
        typedef SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>> csw;
        
        std::vector<std::unique_ptr<ElementWindow>> newWins;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mWindows[m]->isFluid()) {
                if (mWindows.size() > 1) {
                    newWins.push_back(std::make_unique<rfw>(*dynamic_cast<rfw*>(mWindows[m].get())));
                } else {
                    newWins.push_back(std::make_unique<cfw>(*dynamic_cast<cfw*>(mWindows[m].get())));
                }
            } else {
                if (mWindows.size() > 1) {
                    newWins.push_back(std::make_unique<rsw>(*dynamic_cast<rsw*>(mWindows[m].get())));
                } else {
                    newWins.push_back(std::make_unique<csw>(*dynamic_cast<csw*>(mWindows[m].get())));
                }
            }
        }
        return newWins;
    }
    
public:
    void countInfo(std::map<std::string, int> &typeCountElementWindow) const {
        for (auto &win: mWindows) {
            vector_tools::aggregate(typeCountElementWindow, win->typeInfo(), 1);
        }
    }
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << "EleWins:";
        for (int m = 0; m < mWindows.size(); m++) {
            ss << mWindows[m]->costSignature();
            if (m < mWindows.size() - 1) ss << "$";
        }
        return ss.str();
    }
    void createCoordTransform();
    
    const std::shared_ptr<const CoordTransform> getCoordTransform() { // called by moment source
        if (!mTransform) createCoordTransform();
        return mTransform;
    }
    
    void setDomainTag (int tag) {
        mDomainTag = tag;
    }
    int getDomainTag() const {
        return mDomainTag;
    }
    
    /////////////////////////// point ///////////////////////////
    // get point
    Point &getPoint(int ipnt) const;
    const eigen::DRow2 &getPointCoords(int ipnt) const;
    
    int getWindowNr(int m) const {return mWindows[m]->getNr();};
    int getWindowNu_1(int m) const {return mWindows[m]->getNu_1_buffered();};
    int getWindowNu_1_noBuffer(int m) const {return mWindows[m]->getNu_1();};
    int getMaxNu_1() const {
        int max_nu_1 = -1;
        for (auto &win: mWindows) {
            max_nu_1 = std::max({max_nu_1, win->getNu_1_buffered()});
        }
        return max_nu_1;
    };
    int getMaxNu_1_noBuffer() const {
        int max_nu_1 = -1;
        for (auto &win: mWindows) {
            max_nu_1 = std::max({max_nu_1, win->getNu_1()});
        }
        return max_nu_1;
    };
    
    int getQuadTag() const {return mQuadTag;};
    
    std::vector<int> findBoundaryPointsByTag(const std::vector<int> &boundaryMeshTags) const;
    std::vector<int> findBoundaryPointsByCrds(const std::vector<double> &boundaryCrdsRorZ,
                         const std::vector<double> &boundaryCrdsTorS, double distTol) const;
    
    /////////////////////////// time loop ///////////////////////////
    // collect displacement from points
    void
    collectDisplFromPoints(eigen::vec_ar1_CMatPP_RM &displElem) const;
    
    // displacement to stiffness
    void displToStiff(const eigen::vec_ar1_CMatPP_RM &displElem,
                      eigen::vec_ar1_CMatPP_RM &stiffElem) const;
    
    // add stiffness to points
    // allow a derived class to change stiffElem (no const)
    void
    addStiffToPoints(eigen::vec_ar1_CMatPP_RM &stiffElem) const;
    
    // compute stiffness term
    void computeStiff() const;
    
    // measure cost
    double measure(int count) const;
    
    void overlapAndAddStrain() const;
    
    /////////////////////////// source ///////////////////////////
    void preparePressureSource(int m) const {
        mWindows[m]->preparePressureSource();
    }
    void addPressureSource(int m, const eigen::CMatXN &pressure, int nu_1_pressure) const {
        mWindows[m]->addPressureSource(pressure, nu_1_pressure);
    }
    
    void prepareForceSource(int m) const {
        mWindows[m]->prepareForceSource();
    }
    void addForceSource(int m, const eigen::CMatXN3 &force, int nu_1_force) const {
        mWindows[m]->addForceSource(force, nu_1_force);
    }
                        
    void prepareMomentSource(int m) const {
        mWindows[m]->prepareMomentSource();
    }
    void addMomentSource(int m, const eigen::CMatXN6 &moment, int nu_1_moment) const {
        mWindows[m]->addMomentSource(moment, nu_1_moment, mTransform);
    }
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare wavefield output
    void prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops, int winTag,
                                bool enforceCoordTransform) {
        bool needTransform = mWindows[winTag]->prepareWavefieldOutput(chops);
        if (enforceCoordTransform && needTransform) {
            createCoordTransform();
        }                        
    }
    void prepareWavefieldOutput(const channel::solid::ChannelOptions &chops, int winTag,
                                bool enforceCoordTransform) {
        bool needTransform = mWindows[winTag]->prepareWavefieldOutput(chops);
        if (enforceCoordTransform && needTransform) {
            createCoordTransform();
        }                           
    }
    
    bool getMajorityDisplInRTZ(const std::vector<std::tuple<int, double, double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[std::get<0>(winPhi)]->displInRTZ()) {
                nInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            } else {
                nNotInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityNablaInRTZ(const std::vector<std::tuple<int, double, double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[std::get<0>(winPhi)]->nablaInRTZ()) {
                nInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            } else {
                nNotInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityStrainInRTZ(const std::vector<std::tuple<int, double, double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[std::get<0>(winPhi)]->strainInRTZ()) {
                nInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            } else {
                nNotInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityCurlInRTZ(const std::vector<std::tuple<int, double, double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[std::get<0>(winPhi)]->curlInRTZ()) {
                nInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            } else {
                nNotInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityStressInRTZ(const std::vector<std::tuple<int, double, double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[std::get<0>(winPhi)]->stressInRTZ()) {
                nInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            } else {
                nNotInRTZ += mWindows[std::get<0>(winPhi)]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    //// Solid and Fluid         
    // displ field
    void getDisplField(eigen::CMatXN3 &displ, bool needRTZ, int m) const;
    
    //// Fluid only
    // chi field
    void getChiField(eigen::CMatXN &chi, int m) const;
    // pressure field
    void getPressureField(eigen::CMatXN &pressure, int m) const;
    // delta field
    void getDeltaField(eigen::CMatXN &delta, int m) const;
    
    //// Solid only
    // nabla field
    void getNablaField(eigen::CMatXN9 &nabla, bool needRTZ, int m, int &nu) const;
    // strain field
    void getStrainField(eigen::CMatXN6 &strain, bool needRTZ, int m, int &nu) const;
    // curl field
    void getCurlField(eigen::CMatXN3 &curl, bool needRTZ, int m, int &nu) const;
    // stress field
    void getStressField(eigen::CMatXN6 &stress, bool needRTZ, int m, int &nu) const;
    
private:
    // points
    const int mQuadTag; 
    int mDomainTag;
    std::vector<std::unique_ptr<ElementWindow>> mWindows;
    std::array<std::shared_ptr<Point>, spectral::nPEM> mPoints;
    std::shared_ptr<const CoordTransform> mTransform;
    
    const std::shared_ptr<WinInt> mInterpolator;
    eigen::IMatX2 mInterpTags;
    inline static eigen::RMatXN6 sWin1ToWin2 =
    eigen::RMatXN6(0, spectral::nPEM * 6);;
    inline static eigen::RMatXN6 sWin2ToWin1 =
    eigen::RMatXN6(0, spectral::nPEM * 6);;
    inline static std::vector<eigen::RMatXN6> sStrainWindows = std::vector<eigen::RMatXN6>(0);
};

#endif /* FluidElement_hpp */
