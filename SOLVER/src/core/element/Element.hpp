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
            const std::array<std::shared_ptr<Point>, spectral::nPEM> &points);
    
    // copy constructor
    Element(const Element &other);
    
private:

    void expandWorkspace(const int nr) {
        sWin1ToWin2.resize(nr, 1);
        sWin2ToWin1.resize(nr, 1);
    };
    
    std::vector<std::unique_ptr<ElementWindow>> copyWindows() const {
        std::vector<std::unique_ptr<ElementWindow>> newWins;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mWindows[m]->isFluid()) {
                newWins.push_back(std::make_unique<FluidElementWindow>(mWindows[m]->getFluidElementWindow()));
            } else {
                newWins.push_back(std::make_unique<SolidElementWindow>(mWindows[m]->getSolidElementWindow()));
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
    
    FluidElementWindow &getFluidElementWindow(int m) const {
        mWindows[m]->getFluidElementWindow();
    };
    SolidElementWindow &getSolidElementWindow(int m) const {
        mWindows[m]->getSolidElementWindow();
    };
    
    int getWindowNr(int m) const {return mWindows[m]->getNr();};
    int getWindowNu_1(int m) const {return mWindows[m]->getNu_1();};
    int getMaxNu_1() const {
        int max_nu_1 = -1;
        for (auto &win: mWindows) {
            max_nu_1 = std::max({max_nu_1, win->getNu_1()});
        }
        return max_nu_1;
    };
    
    int getQuadTag() const {return mQuadTag;};
    
    void setAlignment(const double tol) {
        if (mWindows.size() == 1) {
            mWindows[0]->setAlignment(true);
            return;
        }
        
        for (int m = 0; m < mWindows.size(); m++) {
            int mnext = (m == mWindows.size() - 1) ? 0 : m + 1;
            eigen::RColX phi_next = mWindows[mnext]->getPhiForInterp(1);
            eigen::RColX phi = mWindows[m]->getPhiForInterp(0);
            if (phi_next.size() == phi.size()) {
                if ((phi_next - phi).norm() < tol) {
                    mWindows[m]->setAlignment(true);
                }
            }
        }
    }
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
    void interpolate(const eigen::RColX &phi_q, eigen::RColX &fun, const eigen::RColX &phi, const int side) const;
    
    /////////////////////////// source ///////////////////////////
    void preparePressureSource(int m) const {
        getFluidElementWindow(m).preparePressureSource();
    }
    void addPressureSource(int m, const eigen::CMatXN &pressure, int nu_1_pressure) const {
        getFluidElementWindow(m).addPressureSource(pressure, nu_1_pressure);
    }
    
    void prepareForceSource(int m) const {
        getSolidElementWindow(m).prepareForceSource();
    }
    void addForceSource(int m, const eigen::CMatXN3 &force, int nu_1_force) const {
        getSolidElementWindow(m).addForceSource(force, nu_1_force);
    }
                        
    void prepareMomentSource(int m) const {
        getSolidElementWindow(m).prepareMomentSource();
    }
    void addMomentSource(int m, const eigen::CMatXN6 &moment, int nu_1_moment) const {
        getSolidElementWindow(m).addMomentSource(moment, nu_1_moment, mTransform);
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
                                  
    }
    
    bool getMajorityDisplInRTZ(const std::vector<std::pair<int,double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[winPhi.first]->displInRTZ()) {
                nInRTZ += mWindows[winPhi.first]->getNu_1();
            } else {
                nNotInRTZ += mWindows[winPhi.first]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityNablaInRTZ(const std::vector<std::pair<int,double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[winPhi.first]->nablaInRTZ()) {
                nInRTZ += mWindows[winPhi.first]->getNu_1();
            } else {
                nNotInRTZ += mWindows[winPhi.first]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityStrainInRTZ(const std::vector<std::pair<int,double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[winPhi.first]->strainInRTZ()) {
                nInRTZ += mWindows[winPhi.first]->getNu_1();
            } else {
                nNotInRTZ += mWindows[winPhi.first]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityCurlInRTZ(const std::vector<std::pair<int,double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[winPhi.first]->curlInRTZ()) {
                nInRTZ += mWindows[winPhi.first]->getNu_1();
            } else {
                nNotInRTZ += mWindows[winPhi.first]->getNu_1();
            }
        }
        return (nInRTZ > nNotInRTZ);
    }
    
    bool getMajorityStressInRTZ(const std::vector<std::pair<int,double>> &windowPhis) const {
        int nInRTZ = 0, nNotInRTZ = 0;
        for (auto &winPhi: windowPhis) {
            if (mWindows[winPhi.first]->stressInRTZ()) {
                nInRTZ += mWindows[winPhi.first]->getNu_1();
            } else {
                nNotInRTZ += mWindows[winPhi.first]->getNu_1();
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
    void getNablaField(eigen::CMatXN9 &nabla, bool needRTZ, int m) const;
    // strain field
    void getStrainField(eigen::CMatXN6 &strain, bool needRTZ, int m) const;
    // curl field
    void getCurlField(eigen::CMatXN3 &curl, bool needRTZ, int m) const;
    // stress field
    void getStressField(eigen::CMatXN6 &stress, bool needRTZ, int m) const;
    
private:
    // points
    const int mQuadTag; 
    int mDomainTag;
    std::vector<std::unique_ptr<ElementWindow>> mWindows;
    std::array<std::shared_ptr<Point>, spectral::nPEM> mPoints;
    std::shared_ptr<const CoordTransform> mTransform;
    
    inline static eigen::RColX sWin1ToWin2;
    inline static eigen::RColX sWin2ToWin1;
};

#endif /* FluidElement_hpp */
