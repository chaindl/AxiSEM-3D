//
//  FluidPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid GLL point
#ifndef Point_hpp
#define Point_hpp

#include "point_time.hpp"
#include "PointWindow.hpp"
#include "WindowSum.hpp"
#include "vector_tools.hpp"
#include <memory>

class Point {
public:
    // constructor
    Point(int meshTag, eigen::DRow2 crds):
    mMeshTag(meshTag), mCoords(crds) {
        // nothing
    };
    
    // destructor
    virtual ~Point() = default;

    void addWindow(std::shared_ptr<PointWindow> pw) {
        mWindows.push_back(pw);
    };
    
    void addWindowSum(std::shared_ptr<WindowSum> pws) {
        mWindowSums.push_back(pws);
        pws->setInFourier();
    };
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << "PointWins:";
        for (int m = 0; m < mWindows.size(); m++) {
            ss << mWindows[m]->costSignature();
            if (m < mWindows.size() - 1) ss << "$";
        }
        return ss.str();
    }
    
    /////////////////////////// time loop ///////////////////////////
    void combineWindows() {
        for (auto win: mWindows) {
            win->transformToPhysical();
        }
        for (auto ws: mWindowSums) {
            ws->overlapAndAddStiff();
        }
    };

    void separateWindows() {
        for (auto ws: mWindowSums) {
            ws->scatterStiffToWindows();
        }
        
    };

    // stiff to accel
    void computeStiffToAccel() {
        for (auto win: mWindows) {
            win->computeStiffToAccel();
            win->transformToFourier();
            win->applyPressureSource();
        }
    };

    void countInfo(std::map<std::string, int> &typeCountPointWindow) {
        for (auto win: mWindows) {
            vector_tools::aggregate(typeCountPointWindow, win->typeInfo(), 1);
        }
    };
    
    // location
    const eigen::DRow2 &getCoords() const {
        return mCoords;
    }
    
    // tag
    int getMeshTag() const {
        return mMeshTag;
    }
    
    void randomDispl() {
        for (auto win: mWindows) {
            win->randomDispl();
        }
    };
    void randomStiff() {
        for (auto win: mWindows) {
            win->randomStiff();
        }
    };
    void resetToZero() {
        for (auto win: mWindows) {
            win->resetToZero();
        }
    };
    
    bool stable() const {
        for (auto win: mWindows) {
            if (!win->stable()) return false;
        }
        return true;
    }
    
    /////////////////////////// domain ///////////////////////////
    // set domain tag
    void setDomainTag(int domainTag) {
        mDomainTag = domainTag;
    }
    
    // get domain tag
    int getDomainTag() const {
        return mDomainTag;
    }
    
    const std::vector<std::shared_ptr<PointWindow>> &getWindows() const {return mWindows;};
    const std::vector<std::shared_ptr<WindowSum>> &getWindowSums() const {return mWindowSums;};

    /////////////////////// wavefield scanning ///////////////////////
    // enable scanning
    // need this abstract for vertex-only option
    void enableScanning() {
        if (mWindowSums[0]->onlyOneWindow()) {
            mWindows[0]->enableScanning(); // this works because in solid-fluid windows the first window is solid
        } else {
            throw std::runtime_error("Point::enableScanning || Wavefield scanning not implemented for multiple windows.");
        }
    };
    
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) {
        mWindows[0]->doScanning(relTolFourierH2, relTolH2, absTolH2, maxNumPeaks);        
    };
    
    int getStartingNr() const {
        return mWindows[0]->getNr();
    };
    
    int reportScanningNr() const {
        return mWindows[0]->reportScanningNr();
    }
                    
private:
   std::vector<std::shared_ptr<PointWindow>> mWindows;
   std::vector<std::shared_ptr<WindowSum>> mWindowSums;
   int mMeshTag;
   int mDomainTag;
   eigen::DRow2 mCoords;
};

#endif
