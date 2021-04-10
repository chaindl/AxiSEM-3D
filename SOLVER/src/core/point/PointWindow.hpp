//
//  Point.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  GLL point

#ifndef PointWindow_hpp
#define PointWindow_hpp

#include "eigen_point.hpp"
#include "Mass.hpp"
#include "bstring.hpp"
#include <iostream>
class Point;
class SolidPointWindow;
class FluidPointWindow;

template <int dims>
struct Fields {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, dims> RMat;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, dims> CMat;
    
    CMat mStiff = CMat(0, dims);
    RMat mStiffR = RMat(0, dims);
    CMat mDispl = CMat(0, dims);
    CMat mVeloc = CMat(0, dims);
    CMat mAccel = CMat(0, dims);
};

class PointWindow {
public:
    // constructor
    PointWindow(const eigen::RMatX2 &windowSumPhi,
          std::unique_ptr<const Mass> &mass,
          const std::shared_ptr<const Point> point):
    mNr(windowSumPhi.rows()), mNu_1(mNr / 2 + 1), mWindowSumPhi(windowSumPhi.col(0)), 
    mWindowSumFrac(windowSumPhi.col(1)), mMass(mass.release()), mPoint(point) {
        // nothing
    }
    
    // destructor
    virtual ~PointWindow() = default;
    
    virtual bool isFluid() const = 0;
    
    // type info
    std::string typeInfo() const {
        return bstring::typeName(*this) + "$" + bstring::typeName(*mMass);
    }
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << typeInfo() << "$" << mNr;
        return ss.str();
    }
    
    void setInFourier(bool inFourier) {
        mInFourier = inFourier;
    }
    
    bool is3D() {
        return mMass->is3D();
    }
    
    virtual Fields<1> &getFluidFields() {
        throw std::runtime_error("PointWindow || requesting fluid fields from solid window.");
    }
    
    virtual Fields<3> &getSolidFields() {
        throw std::runtime_error("PointWindow || requesting solid fields from fluid window.");
    }
    
    virtual const Fields<1> &getFluidFields() const {
        throw std::runtime_error("PointWindow || requesting fluid fields from solid window.");
    }
    
    virtual const Fields<3> &getSolidFields() const {
        throw std::runtime_error("PointWindow || requesting solid fields from fluid window.");
    }
    /////////////////////////// time loop ///////////////////////////
    
    virtual void randomDispl() = 0;
    virtual void randomStiff() = 0;
    virtual void resetToZero() = 0;
    
    // check stability
    virtual bool stable() const = 0;
    
    // stiff to accel
    virtual void transformToPhysical() = 0;
    virtual void computeStiffToAccel() = 0;
    virtual void transformToFourier() = 0;
    virtual void applyPressureSource() {
        //nothing
    };
    virtual void maskNyquist() = 0;
    virtual void applyAxialBC() = 0;
    /////////////////////////// properties ///////////////////////////
    // nr
    int getNr() const {
        return mNr;
    }
    
    // nu
    int getNu_1() const {
        return mNu_1;
    }
    
    const eigen::RColX &getPhiForWindowSum() const {
        return mWindowSumPhi;
    }
    
    int getMeshTag() const;
    
    /////////////////////// wavefield scanning ///////////////////////
    // enable scanning
    // need this abstract for vertex-only option
    virtual void enableScanning() = 0;
    virtual void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                            numerical::Real absTolH2, int maxNumPeaks) = 0;
    virtual int reportScanningNr() const = 0;
    
protected:
    bool mInFourier = false;
    
    // order
    const int mNr;
    const int mNu_1;
    
    // mass
    const std::unique_ptr<const Mass> mMass;
    const std::weak_ptr<const Point> mPoint;
    const eigen::RColX mWindowSumPhi;
    const eigen::RColX mWindowSumFrac;
};

#endif /* Point_hpp */
