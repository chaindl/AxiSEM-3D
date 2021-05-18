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
#include "eigen_element.hpp"
#include "Mass.hpp"
#include "bstring.hpp"
#include <iostream>

class Point;
class TimeScheme;

template <int dims, typename DisplType>
struct Fields {};

template <int dims>
struct Fields<dims, numerical::ComplexR> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, dims> RMat;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, dims> CMat;
    
    CMat mStiff = CMat(0, dims);
    RMat mStiffR = RMat(0, dims);
    CMat mDispl = CMat(0, dims);
    CMat mVeloc = CMat(0, dims);
    CMat mAccel = CMat(0, dims);
    
    CMat& mStiffUpdate = mStiff;
};

template <int dims>
struct Fields<dims, numerical::Real> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, dims> RMat;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, dims> CMat;
    
    CMat mStiff = CMat(0, dims);
    RMat mStiffR = RMat(0, dims);
    RMat mDispl = RMat(0, dims);
    RMat mVeloc = RMat(0, dims);
    RMat mAccel = RMat(0, dims);
    
    RMat& mStiffUpdate = mStiffR;
};

class PointWindow {
public:
    // constructor
    PointWindow(std::unique_ptr<const Mass> &mass, const eigen::RMatX2 &windowSumPhi, const int tag):
    mNr(windowSumPhi.rows()), mNu_1(mNr / 2 + 1), mMass(mass.release()), 
    mWindowSumPhi(windowSumPhi.col(0)), mWindowSumFrac(windowSumPhi.col(1)), mPointMeshTag(tag) {
        // nothing
    }
    
    // destructor
    virtual ~PointWindow() = default;
    
    virtual void checkCompatibility(const TimeScheme &timeScheme) = 0;
    /////////////////////////// properties ///////////////////////////
    
    int getNr() const {
      return mNr;
    }
    int getNu_1() const {return mNu_1;}
    int getMeshTag() const {return mPointMeshTag;};
    bool is3D() {return mMass->is3D();}
    virtual bool isFluid() const = 0;
    virtual bool storesFieldsInFourier() const = 0;
    
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
    
    virtual void setInFourier(bool inFourier) {};

    virtual Fields<1, numerical::Real> &getFluidFieldsR() = 0;
    virtual Fields<3, numerical::Real> &getSolidFieldsR() = 0;
    virtual const Fields<1, numerical::Real> &getFluidFieldsR() const = 0;
    virtual const Fields<3, numerical::Real> &getSolidFieldsR() const = 0;
    virtual Fields<1, numerical::ComplexR> &getFluidFieldsC() = 0;
    virtual Fields<3, numerical::ComplexR> &getSolidFieldsC() = 0;
    virtual const Fields<1, numerical::ComplexR> &getFluidFieldsC() const = 0;
    virtual const Fields<3, numerical::ComplexR> &getSolidFieldsC() const = 0;
    virtual eigen::RMatX3 &getSolidStiffForWindowSum() = 0;
    virtual eigen::RColX &getFluidStiffForWindowSum() = 0;
    
    virtual void collectStiffFromWindowSum(const eigen::CMatX3 &stiff) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding solid stiff to fluid window.");
    };
        
    virtual void collectStiffFromWindowSum(const eigen::RMatX3 &stiff) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding solid stiff to fluid window.");
    };
    
    virtual void collectStiffFromWindowSum(const eigen::RMatX3 &stiff, const eigen::IColX &indices) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding solid stiff to fluid window.");
    };
    
    virtual void collectStiffFromWindowSum(const eigen::CColX &stiff) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding fluid stiff to solid window.");
    };
    
    virtual void collectStiffFromWindowSum(const eigen::RColX &stiff) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding fluid stiff to solid window.");
    };
    
    virtual void collectStiffFromWindowSum(const eigen::RColX &stiff, const eigen::IColX &indices) {
        throw std::runtime_error("PointWindow::collectStiffFromWindowSum || "
                                 "adding solid stiff to fluid window.");
    };

    /////////////////////////// time loop ///////////////////////////
    
    virtual void randomDispl() = 0;
    virtual void randomStiff() = 0;
    virtual void resetToZero() = 0;
    virtual void resetStiffToZero() = 0;
    
    // check stability
    virtual bool stable() const = 0;
    
    // stiff to accel
    virtual void transformToPhysical() = 0;
    virtual void computeStiffToAccel() = 0;
    virtual void transformToFourier() = 0;
    virtual void applyPressureSource() = 0;
    virtual void maskNyquist() = 0;
    virtual void applyAxialBC() = 0;
    
    const eigen::RColX &getPhiForWindowSum() const {
        return mWindowSumPhi;
    }
    
    const eigen::RColX &getFracForWindowSum() const {
        return mWindowSumFrac;
    }
    
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
    const int mPointMeshTag;
    
    // mass
    const std::unique_ptr<const Mass> mMass;
    
    const eigen::RColX mWindowSumPhi;
    const eigen::RColX mWindowSumFrac;
};

#endif /* Point_hpp */
