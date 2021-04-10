//
//  FluidPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid GLL point

#ifndef FluidPointWindow_hpp
#define FluidPointWindow_hpp

#include "PointWindow.hpp"
#include "eigen_element.hpp"
#include "Scanning1D.hpp"
#include <iostream>
class TimeScheme;

struct FluidStore {
    // pressure source (to be added to accel)
    eigen::CColX mPressureSource = eigen::CColX(0, 1);
    // acceleration storage for pressure output
    // NOTE: duplicated in Newmark
    eigen::CColX mPressureStore = eigen::CColX(0, 1);
    // stiffness storage for delta output
    eigen::CColX mDeltaStore = eigen::CColX(0, 1);
};

class FluidPointWindow: public PointWindow {
public:
    // constructor
    FluidPointWindow(const eigen::RMatX2 &windowSumPhi,
               std::unique_ptr<const Mass> &mass,
               const TimeScheme &timeScheme,
               const std::shared_ptr<Point> point);
    
public:
    bool isFluid() const {return true;};
  
    /////////////////////////// measure ///////////////////////////
    // random displ
    void randomDispl() {
        mFields.mDispl.setRandom();
        mFields.mDispl.row(0).imag().setZero();
    }
    
    // random stiff
    void randomStiff() {
        mFields.mStiff.setRandom();
        mFields.mStiff.row(0).imag().setZero();
    }
    
    // reset to zero
    void resetToZero() {
        mFields.mStiff.setZero();
        mFields.mStiffR.setZero();
        mFields.mDispl.setZero();
        mFields.mVeloc.setZero();
        mFields.mAccel.setZero();
        mFluidStore.mPressureSource.setZero();
        mFluidStore.mPressureStore.setZero();
        mFluidStore.mDeltaStore.setZero();
    }
    
    /////////////////////////// time loop ///////////////////////////
    // check stability
    bool stable() const {
        return mFields.mDispl.allFinite();
    }
    
    // stiff to accel
    void transformToPhysical();
    void computeStiffToAccel();
    void transformToFourier();
    void applyPressureSource();
    void maskNyquist();
    void applyAxialBC();
    
    /////////////////////////// window sum  ///////////////////////////

    eigen::RColX getStiffForWindowSum() {return mFields.mStiffR;};
    
    void collectStiffFromWindowSum(const eigen::RColX &stiff) {
        mFields.mStiffR = stiff.cwiseProduct(mWindowSumFrac);
    }
    
    eigen::RColX getStiffForCommR() {return mFields.mStiffR;};
    eigen::CColX getStiffForCommC() {return mFields.mStiff;};
    void collectStiffFromMessaging(const eigen::RColX &stiff) {mFields.mStiffR = stiff;};
    void collectStiffFromMessaging(const eigen::CColX &stiff) {mFields.mStiff = stiff;};
    
    /////////////////////////// element ///////////////////////////
    // scatter displ to element
    void scatterDisplToElementWindow(eigen::vec_ar1_CMatPP_RM &displ,
                               int nu_1_element, int ipol, int jpol) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            displ[alpha][0](ipol, jpol) = mFields.mDispl(alpha);
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = mNu_1; alpha < nu_1_element; alpha++) {
            displ[alpha][0](ipol, jpol) = czero;
        }
    }
    
    // gather stiff from element
    void gatherStiffFromElementWindow(const eigen::vec_ar1_CMatPP_RM &stiff,
                                int ipol, int jpol) {
        // add lower orders only
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            mFields.mStiff(alpha) -= stiff[alpha][0](ipol, jpol);
        }
    }
    
    
    /////////////////////////// source ///////////////////////////
    // prepare pressure source
    void preparePressureSource() {
        mFluidStore.mPressureSource = eigen::CColX::Zero(mNu_1, 1);
    }
    
    // add pressure source
    void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(mNu_1, nu_1_pressure);
        mFluidStore.mPressureSource.topRows(nu_1_min) +=
        pressure.block(0, ipnt, nu_1_min, 1);
    }
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare pressure output
    void preparePressureOutput() {
        // pressure is mAccel, which may not be allocated by the time scheme
        if (mFluidStore.mPressureStore.rows() == 0) {
            mFluidStore.mPressureStore = eigen::CColX::Zero(mNu_1, 1);
        }
    }
    
    // prepare delta output
    void prepareDeltaOutput() {
        // delta is mStiff, but we need to store it because mStiff is set
        // to zero after dividing by mass
        if (mFluidStore.mDeltaStore.rows() == 0) {
            mFluidStore.mDeltaStore = eigen::CColX::Zero(mNu_1, 1);
        }
    }
    
    // scatter pressure to element
    void scatterPressureToElementWindow(eigen::CMatXN &pressure,
                                  int nu_1_element, int ipnt) const {
        // copy lower orders
        pressure.block(0, ipnt, mNu_1, 1) = mFluidStore.mPressureStore;
        
        // mask higher orders
        pressure.block(mNu_1, ipnt, nu_1_element - mNu_1, 1).setZero();
    }
    
    // scatter delta to element
    void scatterDeltaToElementWindow(eigen::CMatXN &delta,
                               int nu_1_element, int ipnt) const {
        // copy lower orders
        delta.block(0, ipnt, mNu_1, 1) = mFluidStore.mDeltaStore;
        
        // mask higher orders
        delta.block(mNu_1, ipnt, nu_1_element - mNu_1, 1).setZero();
    }
    
    
    /////////////////////////// fields ///////////////////////////
    
    // get
    const Fields<1> &getFluidFields() const {
        return mFields;
    }
    
    // set
    Fields<1> &getFluidFields() {
        return mFields;
    }
    
private:
    // fields on a fluid point
    Fields<1> mFields;
    FluidStore mFluidStore;
    
    /////////////////////////// wavefield scanning ///////////////////////////
public:
    // enable scanning
    void enableScanning() {
        if (!mScanningChi) {
            mScanningChi = std::make_unique<Scanning1D>();
        }
    }
    
    // disable scanning
    void disableScanning() {
        if (mScanningChi) {
            mScanningChi.reset();
            mScanningChi = nullptr;
        }
    }
    
    // do scanning
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) {
        if (mScanningChi) {
            mScanningChi->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                     maxNumPeaks, mFields.mDispl);
        }
    }
    
    // report scanning Nr
    int reportScanningNr() const {
        if (mScanningChi) {
            return mScanningChi->reportScanningNr(mNr);
        } else {
            return -1;
        }
    }
    
private:
    // scanning
    std::unique_ptr<Scanning1D> mScanningChi = nullptr;
};

#endif /* FluidPoint_hpp */
