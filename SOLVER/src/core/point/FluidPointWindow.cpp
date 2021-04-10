//
//  FluidPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid GLL point

#include "FluidPointWindow.hpp"
#include "point_time.hpp"
#include "fft.hpp"

// constructor
FluidPointWindow::FluidPointWindow(const eigen::RMatX2 &windowSumPhi, 
                       std::unique_ptr<const Mass> &mass,
                       const TimeScheme &timeScheme,
                       const std::shared_ptr<Point> point):
PointWindow(windowSumPhi, mass, point) {
    // fields
    point_time::createFields(*this, timeScheme);
    // mass
    mMass->checkCompatibility(mNr, false);
    fft::gFFT_1.addNR(mNr);
}

/////////////////////////// time loop ///////////////////////////
void FluidPointWindow::transformToPhysical() {
    if (!mInFourier) {
        fft::gFFT_1.computeC2R(mFields.mStiff, mFields.mStiffR, mNr);
    }
}

void FluidPointWindow::maskNyquist() {
    if (mNr % 2 == 0) mFields.mStiff.bottomRows(1).imag().setZero();
}

// stiff to accel
void FluidPointWindow::computeStiffToAccel() {
    // store stiffness for delta output
    if (mFluidStore.mDeltaStore.rows() > 0) {
        if (mInFourier) {
            mFluidStore.mDeltaStore = mFields.mStiff;
        } else {
            fft::gFFT_1.computeR2C(mFields.mStiffR, mFluidStore.mDeltaStore, mNr);
        }
    }
    
    // stiff to accel in-place
    if (mInFourier) {
        mMass->computeAccel(mFields.mStiff);
    } else {
        mMass->computeAccel(mFields.mStiffR);
    }
}

void FluidPointWindow::transformToFourier() {
    if (!mInFourier) fft::gFFT_1.computeR2C(mFields.mStiffR, mFields.mStiff, mNr);
}

void FluidPointWindow::applyPressureSource() {
    // apply pressure source
    if (mFluidStore.mPressureSource.rows() > 0) {
        // add pressure to new acceleration
        mFields.mStiff += mFluidStore.mPressureSource;
        // zero pressure for the next time step
        mFluidStore.mPressureSource.setZero();
    }
    
    // store acceleration for pressure output
    if (mFluidStore.mPressureStore.rows() > 0) {
        mFluidStore.mPressureStore = mFields.mStiff;
    }
}

void FluidPointWindow::applyAxialBC() {    
    transformToFourier();
    mFields.mStiff.bottomRows(mNu_1 - 1).setZero();
    maskNyquist();
    transformToPhysical();
}
