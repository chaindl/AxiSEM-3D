//
//  SolidPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid GLL point

#include "SolidPointWindow.hpp"
#include "point_time.hpp"
#include "fft.hpp"

// constructor
SolidPointWindow::SolidPointWindow(const eigen::RMatX2 &windowSumPhi, 
                       std::unique_ptr<const Mass> &mass,
                       const TimeScheme &timeScheme,
                       const std::shared_ptr<Point> point):
PointWindow(windowSumPhi, mass, point) {
    // fields
    point_time::createFields(*this, timeScheme);
    // mass
    mMass->checkCompatibility(mNr, false);
    fft::gFFT_3.addNR(mNr);
}

/////////////////////////// time loop ///////////////////////////
void SolidPointWindow::transformToPhysical() {
    // Nyquist
    if (mNr % 2 == 0) {
        mFields.mStiff.bottomRows(1).imag().setZero();
    }
  
    if (!mInFourier) fft::gFFT_3.computeC2R(mFields.mStiff, mFields.mStiffR, mNr);
}

// stiff to accel
void SolidPointWindow::computeStiffToAccel() {
    // stiff to accel in-place
    mMass->computeAccel(mFields.mStiffR);
}

void SolidPointWindow::transformToFourier() {
    if (!mInFourier) fft::gFFT_3.computeR2C(mFields.mStiffR, mFields.mStiff, mNr);
}

void SolidPointWindow::applyPressureSource() {
    throw std::runtime_error("SolidPointWindow::applyPressureSource || "
                                     "Incompatible types: "
                                     "pressure source and solid point window.");
}
