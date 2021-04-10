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
#include <iostream>
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
    if (!mInFourier) {
        fft::gFFT_3.computeC2R(mFields.mStiff, mFields.mStiffR, mNr);
    }
}

void SolidPointWindow::maskNyquist() {
    if (mNr % 2 == 0) mFields.mStiff.bottomRows(1).imag().setZero();
}

// stiff to accel
void SolidPointWindow::computeStiffToAccel() {
    // stiff to accel in-place
    if (mInFourier) {
        mMass->computeAccel(mFields.mStiff);
    } else {
        mMass->computeAccel(mFields.mStiffR);
    }
}

void SolidPointWindow::transformToFourier() {
    if (!mInFourier) fft::gFFT_3.computeR2C(mFields.mStiffR, mFields.mStiff, mNr);
}

void SolidPointWindow::applyAxialBC() {
    static const numerical::ComplexR czero = 0.;
    static const numerical::ComplexR cJ = {0., 1.};
    static const numerical::Real half = .5;
    
    transformToFourier();
        
    // alpha = 0
    mFields.mStiff(0, 0) = czero;
    mFields.mStiff(0, 1) = czero;
    
    // alpha > 0
    if (mNu_1 - 1 >= 1) {
        // alpha = 1
        mFields.mStiff(1, 0) = (mFields.mStiff(1, 0) - cJ * mFields.mStiff(1, 1)) * half;
        mFields.mStiff(1, 1) = cJ * mFields.mStiff(1, 0);
        mFields.mStiff(1, 2) = czero;
        
        // alpha >= 2
        mFields.mStiff.bottomRows(mNu_1 - 2).setZero();
    }
    maskNyquist();
    transformToPhysical();
}
