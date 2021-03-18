//
//  FluidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#include "FluidElementWindow.hpp"
#include "FluidPointWindow.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// measure
#include "timer.hpp"

using spectral::nPEM;

// constructor
FluidElementWindow::FluidElementWindow(std::unique_ptr<const GradQuad> &grad,
                                       std::unique_ptr<const PRT> &prt,
                                       std::unique_ptr<const Acoustic> &acoustic,
                                       const std::array<std::shared_ptr<FluidPointWindow>, spectral::nPEM> &pointWindows, 
                                       std::array<eigen::RMatX2, 2> &overlapPhi):
ElementWindow(grad, prt, overlapPhi), mAcoustic(acoustic.release()),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D() && getNrOverlap(0) + getNrOverlap(1) == 0),
mPointWindows(pointWindows) {
    // construct derived
    constructDerived();
}

// copy constructor
FluidElementWindow::FluidElementWindow(const FluidElementWindow &other):
ElementWindow(other), mAcoustic(std::make_unique<Acoustic>(*(other.mAcoustic))),
mInFourier((mPRT ? mPRT->is1D() : true) && mAcoustic->is1D() && getNrOverlap(0) + getNrOverlap(1) == 0),
mPointWindows(other.mPointWindows) {
    // construct derived
    constructDerived();
}

// construct derived
void FluidElementWindow::constructDerived() {
    // point set
    pointWindowSet(mInFourier);
    
    // check compatibility
    mAcoustic->checkCompatibility(mNr, mInFourier);
    
    // workspace
    if (sStrainSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    }
    
    // report request to FFT
    if (!mInFourier) {
        fft::gFFT_N3.addNR(mNr);
    }
}

// type info
std::string FluidElementWindow::typeInfo() const {
    std::string info = "FluidElementWindow";
    if (mPRT) {
        if (mPRT->is1D()) {
            info = info + "$PRT1D";
        } else {
            info = info + "$PRT3D";
        }
    }
    if (mAcoustic->is1D()) {
        info = info + "$Acoustic1D";
    } else {
        info = info + "$Acoustic3D";
    }
    return info;
}


/////////////////////////// point ///////////////////////////
// get point
PointWindow &FluidElementWindow::getPointWindow(int ipnt) const {
    return *(mPointWindows[ipnt]);
}

int FluidElementWindow::getMaxNrFromPoints() const {
  int nr = 1;
  for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        nr = std::max(mPointWindows[ipnt]->getNr(), nr);
  }
  return nr;
}

void FluidElementWindow::randomPointWindowsDispl() const {
    for (auto pw: mPointWindows) {
        pw->randomDispl();
    }
}

void FluidElementWindow::resetPointWindowsToZero() const {
    for (auto pw: mPointWindows) {
        pw->resetToZero();
    }
}
/////////////////////////// time loop ///////////////////////////
// collect displacement from points
void FluidElementWindow::
collectDisplFromPointWindows(eigen::vec_ar1_CMatPP_RM &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, mNu_1, ipol, jpol);
        }
    }
}

void FluidElementWindow::displToStrain() const {
    collectDisplFromPointWindows(sDisplSpherical_FR);
    mGradQuad->computeGrad3(sDisplSpherical_FR,
                            sStrainSpherical_FR, mNu_1);
}

void FluidElementWindow::transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        transform->transformSPZ_RTZ3(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR, sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
        }
    } else {
       if (!mInFourier) fft::gFFT_N3.computeC2R(sStrainSpherical_FR, sStrainSpherical_CD, mNr);
    }
}

void FluidElementWindow::strainToStress() const {
    if (mPRT) {
        if (mInFourier) {
            mAcoustic->strainToStress_FR(sStrainUndulated_FR,
                                         sStressUndulated_FR, mNu_1);
        } else {
            mAcoustic->strainToStress_CD(sStrainUndulated_CD,
                                         sStressUndulated_CD, mNr);
        }
    } else {
       if (mInFourier) {
           mAcoustic->strainToStress_FR(sStrainSpherical_FR,
                                         sStressSpherical_FR, mNu_1);
       } else {
           mAcoustic->strainToStress_CD(sStrainSpherical_CD,
                                         sStressSpherical_CD, mNr);
       }
    }
}

void FluidElementWindow::transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        if (mInFourier) {
            mPRT->undulatedToSpherical3_FR(sStressUndulated_FR,
                                           sStressSpherical_FR, mNu_1);
        } else {
            mPRT->undulatedToSpherical3_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }        
        transform->transformRTZ_SPZ3(sStressSpherical_FR, mNu_1);
    } else {
       if (!mInFourier) fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
    }
}

void FluidElementWindow::stressToStiffness() const {
    mGradQuad->computeQuad3(sStressSpherical_FR,
                            sStiffSpherical_FR, mNu_1);
    addStiffToPointWindows(sStiffSpherical_FR);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
void FluidElementWindow::
addStiffToPointWindows(eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(sStiffSpherical_FR, ipol, jpol);
        }
    }
}

void FluidElementWindow::getStrainForInterp(eigen::RColX &strain, const int side, const int dim) const {
    int nrows = getNrOverlap(side);
    if (mPRT) {
        if (side == 0) {
            strain.topRows(nrows) = sStrainUndulated_CD.col(dim).topRows(nrows);
        } else if (side == 1) {
            strain.topRows(nrows) = sStrainUndulated_CD.col(dim).bottomRows(nrows);
        }
    } else {
        if (side == 0) {
            strain.topRows(nrows) = sStrainSpherical_CD.col(dim).topRows(nrows);
        } else if (side == 1) {
            strain.topRows(nrows) = sStrainSpherical_CD.col(dim).bottomRows(nrows);
        }
    }
}
    
void FluidElementWindow::addOverlapToStrain(const eigen::RColX &strain, const int side, const int dim) const {
    int nrows = getNrOverlap(side);
    if (mPRT) {
        if (side == 0) {
            sStrainUndulated_CD.col(dim).topRows(nrows - 1) += strain.topRows(nrows).bottomRows(nrows - 1);
            sStrainUndulated_CD.col(dim).topRows(nrows) = sStrainUndulated_CD.col(dim).topRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
        } else if (side == 1) {
            sStrainUndulated_CD.col(dim).bottomRows(nrows - 1) += strain.topRows(nrows).topRows(nrows - 1);
            sStrainUndulated_CD.col(dim).bottomRows(nrows) = sStrainUndulated_CD.col(dim).bottomRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
        }
    } else {
        if (side == 0) {
            sStrainSpherical_CD.col(dim).topRows(nrows - 1) += strain.topRows(nrows).bottomRows(nrows - 1);
            sStrainSpherical_CD.col(dim).topRows(nrows) = sStrainSpherical_CD.col(dim).topRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
        } else if (side == 1) {
            sStrainSpherical_CD.col(dim).bottomRows(nrows - 1) += strain.topRows(nrows).topRows(nrows - 1);
            sStrainSpherical_CD.col(dim).bottomRows(nrows) = sStrainSpherical_CD.col(dim).bottomRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
        }
    }
}


/////////////////////////// source ///////////////////////////
// prepare pressure source
void FluidElementWindow::preparePressureSource() const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->preparePressureSource();
    }
}

// add pressure source
void FluidElementWindow::addPressureSource(const eigen::CMatXN &pressure,
                                     int nu_1_pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->addPressureSource(pressure, nu_1_pressure, ipnt);
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
bool FluidElementWindow::
prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops) {
    // pressure
    if (chops.mNeedBufferP) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPointWindows[ipnt]->preparePressureOutput();
        }
    }
    
    // delta
    if (chops.mNeedBufferD) {
        for (int ipnt = 0; ipnt < nPEM; ipnt++) {
            mPointWindows[ipnt]->prepareDeltaOutput();
        }
    }
    
    // coord
    bool needTransform = false;
    bool needRTZ = (chops.mWCS == channel::WavefieldCS::RTZ);
    if (chops.mNeedBufferU && displInRTZ() != needRTZ) {
        needTransform = true;
    }
    return needTransform;
}

// chi field
void FluidElementWindow::getChiField(eigen::CMatXN &chi) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);
    
    // convert to flattened
    mapPPvsN::PP2N(sDisplSpherical_FR, chi, mNu_1);
}

// displ field
void FluidElementWindow::getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
  
    displToStrain();
    transformStrainToPhysical(transform);
    strainToStress();
    transformStressToFourier(transform);
    
    if (mPRT && !needRTZ) transform->transformRTZ_SPZ3(sStressSpherical_FR, mNu_1);
    if (!mPRT && needRTZ) transform->transformSPZ_RTZ3(sStressSpherical_FR, mNu_1);
    
    // convert to flattened
    mapPPvsN::PP2N(sStressUndulated_FR, displ, mNu_1);
}

// pressure field
void FluidElementWindow::getPressureField(eigen::CMatXN &pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->scatterPressureToElementWindow(pressure, mNu_1, ipnt);
    }
}

// delta field
void FluidElementWindow::getDeltaField(eigen::CMatXN &delta) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->scatterDeltaToElementWindow(delta, mNu_1, ipnt);
    }
}
