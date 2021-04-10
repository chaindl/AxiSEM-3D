//
//  SolidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid spectral element

#include "SolidElementWindow.hpp"
#include "SolidPointWindow.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// typeInfo
#include "bstring.hpp"
// measure
#include "timer.hpp"
#include <iostream>
using spectral::nPEM;

// constructor
SolidElementWindow::SolidElementWindow(std::unique_ptr<const GradQuad> &grad,
                           std::unique_ptr<const PRT> &prt,
                           std::unique_ptr<const Elastic> &elastic,
                           const std::array<std::shared_ptr<SolidPointWindow>,
                           spectral::nPEM> &pointWindows, 
                           std::array<eigen::RMatX2, 2> overlapPhi):
ElementWindow(grad, prt, overlapPhi), mElastic(elastic.release()),
mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D()),
mPointWindows(pointWindows) {
    // construct derived
    constructDerived();
    mSE = false;
}

// copy constructor
SolidElementWindow::SolidElementWindow(const SolidElementWindow &other):
ElementWindow(other), mElastic(other.mElastic->clone()),
mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D()),
mPointWindows(other.mPointWindows) {
    // construct derived
    constructDerived();
}

// construct derived
void SolidElementWindow::constructDerived() {
    // point set
    pointWindowSet(mInFourier);
    
    // check compatibility
    mElastic->checkCompatibility(mNr, mInFourier);

    // workspace
    if (sStrainSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    }
    
    // report request to FFT
    if (!mInFourier) {
        if (mPRT) {
            fft::gFFT_N9.addNR(mNr);
        } else {
            fft::gFFT_N6.addNR(mNr);
        }
    }
}

// type info
std::string SolidElementWindow::typeInfo() const {
    std::string info = "SolidElementWindow";
    if (mPRT) {
        if (mPRT->is1D()) {
            info = info + "$PRT1D";
        } else {
            info = info + "$PRT3D";
        }
    }
    if (mElastic->is1D()) {
        info = info + "$" + bstring::typeName(*mElastic) + "1D";
    } else {
        info = info + "$" + bstring::typeName(*mElastic) + "3D";
    }
    return info;
}


/////////////////////////// point ///////////////////////////
// get point
PointWindow &SolidElementWindow::getPointWindow(int ipnt) const {
    return *(mPointWindows[ipnt]);
}

int SolidElementWindow::getMaxNrFromPoints() const {
  int nr = 1;
  for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        nr = std::max(mPointWindows[ipnt]->getNr(), nr);
  }
  return nr;
}

void SolidElementWindow::randomPointWindowsDispl() const {
    for (auto pw: mPointWindows) {
        pw->randomDispl();
    }
}

void SolidElementWindow::resetPointWindowsToZero() const {
    for (auto pw: mPointWindows) {
        pw->resetToZero();
    }
}

/////////////////////////// time loop ///////////////////////////
// collect displacement from points
void SolidElementWindow::
collectDisplFromPointWindows(eigen::vec_ar3_CMatPP_RM &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, mNu_1, ipol, jpol);
        }
    }
}

void SolidElementWindow::displToStrain() const {
    collectDisplFromPointWindows(sDisplSpherical_FR);
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR, sStrainSpherical_FR, mNu_1);
    } else {
        mGradQuad->computeGrad6(sDisplSpherical_FR, sStrainUndulated_FR, mNu_1);
    }
}

void SolidElementWindow::transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);

        } else {
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
        }
    } else {
        if (mElastic->inRTZ()) {
            transform->transformSPZ_RTZ6(sStrainUndulated_FR, mNu_1);
        }
        if (!mInFourier) fft::gFFT_N6.computeC2R(sStrainUndulated_FR,
                                    sStrainUndulated_CD, mNr);
    }
}

void SolidElementWindow::strainToStress() const {
    if (mPRT) {
        if (mInFourier) {
            mElastic->strainToStress_FR(sStrainUndulated_FR,
                                        sStressUndulated_FR, mNu_1);
        } else {
            mElastic->strainToStress_CD(sStrainUndulated_CD,
                                        sStressUndulated_CD, mNr);
        }
    } else {
       if (mInFourier) {
           mElastic->strainToStress_FR(sStrainUndulated_FR,
                                       sStressUndulated_FR, mNu_1);
       } else {
           mElastic->strainToStress_CD(sStrainUndulated_CD,
                                       sStressUndulated_CD, mNr);
       }
    }
}

void SolidElementWindow::transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        if (mInFourier) {
            if (mStressBuffer->size() > 0) {
                for (int alpha = 0; alpha < mNu_1; alpha++) {
                    for (int idim = 0; idim < 6; idim++) {
                        (*mStressBuffer)[alpha][idim] =
                        sStressUndulated_FR[alpha][idim];
                    }
                }
            }
            mPRT->undulatedToSpherical6_FR(sStressUndulated_FR,
                                           sStressSpherical_FR, mNu_1);
        } else {
            // record stress if needed for output
            if (mStressBuffer->size() > 0) {
                //*** additional support required from gFFT_N6 ***//
                fft::gFFT_N6.computeR2C(sStressUndulated_CD,
                                        *mStressBuffer, mNr);
            }
            mPRT->undulatedToSpherical6_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
            fft::gFFT_N9.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
        transform->transformRTZ_SPZ9(sStressSpherical_FR, mNu_1);   
    } else {
        if (!mInFourier) fft::gFFT_N6.computeR2C(sStressUndulated_CD,
                                    sStressUndulated_FR, mNr);
                                    
        // record stress if needed for output
        // record before rotation to keep RTZ consistent with strain
        if (mStressBuffer->size() > 0) {
            for (int alpha = 0; alpha < mNu_1; alpha++) {
                for (int idim = 0; idim < 6; idim++) {
                    (*mStressBuffer)[alpha][idim] =
                    sStressUndulated_FR[alpha][idim];
                }
            }
        }
        
        // stress to stiffness
        if (mElastic->inRTZ()) {
            transform->transformRTZ_SPZ6(sStressUndulated_FR, mNu_1);
        }
    }
}

void SolidElementWindow::stressToStiffness() const {
    if (mPRT) {
        mGradQuad->computeQuad9(sStressSpherical_FR, sStiffSpherical_FR, mNu_1);
    } else {
        mGradQuad->computeQuad6(sStressUndulated_FR, sStiffSpherical_FR, mNu_1);
    }
    addStiffToPointWindows(sStiffSpherical_FR);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
void SolidElementWindow::
addStiffToPointWindows(eigen::vec_ar3_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(stiffElem, ipol, jpol);
        }
    }
}

void SolidElementWindow::getStrainForInterp(eigen::RColX &strain, const int side, const int dim) const {
    int nrows = getNrOverlap(side);
    if (side == 0) {
        strain.topRows(nrows) = sStrainUndulated_CD.col(dim).topRows(nrows);
    } else if (side == 1) {
        strain.topRows(nrows) = sStrainUndulated_CD.col(dim).bottomRows(nrows);
    }
}
    
void SolidElementWindow::addOverlapToStrain(const eigen::RColX &strain, const int side, const int dim) const {
    int nrows = getNrOverlap(side);
    if (side == 0) {
        sStrainUndulated_CD.col(dim).topRows(nrows - 1) += strain.topRows(nrows).bottomRows(nrows - 1);
        sStrainUndulated_CD.col(dim).topRows(nrows) = sStrainUndulated_CD.col(dim).topRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
    } else if (side == 1) {
        sStrainUndulated_CD.col(dim).bottomRows(nrows - 1) += strain.topRows(nrows).topRows(nrows - 1);
        sStrainUndulated_CD.col(dim).bottomRows(nrows) = sStrainUndulated_CD.col(dim).bottomRows(nrows).cwiseProduct(mOverlapPhi[side].col(1));
    }
}

/////////////////////////// source ///////////////////////////
// prepare force source (force given in SPZ)
void SolidElementWindow::prepareForceSource() const {
    // seems nothing
}

// add force source (force given in SPZ)
void SolidElementWindow::addForceSource(const eigen::CMatXN3 &force,
                                  int nu_1_force) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->addForceSource(force, nu_1_force, ipnt);
    }
}

// prepare moment source (moment tensor given in SPZ)
void SolidElementWindow::prepareMomentSource() const {
    // requires FFT N6 even with PRT
    if (mPRT) {
        if (!(mPRT->is1D())) {
            fft::gFFT_N6.addNR(mNr);
        }
    }
    mSE = true;
}

// add moment source (moment tensor given in SPZ)
void SolidElementWindow::addMomentSource(const eigen::CMatXN6 &moment,
                                   int nu_1_moment, const std::shared_ptr<const CoordTransform> &transform) const {
    // pad source with zeros if source has lower order than element
    // truncate source if source has higher order than element
    int nu_1_coexist = std::min(mNu_1, nu_1_moment);
    
    // multiply moment with -1 to convert it to an internal stress
    mapPPvsN::N2PP(-moment, sStressUndulated_FR, nu_1_coexist);
    
    // mask higher orders
    for (int alpha = nu_1_coexist; alpha < mNu_1; alpha++) {
        for (int idim = 0; idim < 6; idim++) {
            sStressUndulated_FR[alpha][idim].setZero();
        }
    }
    
    // by default, source order does not change without 3D PRT
    int nu_1_source = nu_1_coexist;
    
    // stress to stiffness
    if (mPRT) {
        /////////////// with PRT, strain in 3 * 3 ///////////////
        // change to Voigt convention for strain before rotation
        for (int alpha = 0; alpha < nu_1_source; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                sStressUndulated_FR[alpha][idim] *= (numerical::Real)2.;
            }
        }
        transform->transformSPZ_RTZ6(sStressUndulated_FR, nu_1_source);
        for (int alpha = 0; alpha < nu_1_source; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                sStressUndulated_FR[alpha][idim] *= (numerical::Real).5;
            }
        }
        
        // PRT
        if (mPRT->is1D()) {
            mPRT->undulatedToSpherical6_NoIntegration_FR(sStressUndulated_FR,
                                                         sStressSpherical_FR,
                                                         nu_1_source);
        } else {
            //*** additional support required from gFFT_N6 ***//
            fft::gFFT_N6.computeC2R(sStressUndulated_FR,
                                    sStressUndulated_CD, mNr);
            mPRT->undulatedToSpherical6_NoIntegration_CD(sStressUndulated_CD,
                                                         sStressSpherical_CD,
                                                         mNr);
            fft::gFFT_N9.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
            // source order may increase with 3D PRT
            nu_1_source = mNu_1;
        }
        
        // stress to stiffness
        transform->transformRTZ_SPZ9(sStressSpherical_FR, nu_1_source);
        mGradQuad->computeQuad9_NoIntegration(sStressSpherical_FR,
                                              sStiffSpherical_FR, nu_1_source);
    } else {
        /////////////// without PRT, strain in 6 * 1, Voigt ///////////////
        // stress to stiffness
        mGradQuad->computeQuad6_NoIntegration(sStressUndulated_FR,
                                              sStiffSpherical_FR,
                                              nu_1_source);
    }
    
    // mask higher orders
    for (int alpha = nu_1_source; alpha < mNu_1; alpha++) {
        for (int idim = 0; idim < 3; idim++) {
            sStiffSpherical_FR[alpha][idim].setZero();
        }
    }
    
    // add stiffness to points
    // NOTE: doing the following with addStiffToPoints(sStiffSpherical_FR)
    //       is incorrect if the moment tensor is on the injection boundary
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(sStiffSpherical_FR, ipol, jpol);
        }
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
bool SolidElementWindow::
prepareWavefieldOutput(const channel::solid::ChannelOptions &chops) {
    // buffer
    if (chops.mNeedBufferS) {
        mStressBuffer->resize(mNu_1);
    }
    
    // fft
    if (mPRT && !mInFourier && (chops.mNeedBufferE || chops.mNeedBufferS)) {
        fft::gFFT_N6.addNR(mNr);
    }
    
    // coord
    bool needTransform;
    bool needRTZ = (chops.mWCS == channel::WavefieldCS::RTZ);
    if (chops.mNeedBufferU && displInRTZ() != needRTZ) {
        needTransform = true;
    }
    if (chops.mNeedBufferG && nablaInRTZ() != needRTZ) {
        needTransform = true;
    }
    if (chops.mNeedBufferE && strainInRTZ() != needRTZ) {
        needTransform = true;
    }
    if (chops.mNeedBufferR && curlInRTZ() != needRTZ) {
        needTransform = true;
    }
    if (chops.mNeedBufferS && stressInRTZ() != needRTZ) {
        needTransform = true;
    }
    return needTransform;
}

// displ field
void SolidElementWindow::getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);
    // coord
    if (needRTZ) {
        transform->transformSPZ_RTZ3(sDisplSpherical_FR, mNu_1);
    }
    
    // copy
    mapPPvsN::PP2N(sDisplSpherical_FR, displ, mNu_1);
    //if (!mSE && displ.real().norm() > 0) std::cout << displ << std::endl << std::endl;
}

// nabla field
void SolidElementWindow::getNablaField(eigen::CMatXN9 &nabla, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, mNu_1);
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, mNu_1);
        } else {
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, mNu_1);
            fft::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                    sStrainUndulated9_FR, mNr);
        }
        if (!needRTZ) {
            transform->transformRTZ_SPZ9(sStrainUndulated9_FR, mNu_1);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, mNu_1);
        if (needRTZ) {
            transform->transformSPZ_RTZ9(sStrainUndulated9_FR, mNu_1);
        }
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated9_FR, nabla, mNu_1);
}

// strain field
void SolidElementWindow::getStrainField(eigen::CMatXN6 &strain, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {

    displToStrain();
    if (mPRT) {
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
        } else {
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            //*** additional support required from gFFT_N6 ***//
            fft::gFFT_N6.computeR2C(sStrainUndulated_CD,
                                    sStrainUndulated_FR, mNr);
        }
        if (!needRTZ) {
            // change to Voigt convention for strain before rotation
            for (int alpha = 0; alpha < mNu_1; alpha++) {
                // dims 3, 4, 5 are shear components
                for (int idim = 3; idim < 6; idim++) {
                    sStrainUndulated_FR[alpha][idim] *= (numerical::Real).5;
                }
            }
            transform->transformRTZ_SPZ6(sStrainUndulated_FR, mNu_1);
            for (int alpha = 0; alpha < mNu_1; alpha++) {
                // dims 3, 4, 5 are shear components
                for (int idim = 3; idim < 6; idim++) {
                    sStrainUndulated_FR[alpha][idim] *= (numerical::Real)2.;
                }
            }
        }
    } else {
        if (needRTZ) {
            transform->transformSPZ_RTZ6(sStrainUndulated_FR, mNu_1);
        }
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated_FR, strain, mNu_1);
}

// curl field
void SolidElementWindow::getCurlField(eigen::CMatXN3 &curl, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, mNu_1);
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, mNu_1);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, mNu_1);
        } else {
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, mNu_1);
            fft::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                    sStrainUndulated9_FR, mNr);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, mNu_1);
    }
    
    // nabla to curl, use sStiffSpherical_FR as temp storage
    eigen::vec_ar3_CMatPP_RM &sCurlUndulated_FR = sStiffSpherical_FR;
    for (int alpha = 0; alpha < mNu_1; alpha++) {
        sCurlUndulated_FR[alpha][0] = (sStrainUndulated9_FR[alpha][7] -
                                       sStrainUndulated9_FR[alpha][5]);
        sCurlUndulated_FR[alpha][1] = (sStrainUndulated9_FR[alpha][2] -
                                       sStrainUndulated9_FR[alpha][6]);
        sCurlUndulated_FR[alpha][2] = (sStrainUndulated9_FR[alpha][3] -
                                       sStrainUndulated9_FR[alpha][1]);
    }
    
    // coord
    if (needRTZ && !curlInRTZ()) {
        transform->transformSPZ_RTZ3(sCurlUndulated_FR, mNu_1);
    } else if (!needRTZ && curlInRTZ()) {
        transform->transformRTZ_SPZ3(sCurlUndulated_FR, mNu_1);
    }
    
    // copy
    mapPPvsN::PP2N(sCurlUndulated_FR, curl, mNu_1);
}

// stress field
void SolidElementWindow::getStressField(eigen::CMatXN6 &stress, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // coord
    if (needRTZ && !stressInRTZ()) {
        // change to Voigt convention for strain before rotation
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                (*mStressBuffer)[alpha][idim] *= (numerical::Real)2.;
            }
        }
        transform->transformSPZ_RTZ6(*mStressBuffer, mNu_1);
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                (*mStressBuffer)[alpha][idim] *= (numerical::Real).5;
            }
        }
    } else if (!needRTZ && stressInRTZ()) {
        transform->transformRTZ_SPZ6(*mStressBuffer, mNu_1);
    }
    
    // copy
    mapPPvsN::PP2N(*mStressBuffer, stress, mNu_1);
}
