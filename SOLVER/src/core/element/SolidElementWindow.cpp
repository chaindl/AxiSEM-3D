//
//  SolidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid spectral element

#include "SolidElementWindow.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// typeInfo
#include "bstring.hpp"
// measure
#include "timer.hpp"
#include <iostream>
using spectral::nPEM;

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::constructDerived() {
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

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::constructDerived() {
    // point set
    pointWindowSet(false);
    // check compatibility
    mElastic->checkCompatibility(mNr, false);
    // workspace
    if (sStrainSpherical_CD.rows() < mNr + mFTBufferNr) {
        expandWorkspace(mNr + mFTBufferNr);
        expandWorkspaceBFSM(mNr + mFTBufferNr);
    } 
    // report request to FFT
    if (mPRT) {
        fft::gFFT_N9.addNR(mNr + mFTBufferNr);
    } else {
        fft::gFFT_N6.addNR(mNr + mFTBufferNr);
    }
    fft::gFFT_N3.addNR(mNr + mFTBufferNr);
    
    setUpBufferedFourierTransform();
}

// type info
template <class SolidPointWindow>
std::string SolidElementWindow<SolidPointWindow>::typeInfo() const {
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
template <class SolidPointWindow>
PointWindow &SolidElementWindow<SolidPointWindow>::getPointWindow(int ipnt) const {
    return *(mPointWindows[ipnt]);
}

template <class SolidPointWindow>
int SolidElementWindow<SolidPointWindow>::getMaxNrFromPoints() const {
  int nr = 1;
  for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        nr = std::max(mPointWindows[ipnt]->getNr(), nr);
  }
  return nr;
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::randomPointWindowsDispl() {
    for (auto pw: mPointWindows) {
        pw->randomDispl();
    }
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::resetPointWindowsToZero() {
    for (auto pw: mPointWindows) {
        pw->resetToZero();
    }
}

/////////////////////////// time loop ///////////////////////////
// collect displacement from points
template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::
collectDisplFromPointWindows(eigen::vec_ar3_CMatPP_RM &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, mNu_1, ipol, jpol);
        }
    }
}

// collect displacement from points
template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::
collectDisplFromPointWindows(eigen::RMatXN3 &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, ipol * spectral::nPED + jpol);
        }
    }
}

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::
collectDisplAndBuffer() const {
    collectDisplFromPointWindows(sDisplSpherical_FR);
    sBFSMnorm_PP = sBFSMnorm_zero;
}

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::
collectDisplAndBuffer() const {

    collectDisplFromPointWindows(sDisplSpherical_CD);

    addFTBuffer(sDisplSpherical_CD, sBFSMnorm);

    fft::gFFT_N3.computeR2C(sDisplSpherical_CD, sDisplSpherical_FR, mNr + mFTBufferNr);

    mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 3);

}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::displToStrain() const {
    collectDisplAndBuffer();
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR, sStrainSpherical_FR, sBFSMnorm_PP, mNu_1_buffered);
    } else {
        mGradQuad->computeGrad6(sDisplSpherical_FR, sStrainUndulated_FR, sBFSMnorm_PP, mNu_1_buffered);
    }
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, mNu_1_buffered);
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
        } else {
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr + mFTBufferNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
        }
    } else {
        if (mElastic->inRTZ()) {
            transform->transformSPZ_RTZ6(sStrainUndulated_FR, mNu_1_buffered);
        }
        if (!mInFourier) {
            fft::gFFT_N6.computeC2R(sStrainUndulated_FR,
                                                 sStrainUndulated_CD, mNr + mFTBufferNr);
        }
    }
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform, eigen::RMatXN6 &strain) const {
    transformStrainToPhysical(transform);
    strain.topRows(mNr) = sStrainUndulated_CD.topRows(mNr);
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
strainToStress(const eigen::RMatXN6 &strain) const {
    if (mPRT) {
        if (mInFourier) {
            mElastic->strainToStress_FR(sStrainUndulated_FR,
                                        sStressUndulated_FR, mNu_1);
        } else {
            mElastic->strainToStress_CD(strain,
                                        sStressUndulated_CD, mNr);
        }
    } else {
       if (mInFourier) {
           mElastic->strainToStress_FR(sStrainUndulated_FR,
                                       sStressUndulated_FR, mNu_1);
       } else {
           mElastic->strainToStress_CD(strain,
                                       sStressUndulated_CD, mNr);
       }
    }
}

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::
transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
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
        if (!mInFourier) {
            fft::gFFT_N6.computeR2C(sStressUndulated_CD,
                                    sStressUndulated_FR, mNr);                      
        }
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

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::
transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        // record stress if needed for output
        if (mStressBuffer->size() > 0) {
            //*** additional support required from gFFT_N6 ***//
            fft::gFFT_N6.computeR2C(sStressUndulated_CD,
                                    *mStressBuffer, mNr);
        }
        mPRT->undulatedToSpherical6_CD(sStressUndulated_CD,
                                       sStressSpherical_CD, mNr);
        addFTBuffer(sStressSpherical_CD, sBFSMnorm);
        mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 9);
        fft::gFFT_N9.computeR2C(sStressSpherical_CD,
                                sStressSpherical_FR, mNr + mFTBufferNr);
        transform->transformRTZ_SPZ9(sStressSpherical_FR, mNu_1_buffered);
        transform->transformRTZ_SPZ9(sBFSMnorm_PP);
    } else {
        addFTBuffer(sStressUndulated_CD, sBFSMnorm);
        mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 6);
        fft::gFFT_N6.computeR2C(sStressUndulated_CD,
                                    sStressUndulated_FR, mNr + mFTBufferNr);                      
        // record stress if needed for output
        // record before rotation to keep RTZ consistent with strain
        if (mStressBuffer->size() > 0) {
            for (int alpha = 0; alpha < mNu_1_buffered; alpha++) {
                for (int idim = 0; idim < 6; idim++) {
                    (*mStressBuffer)[alpha][idim] =
                    sStressUndulated_FR[alpha][idim];
                }
            }
        }
        
        // stress to stiffness
        if (mElastic->inRTZ()) {
            transform->transformRTZ_SPZ6(sStressUndulated_FR, mNu_1_buffered);
        }
    }
}

template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
stressToStiffness(const int tag) const {
    const eigen::ar9_RMatPP_RM &BFSMnorm = (mInterpolator) ? sBFSMnorm_PP : sBFSMnorm_zero;
  
    if (mPRT) {
        mGradQuad->computeQuad9(sStressSpherical_FR, sStiffSpherical_FR, BFSMnorm, mNu_1_buffered);
    } else {
        mGradQuad->computeQuad6(sStressUndulated_FR, sStiffSpherical_FR, BFSMnorm, mNu_1_buffered);
    }
    addStiffToPointWindows(sStiffSpherical_FR);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::
addStiffToPointWindows(const eigen::vec_ar3_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(stiffElem, ipol, jpol);
        }
    }
}

template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::
addStiffToPointWindows(const eigen::vec_ar3_CMatPP_RM &stiffElem) const {
    fft::gFFT_N3.computeC2R(stiffElem, sDisplSpherical_CD, mNr + mFTBufferNr);
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(sDisplSpherical_CD, ipol * spectral::nPED + jpol);
        }
    }
}

/////////////////////////// source ///////////////////////////
// prepare force source (force given in SPZ)
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::prepareForceSource() const {
    // seems nothing
}

// add force source (force given in SPZ)
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::addForceSource(const eigen::CMatXN3 &force,
                                  int nu_1_force) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->addForceSource(force, nu_1_force, ipnt);
    }
}

// prepare moment source (moment tensor given in SPZ)
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::prepareMomentSource() const {
    // requires FFT N6 even with PRT
    if (mPRT) {
        if (!(mPRT->is1D())) {
            fft::gFFT_N6.addNR(mNr);
            fft::gFFT_N9.addNR(mNr);
        }
    }
}

// add moment source (moment tensor given in SPZ)
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::addMomentSource(const eigen::CMatXN6 &moment,
                                   int nu_1_moment, const std::shared_ptr<const CoordTransform> &transform) const {
    // pad source with zeros if source has lower order than element
    // truncate source if source has higher order than element
    int nu_1_coexist = std::min(mNu_1_buffered, nu_1_moment);
    
    // multiply moment with -1 to convert it to an internal stress
    mapPPvsN::N2PP(-moment, sStressUndulated_FR, nu_1_coexist);
    
    // mask higher orders
    for (int alpha = nu_1_coexist; alpha < mNu_1_buffered; alpha++) {
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
                                    sStressUndulated_CD, mNr + mFTBufferNr);
            mPRT->undulatedToSpherical6_NoIntegration_CD(sStressUndulated_CD,
                                                         sStressSpherical_CD,
                                                         mNr);
            if (mInterpolator) {
                addFTBuffer(sStressSpherical_CD, sBFSMnorm);
                mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 9);
            } else {
                sBFSMnorm_PP = sBFSMnorm_zero;
            }
                                                         
            fft::gFFT_N9.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr + mFTBufferNr);
            // source order may increase with 3D PRT
            nu_1_source = mNu_1;
        }
        
        // stress to stiffness
        transform->transformRTZ_SPZ9(sStressSpherical_FR, nu_1_source);
        mGradQuad->computeQuad9_NoIntegration(sStressSpherical_FR,
                                              sStiffSpherical_FR, 
                                              sBFSMnorm_PP, nu_1_source);
    } else {
        /////////////// without PRT, strain in 6 * 1, Voigt ///////////////
        // stress to stiffness
        mGradQuad->computeQuad6_NoIntegration(sStressUndulated_FR,
                                              sStiffSpherical_FR,
                                              sBFSMnorm_zero, nu_1_source);
    }
    
    // mask higher orders
    for (int alpha = nu_1_source; alpha < mNu_1; alpha++) {
        for (int idim = 0; idim < 3; idim++) {
            sStiffSpherical_FR[alpha][idim].setZero();
        }
    }
    
    addStiffToPointWindows(sStiffSpherical_FR);
}

/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
template <class SolidPointWindow>
bool SolidElementWindow<SolidPointWindow>::
prepareWavefieldOutput(const channel::solid::ChannelOptions &chops) {
    mPrepared = true;
    // buffer
    if (chops.mNeedBufferS) {
        mStressBuffer->resize(mNu_1);
    }
    
    // fft
    if (mPRT && !mInFourier && (chops.mNeedBufferE || chops.mNeedBufferS)) {
        fft::gFFT_N6.addNR(mNr);
    }
    if (mPRT && !mInFourier && (chops.mNeedBufferG || chops.mNeedBufferR)) {
        fft::gFFT_N9.addNR(mNr);
    }
    if (mInterpolator && chops.mNeedBufferU) {
        fft::gFFT_N3.addNR(mNr);
        mNrPrepped = mNr;
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
template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>::
getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);
    // coord
    if (needRTZ) transform->transformSPZ_RTZ3(sDisplSpherical_FR, mNu_1);
    // copy
    mapPPvsN::PP2N(sDisplSpherical_FR, displ, mNu_1);
}

// displ field
template <>
void SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>::
getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_CD);
    fft::gFFT_N3.computeR2C(sDisplSpherical_CD, sDisplSpherical_FR, mNr);
    // coord
    if (needRTZ) transform->transformSPZ_RTZ3(sDisplSpherical_FR, mNu_1);
    // copy
    mapPPvsN::PP2N(sDisplSpherical_FR, displ, mNu_1);
}

// nabla field
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
getNablaField(eigen::CMatXN9 &nabla, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
    // collect displacement from points
    collectDisplAndBuffer();
    nu = mNu_1_buffered;
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, sBFSMnorm_PP, nu);
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, nu);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, nu);
        } else {
            nu = mNu_1;
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr + mFTBufferNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, mNr);
            fft::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                    sStrainUndulated9_FR, mNr);
        }
        if (!needRTZ) {
            transform->transformRTZ_SPZ9(sStrainUndulated9_FR, nu);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, sBFSMnorm_PP, nu);
        if (needRTZ) {
            transform->transformSPZ_RTZ9(sStrainUndulated9_FR, nu);
        }
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated9_FR, nabla, nu);
}

// strain field
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
getStrainField(eigen::CMatXN6 &strain, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
    displToStrain();
    nu = mNu_1_buffered;
    
    if (mPRT) {
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, nu);
        if (mInFourier) {
            mPRT->sphericalToUndulated6_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, nu);
        } else {
            nu = mNu_1;
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr + mFTBufferNr);
            mPRT->sphericalToUndulated6_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
            //*** additional support required from gFFT_N6 ***//
            fft::gFFT_N6.computeR2C(sStrainUndulated_CD,
                                    sStrainUndulated_FR, mNr);
        }
        if (!needRTZ) {
            // change to Voigt convention for strain before rotation
            for (int alpha = 0; alpha < nu; alpha++) {
                // dims 3, 4, 5 are shear components
                for (int idim = 3; idim < 6; idim++) {
                    sStrainUndulated_FR[alpha][idim] *= (numerical::Real).5;
                }
            }
            transform->transformRTZ_SPZ6(sStrainUndulated_FR, nu);
            for (int alpha = 0; alpha < nu; alpha++) {
                // dims 3, 4, 5 are shear components
                for (int idim = 3; idim < 6; idim++) {
                    sStrainUndulated_FR[alpha][idim] *= (numerical::Real)2.;
                }
            }
        }
    } else {
        if (needRTZ) {
            transform->transformSPZ_RTZ6(sStrainUndulated_FR, nu);
        }
    }
    
    // copy
    mapPPvsN::PP2N(sStrainUndulated_FR, strain, nu);
}

// curl field
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
getCurlField(eigen::CMatXN3 &curl, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
    // collect displacement from points
    collectDisplAndBuffer();
    nu = mNu_1_buffered;
    
    // displ to nabla, use sStressSpherical as temp storage
    eigen::vec_ar9_CMatPP_RM &sStrainUndulated9_FR = sStressSpherical_FR;
    eigen::RMatXN9 &sStrainUndulated9_CD = sStressSpherical_CD;
    if (mPRT) {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainSpherical_FR, sBFSMnorm_PP, nu);
        transform->transformSPZ_RTZ9(sStrainSpherical_FR, nu);
        if (mInFourier) {
            mPRT->sphericalToUndulated9_FR(sStrainSpherical_FR,
                                           sStrainUndulated9_FR, nu);
        } else {
            nu = mNu_1;
            fft::gFFT_N9.computeC2R(sStrainSpherical_FR,
                                    sStrainSpherical_CD, mNr + mFTBufferNr);
            mPRT->sphericalToUndulated9_CD(sStrainSpherical_CD,
                                           sStrainUndulated9_CD, mNr);
            fft::gFFT_N9.computeR2C(sStrainUndulated9_CD,
                                    sStrainUndulated9_FR, mNr);
        }
    } else {
        mGradQuad->computeGrad9(sDisplSpherical_FR,
                                sStrainUndulated9_FR, sBFSMnorm_PP, nu);
    }
    
    // nabla to curl, use sStiffSpherical_FR as temp storage
    eigen::vec_ar3_CMatPP_RM &sCurlUndulated_FR = sStiffSpherical_FR;
    for (int alpha = 0; alpha < nu; alpha++) {
        sCurlUndulated_FR[alpha][0] = (sStrainUndulated9_FR[alpha][7] -
                                       sStrainUndulated9_FR[alpha][5]);
        sCurlUndulated_FR[alpha][1] = (sStrainUndulated9_FR[alpha][2] -
                                       sStrainUndulated9_FR[alpha][6]);
        sCurlUndulated_FR[alpha][2] = (sStrainUndulated9_FR[alpha][3] -
                                       sStrainUndulated9_FR[alpha][1]);
    }
    
    // coord
    if (needRTZ && !curlInRTZ()) {
        transform->transformSPZ_RTZ3(sCurlUndulated_FR, nu);
    } else if (!needRTZ && curlInRTZ()) {
        transform->transformRTZ_SPZ3(sCurlUndulated_FR, nu);
    }
    
    // copy
    mapPPvsN::PP2N(sCurlUndulated_FR, curl, nu);
}

// stress field
template <class SolidPointWindow>
void SolidElementWindow<SolidPointWindow>::
getStressField(eigen::CMatXN6 &stress, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
    nu = mNu_1;
    
    // coord
    if (needRTZ && !stressInRTZ()) {
        // change to Voigt convention for strain before rotation
        for (int alpha = 0; alpha < nu; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                (*mStressBuffer)[alpha][idim] *= (numerical::Real)2.;
            }
        }
        transform->transformSPZ_RTZ6(*mStressBuffer, nu);
        for (int alpha = 0; alpha < nu; alpha++) {
            // dims 3, 4, 5 are shear components
            for (int idim = 3; idim < 6; idim++) {
                (*mStressBuffer)[alpha][idim] *= (numerical::Real).5;
            }
        }
    } else if (!needRTZ && stressInRTZ()) {
        transform->transformRTZ_SPZ6(*mStressBuffer, nu);
    }
    
    // copy
    mapPPvsN::PP2N(*mStressBuffer, stress, nu);
}

template class SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>;
template class SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>;
