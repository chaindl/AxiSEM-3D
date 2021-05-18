
#include "FluidElementWindow.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// typeInfo
#include "bstring.hpp"
// measure
#include "timer.hpp"

using spectral::nPEM;

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::constructDerived() {
    // point set
    pointWindowSet(mInFourier);
    
    // check compatibility
    mAcoustic->checkCompatibility(mNr, mInFourier);
    
    // workspace
    if (sStressSpherical_CD.rows() < mNr) {
        expandWorkspace(mNr);
    } 
    
    // report request to FFT
    if (!mInFourier) {
        fft::gFFT_N3.addNR(mNr);
    }
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::constructDerived() {
    // point set
    pointWindowSet(false);
    
    // check compatibility
    mAcoustic->checkCompatibility(mNr, false);
    
    // workspace
    if (sStressSpherical_CD.rows() < mNr + mFTBufferNr) {
        expandWorkspace(mNr + mFTBufferNr);
        expandWorkspaceBFSM(mNr + mFTBufferNr);
    } 
    
    // report request to FFT
    fft::gFFT_N3.addNR(mNr + mFTBufferNr);
    fft::gFFT_N1.addNR(mNr + mFTBufferNr);
    
    setUpBufferedFourierTransform();
}

// type info
template <class FluidPointWindow>
std::string FluidElementWindow<FluidPointWindow>::typeInfo() const {
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
template <class FluidPointWindow>
PointWindow &FluidElementWindow<FluidPointWindow>::getPointWindow(int ipnt) const {
    return *(mPointWindows[ipnt]);
}

template <class FluidPointWindow>
int FluidElementWindow<FluidPointWindow>::getMaxNrFromPoints() const {
  int nr = 1;
  for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        nr = std::max(mPointWindows[ipnt]->getNr(), nr);
  }
  return nr;
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::randomPointWindowsDispl() {
    for (auto pw: mPointWindows) {
        pw->randomDispl();
    }
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::resetPointWindowsToZero() {
    for (auto pw: mPointWindows) {
        pw->resetToZero();
    }
}
/////////////////////////// time loop ///////////////////////////
// collect displacement from points
template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::
collectDisplFromPointWindows(eigen::vec_ar1_CMatPP_RM &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, mNu_1, ipol, jpol);
        }
    }
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::
collectDisplFromPointWindows(eigen::RMatXN &displElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            scatterDisplToElementWindow(displElem, ipol * spectral::nPED + jpol);
        }
    }
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::
displToStrain() const {
    collectDisplFromPointWindows(sDisplSpherical_FR);

    mGradQuad->computeGrad3(sDisplSpherical_FR,
                            sStrainSpherical_FR, sBFSMnorm_zero, mNu_1_buffered);
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::
displToStrain() const {
    collectDisplFromPointWindows(sDisplSpherical_CD);
    addFTBuffer(sDisplSpherical_CD, sBFSMnorm);
    fft::gFFT_N1.computeR2C(sDisplSpherical_CD, sDisplSpherical_FR, mNr + mFTBufferNr);
    mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 1);

    mGradQuad->computeGrad3(sDisplSpherical_FR,
                            sStrainSpherical_FR, sBFSMnorm_PP, mNu_1_buffered);
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::
transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        transform->transformSPZ_RTZ3(sStrainSpherical_FR, mNu_1_buffered);
        if (mInFourier) {
            mPRT->sphericalToUndulated3_FR(sStrainSpherical_FR,
                                           sStrainUndulated_FR, mNu_1);
        } else {
            fft::gFFT_N3.computeC2R(sStrainSpherical_FR, sStrainSpherical_CD, mNr + mFTBufferNr);
            mPRT->sphericalToUndulated3_CD(sStrainSpherical_CD,
                                           sStrainUndulated_CD, mNr);
        }
    } else {
        if (!mInFourier) fft::gFFT_N3.computeC2R(sStrainSpherical_FR, sStrainSpherical_CD, mNr + mFTBufferNr);
    }
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::
transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform, eigen::RMatXN6 &strain) const {
    transformStrainToPhysical(transform);
    if (!mInFourier) {
        if (mPRT) {
            strain.block(0, 0, spectral::nPEM * 3, mNr) = sStrainUndulated_CD.topRows(mNr);
        } else {
            strain.block(0, 0, spectral::nPEM * 3, mNr) = sStrainSpherical_CD.topRows(mNr);
        }
    }
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::strainToStress() const {
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

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::strainToStress(const eigen::RMatXN6 &strain) const {
    if (!mInFourier) {
        if (mPRT) {
          sStrainUndulated_CD.topRows(mNr) = strain.block(0, 0, spectral::nPEM * 3, mNr);
        } else {
          sStrainSpherical_CD.topRows(mNr) = strain.block(0, 0, spectral::nPEM * 3, mNr);
        }
    }
    strainToStress();
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::
transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
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
        if (!mInFourier) {
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
    }
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::
transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const {
    if (mPRT) {
        mPRT->undulatedToSpherical3_CD(sStressUndulated_CD,
                                           sStressSpherical_CD, mNr);
        addFTBuffer(sStressSpherical_CD, sBFSMnorm);
        mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 3);
        fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr + mFTBufferNr);     
        transform->transformRTZ_SPZ3(sStressSpherical_FR, mNu_1_buffered);
        transform->transformRTZ_SPZ3(sBFSMnorm_PP);
    } else {
        addFTBuffer(sStressSpherical_CD, sBFSMnorm);
        mapPPvsN::N2PP_buf(sBFSMnorm, sBFSMnorm_PP, 3);
        fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr + mFTBufferNr);
    }
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::
transformStressToFourier_noBuffer(const std::shared_ptr<const CoordTransform> &transform) const {
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
        if (!mInFourier) {
            fft::gFFT_N3.computeR2C(sStressSpherical_CD,
                                    sStressSpherical_FR, mNr);
        }
    }
}

template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::stressToStiffness(const int tag) const {
    const eigen::ar3_RMatPP_RM &BFSMnorm = (mInterpolator) ? sBFSMnorm_PP : sBFSMnorm_zero;
    
    mGradQuad->computeQuad3(sStressSpherical_FR,
                            sStiffSpherical_FR, BFSMnorm, mNu_1_buffered);
    addStiffToPointWindows(sStiffSpherical_FR);
}

// add stiffness to points
// allow a derived class to change stiffElem (no const)
template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::
addStiffToPointWindows(const eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(stiffElem, ipol, jpol);
        }
    }
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::
addStiffToPointWindows(const eigen::vec_ar1_CMatPP_RM &stiffElem) const {
    fft::gFFT_N1.computeC2R(stiffElem, sDisplSpherical_CD, mNr + mFTBufferNr);
    for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
            mPointWindows[ipol * spectral::nPED + jpol]->
            gatherStiffFromElementWindow(sDisplSpherical_CD, ipol * spectral::nPED + jpol);
        }
    }
}


/////////////////////////// source ///////////////////////////
// prepare pressure source
template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::preparePressureSource() const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->preparePressureSource();
    }
}

// add pressure source
template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::addPressureSource(const eigen::CMatXN &pressure,
                                     int nu_1_pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->addPressureSource(pressure, nu_1_pressure, ipnt);
    }
}


/////////////////////////// wavefield output ///////////////////////////
// prepare wavefield output
template <class FluidPointWindow>
bool FluidElementWindow<FluidPointWindow>::
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
    
    // unless recording wavefields, we only use transform with buffered nr
    if (chops.mNeedBufferX && mInterpolator) fft::gFFT_N1.addNR(mNr);
    if (chops.mNeedBufferU && mInterpolator) fft::gFFT_N3.addNR(mNr);

    return needTransform;
}

// chi field
template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>::getChiField(eigen::CMatXN &chi) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_FR);

    // convert to flattened
    mapPPvsN::PP2N(sDisplSpherical_FR, chi, mNu_1);
}

template <>
void FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>::getChiField(eigen::CMatXN &chi) const {
    // collect displacement from points
    collectDisplFromPointWindows(sDisplSpherical_CD);
    fft::gFFT_N1.computeR2C(sDisplSpherical_CD, sDisplSpherical_FR, mNr);
    
    // convert to flattened
    mapPPvsN::PP2N(sDisplSpherical_FR, chi, mNu_1);
}

// displ field
template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
  
    displToStrain();
    transformStrainToPhysical(transform);
    strainToStress();
    transformStressToFourier_noBuffer(transform);
    
    if (mPRT && !needRTZ) transform->transformRTZ_SPZ3(sStressSpherical_FR, mNu_1);
    if (!mPRT && needRTZ) transform->transformSPZ_RTZ3(sStressSpherical_FR, mNu_1);
    
    // convert to flattened
    mapPPvsN::PP2N(sStressSpherical_FR, displ, mNu_1);
}

// pressure field
template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::getPressureField(eigen::CMatXN &pressure) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->scatterPressureToElementWindow(pressure, mNu_1, ipnt);
    }
}

// delta field
template <class FluidPointWindow>
void FluidElementWindow<FluidPointWindow>::getDeltaField(eigen::CMatXN &delta) const {
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        mPointWindows[ipnt]->scatterDeltaToElementWindow(delta, mNu_1, ipnt);
    }
}

template class FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>;
template class FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>;
