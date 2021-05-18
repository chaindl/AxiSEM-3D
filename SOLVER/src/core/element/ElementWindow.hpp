//
//  Element.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element

#ifndef ElementWindow_hpp
#define ElementWindow_hpp

// components
#include "CoordTransform.hpp"
#include "GradientQuadrature.hpp"
#include "PRT.hpp"
#include "eigen_point.hpp"
#include "channel.hpp"
#include "WindowInterpolator.hpp"

#include <iostream>
// Real-typed GradientQuadrature
typedef GradientQuadrature<numerical::Real> GradQuad;
typedef WindowInterpolator<numerical::Real> WinInt;

// point
class PointWindow;

class ElementWindow {
public:
    // constructor
    ElementWindow(std::unique_ptr<const GradQuad> &gq, std::unique_ptr<const PRT> &prt, 
            std::array<eigen::RMatX2, 2> overlapPhi, std::shared_ptr<WinInt> interpolator):
    mGradQuad(gq.release()), mPRT(prt.release()), mOverlapPhi(overlapPhi), mInterpolator(interpolator) {
        
        if (interpolator) mFTBufferNr = interpolation::SplineOrder;
    }
    
    // copy constructor
    ElementWindow(const ElementWindow &other):
    mGradQuad(std::make_unique<GradQuad>(*(other.mGradQuad))),
    mPRT(other.mPRT == nullptr ? nullptr :
         std::make_unique<PRT>(*(other.mPRT))), mOverlapPhi(other.mOverlapPhi),
    mInterpolator(other.mInterpolator) {
        if (other.mFTBufferSplineTag > 0) setUpBufferedFourierTransform();
    }
        
    void pointWindowSet(bool elemInFourier) {
        // order
        mNr = getMaxNrFromPoints();
        mNu_1 = mNr / 2 + 1;
        mNu_1_buffered = (mNr + mFTBufferNr) / 2 + 1;
        
        // check compatibility
        mGradQuad->checkCompatibility(mNr);
        if (mPRT) {
            mPRT->checkCompatibility(mNr, elemInFourier);
        }
    }
    
    // destructor
    virtual ~ElementWindow() = default;
    
    // type info
    virtual std::string typeInfo() const = 0;
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << typeInfo() << "$" << mNr;
        if (mGradQuad->axial()) {
            ss << "$Axial";
        }
        return ss.str();
    }
    
    int getNrOverlap(const int side) const {return mOverlapPhi[side].rows();};
    
    const eigen::RMatX2 &getPhiAndFuncForInterp(const int side) const {return mOverlapPhi[side];};
    virtual int getDimStrain() const = 0;

    /////////////////////////// point ///////////////////////////
    // get point
    virtual PointWindow &getPointWindow(int ipnt) const = 0;
    
    // point set
    void pointSet(bool elemInFourier);
    
    bool isFluid() const {return false;};
    
    // get nr
    int getNr() const {
        return mNr;
    }
    
    // get nu
    int getNu_1() const {
        return mNu_1;
    }
    
    // get nu
    int getNu_1_buffered() const {
        return mNu_1_buffered;
    }
    
    bool elasticInRTZ() const {return false;};
    
    virtual void randomPointWindowsDispl() = 0;
    virtual void resetPointWindowsToZero() = 0;
        
    /////////////////////////// time loop ///////////////////////////
    // displacement to stiffness
    virtual void displToStrain() const = 0;
    virtual void transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform, eigen::RMatXN6 &strain) const = 0;
    virtual void strainToStress(const eigen::RMatXN6 &strain) const = 0;
    virtual void transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const = 0;
    virtual void stressToStiffness(const int tag) const = 0;
    
    template<class Mat>
    void addFTBuffer(Mat &field, eigen::RRowX &BFSMnorm) const {
        BFSMnorm.leftCols(field.cols()) = (field.row(mNr - 1).array() - field.row(0).array()) / (mNr - 1);
        field.topRows(mNr).noalias() -= eigen::RColX::LinSpaced(mNr, 0, mNr - 1) * BFSMnorm.leftCols(field.cols());
        
        mInterpolator->newDataWithIndexing(mFTBufferSplineTag, field, mFTBufferIdx);
        
        mInterpolator->interpolate(mFTBufferSplineTag, field, interpolation::BufferQuery, mNr, field.cols());
    }
    
    /////////////////////////// source ///////////////////////////
    virtual void prepareForceSource() const {
        throw std::runtime_error("ElementWindow::prepareForceSource || force source must be in solid.");
    }
    virtual void addForceSource(const eigen::CMatXN3 &force, int nu_1_force) const {
        throw std::runtime_error("ElementWindow::addForceSource || force source must be in solid.");
    }
    virtual void prepareMomentSource() const {
        throw std::runtime_error("ElementWindow::prepareMomentSource || moment source must be in solid.");
    }
    virtual void addMomentSource(const eigen::CMatXN6 &moment, int nu_1_moment,
                         const std::shared_ptr<const CoordTransform> &transform) const {
        throw std::runtime_error("ElementWindow::addMomentSource || moment source must be in solid.");
    }
    virtual void preparePressureSource() const {
        throw std::runtime_error("ElementWindow::prepareMomentSource || pressure source must be in fluid.");
    }
    virtual void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure) const {
        throw std::runtime_error("ElementWindow::addMomentSource || pressure source must be in fluid.");
    }
    
    /////////////////////////// output ///////////////////////////
    
    virtual bool prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops) {
        throw std::runtime_error("ElementWindow::prepareWavefieldOutput || using fluid options in solid.");
    }
    virtual bool prepareWavefieldOutput(const channel::solid::ChannelOptions &chops) {
        throw std::runtime_error("ElementWindow::prepareWavefieldOutput || using solid options in fluid.");
    }
                                                          
    // checks for coord transform
    virtual bool displInRTZ() const = 0;
    virtual bool nablaInRTZ() const {
        throw std::runtime_error("ElementWindow::nablaInRTZ || unavailable in fluid.");
    }
    virtual bool strainInRTZ() const {
        throw std::runtime_error("ElementWindow::strainInRTZ || unavailable in fluid.");
    }
    virtual bool curlInRTZ() const {
        throw std::runtime_error("ElementWindow::curlInRTZ || unavailable in fluid.");
    }
    virtual bool stressInRTZ() const {
        throw std::runtime_error("ElementWindow::stressInRTZ || unavailable in fluid.");
    }
    
    // record
    // solid and fluid
    virtual void getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const = 0;
    // fluid only
    virtual void getChiField(eigen::CMatXN &chi) const {
        throw std::runtime_error("ElementWindow::getChiField || unavailable in solid.");
    }
    virtual void getPressureField(eigen::CMatXN &pressure) const {
        throw std::runtime_error("ElementWindow::getPressureField || unavailable in solid.");
    }
    virtual void getDeltaField(eigen::CMatXN &delta) const {
        throw std::runtime_error("ElementWindow::getDeltaField || unavailable in solid.");
    }
    // solid only
    virtual void getNablaField(eigen::CMatXN9 &nabla, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
        throw std::runtime_error("ElementWindow::getNableField || unavailable in fluid.");
    }
    virtual void getStrainField(eigen::CMatXN6 &strain, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
        throw std::runtime_error("ElementWindow::getStrainField || unavailable in fluid.");
    }
    virtual void getCurlField(eigen::CMatXN3 &curl, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
        throw std::runtime_error("ElementWindow::getCurlField || unavailable in fluid.");
    }
    virtual void getStressField(eigen::CMatXN6 &stress, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const {
        throw std::runtime_error("ElementWindow::getStressField || unavailable in fluid.");
    }
    
protected:
    void setUpBufferedFourierTransform() {
        mFTBufferIdx = eigen::IColX::Zero(mFTBufferNr, 1);
        
        eigen::RColX knots(mFTBufferNr);
        knots << eigen::RColX::LinSpaced(mFTBufferNr / 2, 0, mFTBufferNr / 2 - 1), 
                 eigen::RColX::LinSpaced(mFTBufferNr / 2, 3 * mFTBufferNr / 2, 4 * mFTBufferNr / 2 - 1);
                 
        knots.array() /= 2 * mFTBufferNr - 1;
        for (int n = 0; n < mFTBufferNr / 2; n++) {
            mFTBufferIdx(n) = mNr - mFTBufferNr / 2 + n;
            mFTBufferIdx(mFTBufferNr / 2 + n) = n;
        }
        mFTBufferSplineTag = mInterpolator->addSplineFitting(knots, (isFluid() ? 3 * spectral::nPEM : 6 * spectral::nPEM));
    }

private:    
    virtual int getMaxNrFromPoints() const = 0;
    
    eigen::IColX mFTBufferIdx;
    int mFTBufferSplineTag = -1;
    
protected:
    const std::shared_ptr<WinInt> mInterpolator;
    
    // gradient operator
    const std::unique_ptr<const GradQuad> mGradQuad = nullptr;
    
    // particle relabelling
    const std::unique_ptr<const PRT> mPRT = nullptr;
    
    // order
    int mNr = 0;
    int mFTBufferNr = 0;
    int mNu_1 = 0;
    int mNu_1_buffered = 0;
    
    // this is phi (scaled) and frac on left side <0> and phi and frac on right side <1>
    std::array<eigen::RMatX2, 2> mOverlapPhi;
    
    mutable eigen::RMatXX tester1 = eigen::RMatXX::Zero(8, 6 * spectral::nPEM);
    mutable eigen::IColX tester2 = eigen::IColX::Zero(8);
};

#endif /* Element_hpp */
