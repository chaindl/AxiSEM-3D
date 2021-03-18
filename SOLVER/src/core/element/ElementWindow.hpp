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
// Real-typed GradientQuadrature
typedef GradientQuadrature<numerical::Real> GradQuad;

// point
class PointWindow;
class FluidElementWindow;
class SolidElementWindow;

class ElementWindow {
public:
    // constructor
    ElementWindow(std::unique_ptr<const GradQuad> &gq,
            std::unique_ptr<const PRT> &prt, std::array<eigen::RMatX2, 2> overlapPhi):
    mGradQuad(gq.release()), mPRT(prt.release()), mOverlapPhi(overlapPhi) {
        // nothing
    }
    
    // copy constructor
    ElementWindow(const ElementWindow &other):
    mGradQuad(std::make_unique<GradQuad>(*(other.mGradQuad))),
    mPRT(other.mPRT == nullptr ? nullptr :
         std::make_unique<PRT>(*(other.mPRT))) {
        // nothing 
    }
    
    void pointWindowSet(bool elemInFourier) {
        // order
        mNr = getMaxNrFromPoints();
        mNu_1 = mNr / 2 + 1;
        
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
    
    int getNrOverlap(int side) const {return mOverlapPhi[side].rows();};
    
    eigen::RColX getPhiForInterp(const int side) const {return mOverlapPhi[side].col(0);};
    virtual void getStrainForInterp(eigen::RColX &strain, const int side, const int dim) const = 0;
    virtual int getDimStrain() const = 0;
    void setAlignment(bool aligned) {mAligned = aligned;};
    /////////////////////////// point ///////////////////////////
    // get point
    virtual PointWindow &getPointWindow(int ipnt) const = 0;
    
    bool overlapIsAligned() const {return mAligned;};
    
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
    
    bool elasticInRTZ() const {return false;};
    
    virtual void randomPointWindowsDispl() const = 0;
    virtual void resetPointWindowsToZero() const = 0;
    
    virtual FluidElementWindow &getFluidElementWindow() const {
        throw std::runtime_error("ElementWindow::getFluidElementWindow || "
                                 "requested fluid window from solid medium");
    };
    virtual SolidElementWindow &getSolidElementWindow() const {
        throw std::runtime_error("ElementWindow::getSolidElementWindow || "
                                 "requested solid window from fluid medium");
    };
    
        
    /////////////////////////// time loop ///////////////////////////
    // displacement to stiffness
    virtual void displToStrain() const = 0;
    virtual void transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const = 0;
    virtual void strainToStress() const = 0;
    virtual void transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const = 0;
    virtual void stressToStiffness() const = 0;
    virtual void addOverlapToStrain(const eigen::RColX &strain, const int side, const int dim) const = 0;
    
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
    virtual void getNablaField(eigen::CMatXN9 &nabla, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
        throw std::runtime_error("ElementWindow::getNableField || unavailable in fluid.");
    }
    virtual void getStrainField(eigen::CMatXN6 &strain, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
        throw std::runtime_error("ElementWindow::getStrainField || unavailable in fluid.");
    }
    virtual void getCurlField(eigen::CMatXN3 &curl, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
        throw std::runtime_error("ElementWindow::getCurlField || unavailable in fluid.");
    }
    virtual void getStressField(eigen::CMatXN6 &stress, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const {
        throw std::runtime_error("ElementWindow::getStressField || unavailable in fluid.");
    }
    
private:
    virtual int getMaxNrFromPoints() const = 0;
    
protected:
    // gradient operator
    const std::unique_ptr<const GradQuad> mGradQuad;
    
    // particle relabelling
    const std::unique_ptr<const PRT> mPRT;
    
    // order
    int mNr = 0;
    int mNu_1 = 0;
    
    // this is phi and frac on left side <0> and phi and frac on right side <1>
    std::array<eigen::RMatX2, 2> mOverlapPhi;
    bool mAligned;
};

#endif /* Element_hpp */
