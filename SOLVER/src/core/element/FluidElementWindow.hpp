//
//  FluidElement.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#ifndef FluidElementWindow_hpp
#define FluidElementWindow_hpp

#include "ElementWindow.hpp"

// point
#include <array>
class PointWindow;
class FluidPointWindow;

// material
#include "Acoustic.hpp"

// output
#include "channel.hpp"

class FluidElementWindow: public ElementWindow {
public:
    // constructor
    FluidElementWindow(std::unique_ptr<const GradQuad> &grad,
                       std::unique_ptr<const PRT> &prt,
                       std::unique_ptr<const Acoustic> &acoustic,
                       const std::array<std::shared_ptr<FluidPointWindow>, spectral::nPEM> &pointWindows, 
                       std::array<eigen::RMatX2, 2> &overlapPhi);
    
    // copy constructor
    FluidElementWindow(const FluidElementWindow &other);
    
private:
    // construct derived
    void constructDerived();
    
public:
    /////////////////////////// info //////////////////////////
    std::string typeInfo() const;
    std::string mediumInfo() const {return "FLUID";};
    int getMaxNrFromPoints() const;
    bool isFluid() const {return true;};
    int getDimStrain() const {return spectral::nPEM * 3;};
    
    /////////////////////////// pointer access //////////////////////////
    FluidElementWindow &getFluidElementWindow() {
        return *this;
    };
    PointWindow &getPointWindow(int inpt) const;
    
    /////////////////////////// time loop ///////////////////////////
private:
    // collect displacement from points (only called internally by displToStrain)
    void collectDisplFromPointWindows(eigen::vec_ar1_CMatPP_RM &displElem) const;
    // scatter stiffness to points (only called internally by stressToStiffness)
    void addStiffToPointWindows(eigen::vec_ar1_CMatPP_RM &stiffElem) const;
    
public:
    // displacement to stiffness
    void displToStrain() const;
    void transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const;
    void getStrainForInterp(eigen::RColX &strain, const int side, const int dim) const;
    void addOverlapToStrain(const eigen::RColX &strain, const int side, const int dim) const;
    void strainToStress() const;
    void transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const;
    void stressToStiffness() const;
    
    // for measuring cost
    void randomPointWindowsDispl() const;
    void resetPointWindowsToZero() const;
  
    /////////////////////////// source ///////////////////////////
    
    void preparePressureSource() const;
    void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure) const;
    
    /////////////////////////// wavefield output ///////////////////////////
    
    bool prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops);
    
    // get fields
    void getChiField(eigen::CMatXN &chi) const;
    void getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const;
    void getPressureField(eigen::CMatXN &pressure) const;
    void getDeltaField(eigen::CMatXN &delta) const;
    
    // displ crd
    bool displInRTZ() const {
        return (bool)mPRT;
    }
    
private:
    // material
    const std::unique_ptr<const Acoustic> mAcoustic;
    
    // 1D element window without overlap in Fourier space
    const bool mInFourier;
    
    // points
    std::array<std::shared_ptr<FluidPointWindow>, spectral::nPEM> mPointWindows;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // expand workspace
    static void expandWorkspace(int maxNr) {
        int maxNu_1 = maxNr / 2 + 1;
        // Fourier
        sDisplSpherical_FR.resize(maxNu_1);
        sStrainSpherical_FR.resize(maxNu_1);
        sStrainUndulated_FR.resize(maxNu_1);
        sStressUndulated_FR.resize(maxNu_1);
        sStressSpherical_FR.resize(maxNu_1);
        sStiffSpherical_FR.resize(maxNu_1);
        // cardinal
        sStrainSpherical_CD.resize(maxNr, spectral::nPEM * 3);
        sStrainUndulated_CD.resize(maxNr, spectral::nPEM * 3);
        sStressUndulated_CD.resize(maxNr, spectral::nPEM * 3);
        sStressSpherical_CD.resize(maxNr, spectral::nPEM * 3);
    }
    
    // workspace
    // Fourier
    inline static eigen::vec_ar1_CMatPP_RM sDisplSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStrainSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStrainUndulated_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStressUndulated_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStressSpherical_FR;
    inline static eigen::vec_ar1_CMatPP_RM sStiffSpherical_FR;
    // cardinal
    inline static eigen::RMatXN3 sStrainSpherical_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStrainUndulated_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStressUndulated_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStressSpherical_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
};

#endif /* FluidElement_hpp */
