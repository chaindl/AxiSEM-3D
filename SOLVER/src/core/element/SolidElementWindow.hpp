//
//  SolidElement.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid spectral element

#ifndef SolidElementWindow_hpp
#define SolidElementWindow_hpp

#include "ElementWindow.hpp"
#include "SFRCPointWindow.hpp"
// output
#include "mapPPvsN.hpp"
// measure
#include "timer.hpp"
// point
#include <array>
class PointWindow;

// material
#include "Elastic.hpp"

// output
#include "channel.hpp"

template <class SolidPointWindow>
class SolidElementWindow: public ElementWindow {
public:
    // constructor
    SolidElementWindow(std::unique_ptr<const GradQuad> &grad,
                 std::unique_ptr<const PRT> &prt,
                 std::unique_ptr<const Elastic> &elastic,
                 const std::array<std::shared_ptr<SolidPointWindow>, spectral::nPEM> &pointWindows, 
                 std::array<eigen::RMatX2, 2> overlapPhi,
                 std::shared_ptr<WinInt> interpolator):
    ElementWindow(grad, prt, overlapPhi, interpolator), mElastic(elastic.release()),
    mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D() && !interpolator),
    mPointWindows(pointWindows) {
        // construct derived
        constructDerived();
    };
    
    // copy constructor
    SolidElementWindow(const SolidElementWindow &other):
    ElementWindow(other), mElastic(other.mElastic->clone()),
    mInFourier((mPRT ? mPRT->is1D() : true) && mElastic->is1D() && !mInterpolator),
    mPointWindows(other.mPointWindows) {
        // construct derived
        constructDerived();
    };
    
private:
    // construct derived
    void constructDerived();
    
public:
    /////////////////////////// info //////////////////////////
    std::string typeInfo() const;
    std::string mediumInfo() const {return "SOLID";};
    int getMaxNrFromPoints() const;
    int getDimStrain() const {return spectral::nPEM * 6;};
    bool elasticInRTZ() const {return mElastic->inRTZ();};
    
    /////////////////////////// pointer access //////////////////////////
    PointWindow &getPointWindow(int inpt) const;
    
    /////////////////////////// time loop ///////////////////////////
private:
    // collect displacement from points (only called internally by displToStrain)
    void collectDisplAndBuffer() const;
    void collectDisplFromPointWindows(eigen::vec_ar3_CMatPP_RM &displElem) const;
    void collectDisplFromPointWindows(eigen::RMatXN3 &displElem) const;
    // scatter stiffness to points (only called internally by stressToStiffness)
    void addStiffToPointWindows(const eigen::vec_ar3_CMatPP_RM &stiffElem) const;
    void addStiffToPointWindows(const eigen::RMatX3 &stiffElem) const;
    void transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform) const;
    
public:
    // displacement to stiffness
    void displToStrain() const;
    void transformStrainToPhysical(const std::shared_ptr<const CoordTransform> &transform, eigen::RMatXN6 &strain) const;
    void strainToStress(const eigen::RMatXN6 &strain) const;
    void transformStressToFourier(const std::shared_ptr<const CoordTransform> &transform) const;
    void stressToStiffness(const int tag) const;
    // for measuring cost
    void randomPointWindowsDispl();
    void resetPointWindowsToZero();
    
    /////////////////////////// source ///////////////////////////
    // force given in SPZ
    void prepareForceSource() const;
    void addForceSource(const eigen::CMatXN3 &force,
                        int nu_1_force) const;
                        
    // moment tensor given in SPZ
    void prepareMomentSource() const;
    void addMomentSource(const eigen::CMatXN6 &moment,
                         int nu_1_moment,
                         const std::shared_ptr<const CoordTransform> &transform) const;
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare wavefield output
    bool prepareWavefieldOutput(const channel::solid::ChannelOptions &chops);
    
    // get fields
    void getDisplField(eigen::CMatXN3 &displ, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ) const;
    void getNablaField(eigen::CMatXN9 &nabla, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const;
    void getStrainField(eigen::CMatXN6 &strain, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const;
    void getCurlField(eigen::CMatXN3 &curl, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const;
    void getStressField(eigen::CMatXN6 &stress, const std::shared_ptr<const CoordTransform> &transform, bool needRTZ, int &nu) const;
    
    // displ crd
    bool displInRTZ() const {
        return false;
    }
    
    // nabla crd
    bool nablaInRTZ() const {
        return (bool)mPRT;
    }
    
    // strain crd
    bool strainInRTZ() const {
        return nablaInRTZ();
    }
    
    // curl crd
    bool curlInRTZ() const {
        return nablaInRTZ();
    }
    
    // stress crd
    bool stressInRTZ() const {
        return strainInRTZ() || mElastic->inRTZ();
    }
    
private:
    // material
    const std::unique_ptr<const Elastic> mElastic;
    
    // 1D element in Fourier space
    const bool mInFourier;

    // points
    std::array<std::shared_ptr<SolidPointWindow>, spectral::nPEM> mPointWindows;
    
    // stress buffer
    // stress cannot be recomputed with attenuation
    const std::unique_ptr<eigen::vec_ar6_CMatPP_RM> mStressBuffer =
    std::make_unique<eigen::vec_ar6_CMatPP_RM>();
    
    
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
        sStrainSpherical_CD.resize(maxNr, spectral::nPEM * 9);
        sStrainUndulated_CD.resize(maxNr, spectral::nPEM * 6);
        sStressUndulated_CD.resize(maxNr, spectral::nPEM * 6);
        sStressSpherical_CD.resize(maxNr, spectral::nPEM * 9);
    }
    
        // expand workspace
    static void expandWorkspaceBFSM(int maxNr) {
        int maxNu_1 = maxNr / 2 + 1;
        
        sDisplSpherical_CD.resize(maxNr, spectral::nPEM * 3);
        sBFSMnorm = eigen::RRowX::Zero(1, spectral::nPEM * 9);
    } 

    // workspace
    // Fourier
    inline static eigen::vec_ar3_CMatPP_RM sDisplSpherical_FR;
    inline static eigen::vec_ar9_CMatPP_RM sStrainSpherical_FR;
    inline static eigen::vec_ar6_CMatPP_RM sStrainUndulated_FR;
    inline static eigen::vec_ar6_CMatPP_RM sStressUndulated_FR;
    inline static eigen::vec_ar9_CMatPP_RM sStressSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStiffSpherical_FR;
    // cardinal
    inline static eigen::RMatXN3 sDisplSpherical_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN9 sStrainSpherical_CD =
    eigen::RMatXN9(0, spectral::nPEM * 9);
    inline static eigen::RMatXN6 sStrainUndulated_CD =
    eigen::RMatXN6(0, spectral::nPEM * 6);
    inline static eigen::RMatXN6 sStressUndulated_CD =
    eigen::RMatXN6(0, spectral::nPEM * 6);
    inline static eigen::RMatXN9 sStressSpherical_CD =
    eigen::RMatXN9(0, spectral::nPEM * 9);
    
    inline static eigen::RRowX sBFSMnorm;
    inline static eigen::ar9_RMatPP_RM sBFSMnorm_PP;
    const inline static eigen::ar9_RMatPP_RM sBFSMnorm_zero = {eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero(),
                                                               eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero(),
                                                               eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero(), eigen::RMatPP_RM::Zero()};
                                                               
    bool mPrepared = false;
    int mNrPrepped = -1;
};

#endif /* SolidElement_hpp */
