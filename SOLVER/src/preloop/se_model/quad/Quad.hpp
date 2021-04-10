//
//  Quad.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  Quadrilateral
//  generator of Element in core

#ifndef Quad_hpp
#define Quad_hpp

// eigen
#include "eigen_sem.hpp"

// components
#include "Mapping.hpp"
#include "Material.hpp"
#include "Undulation.hpp"
#include "OceanLoad.hpp"

// external
class ExodusMesh;
class ABC;
class LocalMesh;
class GLLPoint;
class NrField;

// release
class Domain;
class AttBuilder;
class SolidElement;
class FluidElement;
#include "GradientQuadrature.hpp"

// measure
class Element;

class Quad {
public:
    // constructor
    Quad(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
         const NrField &nrField, int localTag);
    
    // finishing 3D properties
    void finishing3D() {
        for (int m = 0; m < getM(); m++) {
            mUndulation[m]->finishing3D();
        }
    }
    
    // finished 3D properties
    void finished3D() {
        for (int m = 0; m < getM(); m++) {
            double winFrac = 2 * numerical::dPi / (std::get<0>(*mWindows[m])(3) - std::get<0>(*mWindows[m])(0));
            mUndulation[m]->finished3D(*this, winFrac);
            mMaterial[m]->finished3D();
        }
    }
    
    // setup GLL
    void setupGLL(const ABC &abc, const LocalMesh &localMesh,
                  std::vector<GLLPoint> &GLLPoints) const;
    
    // compute dt
    double computeDt(double courant, const ABC &abc) const;
    
    // get nodal sz
    const eigen::DMat24 &getNodalSZ() const {
        return mMapping->getNodalSZ();
    }
    
    // release to domain
    void release(const LocalMesh &localMesh,
                 const std::vector<GLLPoint> &GLLPoints,
                 const std::unique_ptr<const AttBuilder> &attBuilder,
                 Domain &domain);
    
    // inverse mapping: (s,z) -> (ξ,η)
    // return true if (s,z) is inside this element
    bool inverseMapping(const eigen::DCol2 &sz, eigen::DCol2 &xieta,
                        double maxIter = 10, double tolerance = 1e-9) const {
        return mMapping->inverseMapping(sz, xieta, maxIter, tolerance);
    }
    
private:
    //////////////////////// interal ////////////////////////
    // compute integral factor
    eigen::DRowN computeIntegralFactor
    (eigen::DMat2N &sz, std::array<eigen::DMat22, spectral::nPEM> &J) const;
    
    // get normal
    eigen::DMat2P computeNormal1D(int edge, const eigen::DMat2N &sz,
                         const std::array<eigen::DMat22, spectral::nPEM> &J,
                         std::vector<int> &ipnts) const;
    
    // weights for CG4 attenuation
    eigen::DRow4 computeWeightsCG4(const eigen::DMatPP_RM &ifPP) const;
    
public:
    //////////////////////// get ////////////////////////
    // get global tag
    inline int getGlobalTag() const {
        return mGlobalTag;
    }
    
    // get point nr
    inline const eigen::IRowN getPointNr(const int m) const {
        return std::get<1>(*mWindows[m]);
    }
    
    eigen::DMatXX computeRelativeWindowPhis(const std::vector<double> &phis) const;
    eigen::DColX computeWindowPhi(int m, int ipnt, bool keep_unwrapped) const; 
    eigen::DColX computeWindowFraction(eigen::DColX phi, int m, bool relative_phi) const;
    int getM() const {return mWindows.size();};

    // get point sz
    eigen::DMat2N getPointSZ() const {
        eigen::DMat2N pointSZ;
        const eigen::DMat2N &xieta = spectrals::getXiEtaElement(axial());
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            pointSZ.col(ipnt) = mMapping->mapping(xieta.col(ipnt));
        }
        return pointSZ;
    }
    
    bool isFluid(double phi) const {
        for (int m = 0; m < mWindows.size(); m++) {
            if (std::get<0>(*mWindows[m])(0) < std::get<0>(*mWindows[m])(3)) {
                if (std::get<0>(*mWindows[m])(0) <= phi && phi <= std::get<0>(*mWindows[m])(3)) {
                    return mFluid3D[m];
                }
            } else {
                if (std::get<0>(*mWindows[m])(3) <= phi || phi <= std::get<0>(*mWindows[m])(0)) {
                    return mFluid3D[m];
                }
            }
        }
    }
    
    bool containsMedium(bool fluid) const {
        for (int m = 0; m < mWindows.size(); m++) {
            if (mFluid3D[m] == fluid) return true;
        }
        return false;
    }
    
    // get undulation
    eigen::arN_DColX getUndulation(const int m) const {
        return mUndulation[m]->getPointwise();
    }
    
    // get material
    std::unique_ptr<Material> &getMaterialPtr(const int m) {
        return mMaterial[m];
    }
    
    // get undulation
    std::unique_ptr<Undulation> &getUndulationPtr(const int m) {
        return mUndulation[m];
    }
    
    // get oceanload
    std::unique_ptr<OceanLoad> &getOceanLoadPtr(const int m) {
        return mOceanLoad[m];
    }
    
    // surface edge
    int getSurfaceEdge() const {
        return mEdgesOnBoundary.at("TOP");
    }
    
    int getSFEdge() const {
        mEdgesOnBoundary.at("SOLID_FLUID");
    }
    
    // axial
    bool axial() const {
        return mEdgesOnBoundary.at("LEFT") != -1;
    }
    
    template <typename floatT>
    std::unique_ptr<const GradientQuadrature<floatT>>
    createGradient(eigen::DMat2N &sz, eigen::DMatPP_RM &ifPP, double winFrac) const {
        // integral factor
        static std::array<eigen::DMat22, spectral::nPEM> J;
        const eigen::DRowN &ifact = computeIntegralFactor(sz, J);
        // save integral factor for later use (CG4)
        ifPP = Eigen::Map<const eigen::DMatPP_RM>(ifact.data());
        
        // compute Jacobian and s^-1
        static eigen::DMatPP_RM dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
        for (int ipol = 0; ipol < spectral::nPED; ipol++) {
            for (int jpol = 0; jpol < spectral::nPED; jpol++) {
                int ipnt = ipol * spectral::nPED + jpol;
                // Jacobian
                double detJ = J[ipnt].determinant();
                dsdxii(ipol, jpol) = J[ipnt](0, 0) / detJ;
                dsdeta(ipol, jpol) = -J[ipnt](0, 1) / detJ;
                dzdxii(ipol, jpol) = -J[ipnt](1, 0) / detJ;
                dzdeta(ipol, jpol) = J[ipnt](1, 1) / detJ;
                // s^-1, use zero when s=0
                double s = sz(0, ipnt);
                inv_s(ipol, jpol) = (axial() && ipol == 0) ? 0. : 1. / s;
            }
        }
        
        // return
        return std::make_unique<GradientQuadrature<floatT>>
        (dsdxii, dsdeta, dzdxii, dzdeta, inv_s, axial(), ifPP, winFrac);
    }
    
    // get element
    std::shared_ptr<Element> getElement() const {
        return mElement;
    }
    
    //////////////////////////////////////////////////////
    //////////////////////// data ////////////////////////
    //////////////////////////////////////////////////////
private:
    // tags
    const int mLocalTag;
    const int mGlobalTag;
    
    // solid-fluid
    const bool mFluid1D;
    std::vector<bool> mFluid3D;
    
    // model boundary
    std::map<std::string, int> mEdgesOnBoundary;
    
    /////////////// nr field ///////////////
    // combined container all the window info:
    // <0> phi: start of window, end of first overlap, start of second overlap, end of window
    // <1> nr
    // <2> tags for GLL access
    // <3> overlap with next window (otherwise unclear for windows with only 1 point overlap)
    std::vector<std::unique_ptr<std::tuple<eigen::DRow4, eigen::IRowN, eigen::IRowN, bool>>> mWindows;
    
    ////////////// components //////////////
    // mapping
    std::unique_ptr<Mapping> mMapping = nullptr;
    // material
    std::vector<std::unique_ptr<Material>> mMaterial;
    // undulation
    std::vector<std::unique_ptr<Undulation>> mUndulation;
    // ocean load
    std::vector<std::unique_ptr<OceanLoad>> mOceanLoad;
    
    // element pointers after release
    std::shared_ptr<Element> mElement = nullptr;
};

#endif /* Quad_hpp */
