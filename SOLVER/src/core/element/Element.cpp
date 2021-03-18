//
//  FluidElement.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#include "Element.hpp"
#include "Point.hpp"
#include "fft.hpp"
// output
#include "mapPPvsN.hpp"
// measure
#include "timer.hpp"

#include "vicinity.hpp"
#include "geodesy.hpp"
#include "CoordTransformCartesian.hpp"
#include "CoordTransformSpherical.hpp"

using spectral::nPEM;

// constructor
Element::Element(const int &quadTag, std::vector<std::unique_ptr<ElementWindow>> &windows,
      const std::array<std::shared_ptr<Point>, spectral::nPEM> &points):
mQuadTag(quadTag), mWindows(std::move(windows)), mPoints(points) {

    int maxNrOverlap = -1;
    bool needTransform = false;
    for (auto &win: mWindows) {
        maxNrOverlap = std::max({maxNrOverlap, win->getNrOverlap(0), win->getNrOverlap(1)});
        if (win->elasticInRTZ()) needTransform = true;
    }
    expandWorkspace(maxNrOverlap);
    if (needTransform) createCoordTransform();
}

// copy constructor
Element::Element(const Element &other):
mQuadTag(other.mQuadTag), mWindows(other.copyWindows()), mPoints(other.mPoints) {

    int maxNrOverlap = -1;
    for (auto &win: mWindows) {
        maxNrOverlap = std::max({maxNrOverlap, win->getNrOverlap(0), win->getNrOverlap(1)});
    }
    expandWorkspace(maxNrOverlap);
}

/////////////////////////// point ///////////////////////////
// get point
Point &Element::getPoint(int ipnt) const {
    return *(mPoints[ipnt]);
}

const eigen::DRow2 &Element::getPointCoords(int ipnt) const {
    return mPoints[ipnt]->getCoords();
}

// find boundary points by tag
std::vector<int> Element::
findBoundaryPointsByTag(const std::vector<int> &boundaryMeshTags) const {
    // point indices on an element boundary
    const std::vector<int> &allSidePoints = vicinity::constants::gEdgeIPntAll;
    
    // result
    std::vector<int> pointsFound;
    pointsFound.reserve(allSidePoints.size());
    
    // check mesh tag
    for (int ipnt: allSidePoints) {
        if (vector_tools::findSortedUnique(boundaryMeshTags,
                                           getPoint(ipnt).getMeshTag())) {
            // found a point on the boundary
            pointsFound.push_back(ipnt);
        }
    }
    
    // already sorted and unique
    return pointsFound;
}

// find boundary points by crds
std::vector<int> Element::
findBoundaryPointsByCrds(const std::vector<double> &boundaryCrdsRorZ,
                         const std::vector<double> &boundaryCrdsTorS,
                         double distTol) const {
    // point indices on an element boundary
    const std::vector<int> &allSidePoints = vicinity::constants::gEdgeIPntAll;
    
    // collect coords
    eigen::DColX crdRorZ(nPEM);
    eigen::DColX crdTorS(nPEM);
    for (int ipnt = 0; ipnt < nPEM; ipnt++) {
        const eigen::DRow2 &sz = getPoint(ipnt).getCoords();
        if (geodesy::isCartesian()) {
            crdRorZ(ipnt) = sz(1);
            crdTorS(ipnt) = sz(0);
        } else {
            crdRorZ(ipnt) = sz.norm();
            crdTorS(ipnt) = (crdRorZ(ipnt) < numerical::dEpsilon) ? 0. :
            acos(sz(1) / crdRorZ(ipnt));
        }
    }
    
    // result
    std::vector<int> pointsFound;
    
    // check r or z
    for (double bRorZ: boundaryCrdsRorZ) {
        for (int ipnt: allSidePoints) {
            if (std::abs(crdRorZ[ipnt] - bRorZ) < distTol) {
                pointsFound.push_back(ipnt);
            }
        }
    }
    
    // check theta or s
    for (double bTorS: boundaryCrdsTorS) {
        for (int ipnt: allSidePoints) {
            double dist = std::abs(crdTorS[ipnt] - bTorS);
            if (!geodesy::isCartesian()) {
                // change angular distance to meter
                dist *= crdRorZ(ipnt);
            }
            if (dist < distTol) {
                pointsFound.push_back(ipnt);
            }
        }
    }
    
    // sort and unique points
    vector_tools::sortUnique(pointsFound);
    return pointsFound;
}

/////////////////////////// crd transform ///////////////////////////
// create coordinate transform
// internally called by setToRTZ_ByMaterialOrPRT
// externally called by prepareWavefieldOutput
void Element::createCoordTransform() {
    if (mTransform) return;
    // difference between Cartesian and spherical
    if (geodesy::isCartesian()) {
        mTransform = std::make_unique<const CoordTransformCartesian>();
    } else {
        // compute theta
        eigen::DMatPP_RM theta;
        for (int ipol = 0; ipol < spectral::nPED; ipol++) {
            for (int jpol = 0; jpol < spectral::nPED; jpol++) {
                int ipnt = ipol * spectral::nPED + jpol;
                const eigen::DRow2 &sz = getPoint(ipnt).getCoords();
                const eigen::DRow2 &rt = geodesy::sz2rtheta(sz, true);
                theta(ipol, jpol) = rt(1);
            }
        }
        mTransform = std::make_shared<const CoordTransformSpherical>(theta);
    }
}

/////////////////////////// time loop ///////////////////////////
// displacement to stiffness
void Element::computeStiff() const {
               
    for (int m = 0; m < mWindows.size(); m++) {
        mWindows[m]->displToStrain();
        mWindows[m]->transformStrainToPhysical(mTransform);
    }
    
    overlapAndAddStrain();
    
    for (int m = 0; m < mWindows.size(); m++) {
        mWindows[m]->strainToStress();
        mWindows[m]->transformStressToFourier(mTransform);
        mWindows[m]->stressToStiffness();
    }
}

void Element::overlapAndAddStrain() const {
    if (mWindows.size() == 1) return;
    
    for (int m1 = 0; m1 < mWindows.size(); m1++) {
        if (mWindows[m1]->getNrOverlap(1) == 0) continue;
        
        int m2 = (m1 + 1) % mWindows.size();
        int dim = mWindows[m1]->getDimStrain(); // 9N for solid, 3N for fluid
        if (dim != mWindows[m2]->getDimStrain()) {
            throw std::runtime_error("Element::overlapAndAddStrain || Attemting overlap of solid and fluid window.");
        }
        
        for (int i = 0; i < dim; i++) {
            mWindows[m1]->getStrainForInterp(sWin1ToWin2, 1, i);
            mWindows[m2]->getStrainForInterp(sWin2ToWin1, 0, i);
            
            if (!mWindows[m1]->overlapIsAligned()) {
                interpolate(mWindows[m2]->getPhiForInterp(0), sWin1ToWin2, 
                                          mWindows[m1]->getPhiForInterp(1), 0);
                interpolate(mWindows[m1]->getPhiForInterp(1), sWin2ToWin1, 
                                          mWindows[m2]->getPhiForInterp(0), 1);
            }
                            
            mWindows[m1]->addOverlapToStrain(sWin2ToWin1, 1, i);
            mWindows[m2]->addOverlapToStrain(sWin1ToWin2, 0, i);
        }
    }
}

void Element::interpolate(const eigen::RColX &phi_q, eigen::RColX &fun, const eigen::RColX &phi, const int side) const {
        double dphi = phi(1) - phi(0);
        if (dphi < 0) dphi += 2 * numerical::dPi;
        
        boost::math::interpolators::cardinal_cubic_b_spline<numerical::Real> spline(fun.topRows(phi.rows()).begin(), 
                                                                                    fun.topRows(phi.rows()).end(), 
                                                                                    phi(0), dphi);
        
        double phi_with_wrapping = -1.;
        for (int i = side; i < phi_q.rows() + side - 1; i++) {
            phi_with_wrapping = (phi_with_wrapping > phi_q(i)) ? phi_q(i) + 2 * numerical::dPi : phi_q(i);
            fun(i) = spline(phi_with_wrapping);  
        }
}

// measure cost
double Element::measure(int count) const {
    // random displacement
    for (auto &win: mWindows) {
        win->randomPointWindowsDispl();
    }
    // measure
    SimpleTimer tm;
    tm.start();
    for (int irep = 0; irep < count; irep++) {
        computeStiff();
    }
    tm.pause();
    // reset to zero
    for (auto &win: mWindows) {
        win->resetPointWindowsToZero();
    }
    // return total, not divided by count
    return tm.elapsedTotal();
}

//// Solid and fluid
// displ field
void Element::getDisplField(eigen::CMatXN3 &displ, bool needRTZ, int m) const {
    mWindows[m]->getDisplField(displ, mTransform, needRTZ);
}

//// Fluid only
// chi field
void Element::getChiField(eigen::CMatXN &chi, int m) const {
    mWindows[m]->getChiField(chi);
}

// pressure field
void Element::getPressureField(eigen::CMatXN &pressure, int m) const {
    mWindows[m]->getPressureField(pressure);
}

// delta field
void Element::getDeltaField(eigen::CMatXN &delta, int m) const {
    mWindows[m]->getDeltaField(delta);
}

//// Solid only
// nabla field
void Element::getNablaField(eigen::CMatXN9 &nabla, bool needRTZ, int m) const {
    mWindows[m]->getNablaField(nabla, mTransform, needRTZ);
}

// strain field
void Element::getStrainField(eigen::CMatXN6 &strain, bool needRTZ, int m) const {
    mWindows[m]->getStrainField(strain, mTransform, needRTZ);
}

// curl field
void Element::getCurlField(eigen::CMatXN3 &curl, bool needRTZ, int m) const {
    mWindows[m]->getCurlField(curl, mTransform, needRTZ);
}

// stress field
void Element::getStressField(eigen::CMatXN6 &stress, bool needRTZ, int m) const {
    mWindows[m]->getStressField(stress, mTransform, needRTZ);
}
