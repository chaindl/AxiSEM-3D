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
#include "WindowInterpolator.hpp"
#include "vicinity.hpp"
#include "geodesy.hpp"
#include "CoordTransformCartesian.hpp"
#include "CoordTransformSpherical.hpp"
#include<iostream>
using spectral::nPEM;
using interpolation::NrExtend;
// constructor
Element::Element(const int &quadTag, std::vector<std::unique_ptr<ElementWindow>> &windows,
      const std::array<std::shared_ptr<Point>, spectral::nPEM> &points,
      const std::shared_ptr<WinInt> &interpolator,
      const eigen::IMatX2 &interpTags):
mQuadTag(quadTag), mWindows(std::move(windows)), mPoints(points), 
mInterpolator(interpolator), mInterpTags(interpTags) {
    
    int maxNrOverlap = 0;
    bool needTransform = false;
    
    for (int m = 0; m < mWindows.size(); m++) {
        maxNrOverlap = std::max({maxNrOverlap, mWindows[m]->getNrOverlap(0), mWindows[m]->getNrOverlap(1)});
        if (mWindows[m]->elasticInRTZ()) needTransform = true;
    }
    expandWorkspace(maxNrOverlap);
    if (needTransform) createCoordTransform();
}

// copy constructor
Element::Element(const Element &other):
mQuadTag(other.mQuadTag), mWindows(other.copyWindows()), mPoints(other.mPoints),
mInterpolator(other.mInterpolator), mInterpTags(other.mInterpTags) {

    int maxNrOverlap = 0;
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
        mWindows[m]->transformStrainToPhysical(mTransform, sStrainWindows[m]);
    }

    overlapAndAddStrain();
    
    for (int m = 0; m < mWindows.size(); m++) {
        mWindows[m]->strainToStress(sStrainWindows[m]);
        mWindows[m]->transformStressToFourier(mTransform);
        mWindows[m]->stressToStiffness(mQuadTag);
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
        
        Eigen::Ref<eigen::RMatXX> ol_win1 = sStrainWindows[m1].block(
            mWindows[m1]->getNr() - mWindows[m1]->getNrOverlap(1), 0, mWindows[m1]->getNrOverlap(1), dim);
        Eigen::Ref<eigen::RMatXX> ol_win2 = sStrainWindows[m2].block(
            0, 0, mWindows[m2]->getNrOverlap(0), dim);
        Eigen::Ref<eigen::RMatXX> ol_win1_to_win2 = sWin1ToWin2.block(
            0, 0, mWindows[m2]->getNrOverlap(0), dim);
        Eigen::Ref<eigen::RMatXX> ol_win2_to_win1 = sWin2ToWin1.block(
            0, 0, mWindows[m1]->getNrOverlap(1), dim);    
        
        if (mInterpTags(m1, 0) >= 0) {
            mInterpolator->newDataAsBlock(mInterpTags(m2, 0), sStrainWindows[m1], 
                mWindows[m1]->getNr() - mWindows[m1]->getNrOverlap(1) - NrExtend, 
                0, mWindows[m1]->getNrOverlap(1) + NrExtend, dim);
            mInterpolator->interpolate(mInterpTags(m2, 0), ol_win1_to_win2,
                mWindows[m2]->getPhiAndFuncForInterp(0).col(0), 0, dim);
            mInterpolator->newDataAsBlock(mInterpTags(m1, 1), sStrainWindows[m2],
                0, 0, mWindows[m2]->getNrOverlap(0) + NrExtend, dim);    
            mInterpolator->interpolate(mInterpTags(m1, 1), ol_win2_to_win1,
                mWindows[m1]->getPhiAndFuncForInterp(1).col(0), 0, dim);
        } else {
            ol_win1_to_win2 = ol_win1;
            ol_win2_to_win1 = ol_win2;
        }
        ol_win1 = mWindows[m1]->getPhiAndFuncForInterp(1).col(1).asDiagonal() * ol_win1
                  + (1. - mWindows[m1]->getPhiAndFuncForInterp(1).col(1).array()).matrix().asDiagonal() * ol_win2_to_win1;
        ol_win2 = mWindows[m2]->getPhiAndFuncForInterp(0).col(1).asDiagonal() * ol_win2
                  + (1. - mWindows[m2]->getPhiAndFuncForInterp(0).col(1).array()).matrix().asDiagonal() * ol_win1_to_win2;
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
    Eigen::internal::set_is_malloc_allowed(false);
    for (int irep = 0; irep < count; irep++) {
        computeStiff();
    }
    Eigen::internal::set_is_malloc_allowed(true);
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
void Element::getNablaField(eigen::CMatXN9 &nabla, bool needRTZ, int m, int &nu) const {
    mWindows[m]->getNablaField(nabla, mTransform, needRTZ, nu);
}

// strain field
void Element::getStrainField(eigen::CMatXN6 &strain, bool needRTZ, int m, int &nu) const {
    mWindows[m]->getStrainField(strain, mTransform, needRTZ, nu);
}

// curl field
void Element::getCurlField(eigen::CMatXN3 &curl, bool needRTZ, int m, int &nu) const {
    mWindows[m]->getCurlField(curl, mTransform, needRTZ, nu);
}

// stress field
void Element::getStressField(eigen::CMatXN6 &stress, bool needRTZ, int m, int &nu) const {
    mWindows[m]->getStressField(stress, mTransform, needRTZ, nu);
}
