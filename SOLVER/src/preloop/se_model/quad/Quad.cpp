//
//  Quad.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  Quadrilateral
//  generator of Element in core

#include "Quad.hpp"

// mesh
#include "ExodusMesh.hpp"
#include "ABC.hpp"
#include "LocalMesh.hpp"

// nr field and windows
#include "NrField.hpp"
#include "window_tools.hpp"
#include "numerical.hpp"

// mapping
#include "MappingLinear.hpp"
#include "MappingSpherical.hpp"
#include "MappingSemiSpherical.hpp"

// GLL
#include "GLLPoint.hpp"
#include "vicinity.hpp"

// dt
#include "geodesy.hpp"
#include <iostream>
// release
#include "Domain.hpp"
#include "Element.hpp"
#include "FluidElementWindow.hpp"
#include "SolidElementWindow.hpp"

// constructor
Quad::Quad(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
           const NrField &nrField, int localTag):
mLocalTag(localTag), mGlobalTag(localMesh.mL2G_Element[localTag]),
mFluid1D(localMesh.mIsElementFluid(localTag)) {
    // model boundary
    mEdgesOnBoundary.insert({"LEFT", exodusMesh.getLeftSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"RIGHT", exodusMesh.getRightSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"BOTTOM", exodusMesh.getBottomSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"TOP", exodusMesh.getTopSide(mGlobalTag)});
    mEdgesOnBoundary.insert({"SOLID_FLUID",
        exodusMesh.getSide("solid_fluid_boundary", mGlobalTag)});
    
    // get coords form nodal (s,z)
    static eigen::DMat24 nodalSZ;
    for (int ivt = 0; ivt < 4; ivt++) {
        int inode = localMesh.mConnectivity(mLocalTag, ivt);
        nodalSZ.col(ivt) = localMesh.mNodalCoords.row(inode).transpose();
    }
    
    // mapping
    // shape type
    int gtype = localMesh.mGeometryType(mLocalTag);
    if (gtype < .5) {
        // gtype = 0.0, spherical
        mMapping = std::make_unique<MappingSpherical>(nodalSZ);
    } else if (gtype < 1.5) {
        // gtype = 1.0, linear
        mMapping = std::make_unique<MappingLinear>(nodalSZ);
    } else {
        // gtype = 2.0, semi-spherical
        mMapping = std::make_unique<MappingSemiSpherical>(nodalSZ);
    }
    
    // Nr Windows
    // gather Nr from nodes
    std::vector<std::vector<std::pair<double, double>>> nodalNr(4);
    for (int ivt = 0; ivt < 4; ivt++) {
        int inode = localMesh.mConnectivity(mLocalTag, ivt);
        nodalNr[ivt] = localMesh.mNodalNr[inode];
    }
    // combine localised Nr sections from nodes into one consistent set of sections
    // and add overlap
    std::pair<eigen::DMatX2, eigen::IMatX4> nodalNrCombined = nrField.makeElementalNrWindows(nodalNr, nodalSZ.row(0).mean());
    
    for (int m = 0; m < nodalNrCombined.first.rows(); m++) {
        int m_prev = window_tools::prevWin(m, nodalNrCombined.first.rows());
        int m_next = window_tools::nextWin(m, nodalNrCombined.first.rows());
        // interpolate pointwise nr from nodal nr
        eigen::DRow4 nodalWinNr = nodalNrCombined.second.row(m).cast<double>();
        eigen::IRowN interpPointNr = spectrals::interpolateGLL(nodalWinNr, axial()).array().round().cast<int>();
        // divide into regular windows and apply lucky number
        nrField.finalizeNrWindows(mWindows, 
            nodalNrCombined.first.row(m), interpPointNr, 
            nodalNrCombined.first(m_prev, 1), 
            nodalNrCombined.first(m_next, 0), 
            nodalSZ.row(0).mean(), nodalSZ.maxCoeff() == geodesy::getOuterRadius());
    }
    
    // can't have varying nr for points with physical storage
    // TODO: incorporate this in finaliseNrWindows() to save cost
    if (getM() > 1) {
        for (int m = 0; m < getM(); m++) {
            int maxNr = std::get<1>(*mWindows[m]).maxCoeff();
            std::get<1>(*mWindows[m]) = eigen::IRowN::Constant(1, spectral::nPEM, maxNr);
        }
    }
    
    mMaterial.reserve(getM());
    mUndulation.reserve(getM());
    mOceanLoad.reserve(getM());
    
    for (int m = 0; m < getM(); m++) {
        // window - based components
        // 1) material
        mMaterial.push_back(std::make_unique<Material>(exodusMesh, getNodalSZ(), axial()));
        // 3) undulation
        mUndulation.push_back(std::make_unique<Undulation>());
        // 4) ocean load
        mOceanLoad.push_back(std::make_unique<OceanLoad>());
    }
   
    mFluid3D.resize(getM(), mFluid1D); // placeholder for implementing discontinuous fluids
}

// setup GLL
void Quad::setupGLL(const ABC &abc, const LocalMesh &localMesh,
                    std::vector<GLLPoint> &GLLPoints) const {
    // compute mass
    static eigen::DMat2N sz;
    static std::array<eigen::DMat22, spectral::nPEM> J;
    const eigen::DRowN &ifact = computeIntegralFactor(sz, J);
    // setup tags, nr, coords
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        // target GLL
        int igll = localMesh.mElementGLL(mLocalTag, ipnt);
        // setup tags and inplane coords
        GLLPoints[igll].setup(localMesh.mL2G_GLL(igll), // global tag
                              sz.col(ipnt), // sz
                              mMapping->getMinEdgeLength() / 1000.); // tol
    }

    // setup point windows
    eigen::IMatXN winTags = eigen::IMatXN::Constant(getM(), spectral::nPEM, -1);
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        int igll = localMesh.mElementGLL(mLocalTag, ipnt);
        std::vector<eigen::DMatX2> wins_GLL(getM());
        for (int m = 0; m < getM(); m++) {
            wins_GLL[m] = eigen::DMatX2::Zero(std::get<1>(*mWindows[m])(ipnt) , 2);
            eigen::DColX phi = computeWindowPhi(m, ipnt, true);
            wins_GLL[m].col(1) = computeWindowFunction(phi, m, false);
            window_tools::wrapPhi(phi);
            wins_GLL[m].col(0) = phi;
        //     if (mEdgesOnBoundary.at("TOP") != -1) {
        //         std::vector<int> ipnts_edge = vicinity::constants::gEdgeIPnt[mEdgesOnBoundary.at("TOP")];
        //         for (int &ie: ipnts_edge) {
        //             if (ie != ipnt) continue;
        //             for (int alpha = 0; alpha < wins_GLL[m].rows(); alpha++) {
        //                 std::cout << sz(0, ipnt) << " " << wins_GLL[m](alpha, 0) / numerical::dDegree << " " << m << " " << wins_GLL[m](alpha, 1) << std::endl;
        //             }
        //         }
        //     }
        }
        winTags.col(ipnt) = GLLPoints[igll].addWindows(wins_GLL, getM() == 1);
    }
    for (int m = 0; m < getM(); m++) {
        std::get<2>(*mWindows[m]) = winTags.row(m);
    }
    
    // add mass
    for (int m = 0; m < getM(); m++) {
        const eigen::arN_DColX &J_PRT = mUndulation[m]->getMassJacobian(sz);
        const eigen::arN_DColX &mass = mMaterial[m]->getMass(ifact, J_PRT, mFluid3D[m]);
        
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].addMass(std::get<2>(*mWindows[m])(ipnt), mass[ipnt], mFluid3D[m]);
        }
    }
    
    ////////////////////////////////////
    //////////// boundaries ////////////
    ////////////////////////////////////
    
    // solid-fluid
    if (mEdgesOnBoundary.at("SOLID_FLUID") != -1 && 
        !std::all_of(mFluid3D.begin(), mFluid3D.end(), [](bool b) { return b; })) {
        
        // compute normals
        static std::vector<int> ipnts;
        eigen::DMat2P normal1D = computeNormal1D(mEdgesOnBoundary.at("SOLID_FLUID"), sz, J, ipnts);
        for (int m = 0; m < getM(); m++) {
            if (mFluid3D[m]) continue;
            // add normals to points
            for (int ip = 0; ip < spectral::nPED; ip++) {
                int igll = localMesh.mElementGLL(mLocalTag, ipnts[ip]);
                // normal must point from fluid to solid
                GLLPoints[igll].addNormalSF(std::get<2>(*mWindows[m])(ipnts[ip]),
                    -mUndulation[m]->computeNormal3D(normal1D.col(ip), sz, ipnts[ip]));
            }
        }
    }
    
    // clayton ABC
    if (abc.clayton()) {
        for (const std::string &key: abc.getBoundaryKeys()) {
            if (mEdgesOnBoundary.at(key) == -1) {
                continue;
            }
            // compute normals
            static std::vector<int> ipnts;
            eigen::DMat2P normal1D = computeNormal1D(mEdgesOnBoundary.at(key), sz, J, ipnts);
            for (int m = 0; m < getM(); m++) {
                // get properties from material
                eigen::arN_DColX rho, vp, vs;
                mMaterial[m]->getPointwiseRhoVpVs(rho, vp, vs);
                // add ABCs to points
                for (int ip = 0; ip < spectral::nPED; ip++) {
                    int ipnt = ipnts[ip];
                    int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                    GLLPoints[igll].addClaytonABC(std::get<2>(*mWindows[m])(ipnt), key, mFluid3D[m], 
                                    mUndulation[m]->computeNormal3D(normal1D.col(ip), sz, ipnts[ip]),
                                    rho[ipnt], vp[ipnt], vs[ipnt]);
                }
            }
        }
    }
    
    // sponge ABC
    if (abc.sponge()) {
        // calculate distances from edge outside window loop
        eigen::DRowN distToOuter = eigen::DRowN::Ones();
        eigen::DRowN span = eigen::DRowN::Zero();
        
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            for (const std::string &key: abc.getBoundaryKeys()) {
                /////////// pattern ///////////
                // geometry
                const auto &outerSpan = abc.getSpongeOuterSpan(key);
                double outer = std::get<0>(outerSpan);
                double span_temp = std::get<1>(outerSpan);
                // coord of me
                double coord = 0.;
                double r;
                if (geodesy::isCartesian()) {
                    coord = (key == "RIGHT" ? sz(0,ipnt) : sz(1,ipnt));
                    r = 0.;
                } else {
                    const eigen::DCol2 &rt =
                    geodesy::sz2rtheta(sz.col(ipnt).eval(), false);
                    coord = (key == "RIGHT" ? rt(1) : rt(0));
                    r = rt(0);
                }
                // get distance
                double distToOuter_temp = 1. / span_temp * (outer - coord);
                
                // store if within sponge
                if (distToOuter_temp < distToOuter(ipnt)) {
                    distToOuter(ipnt) = distToOuter_temp;
                    span(ipnt) = span_temp;
                    
                    if (!geodesy::isCartesian() && key == "RIGHT") {
                        span(ipnt) *= r;
                    }
                }
                // point is inside the inner boundary, skip;
                // point is closer to other boundary, skip;
                // there is no need to check distToOuter < 0.
            }
        }
      
        if (distToOuter.minCoeff() < 1) { // if any points inside sponge
            for (int m = 0; m < getM(); m++) {
                // get material
                eigen::arN_DColX rho, vp, vs;
                mMaterial[m]->getPointwiseRhoVpVs(rho, vp, vs);
                
                // loop over gll
                for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                    if (distToOuter(ipnt) == 1) continue;
                    // regularize 1D/3D
                    op1D_3D::regularize1D<eigen::DColX>({
                        std::ref(rho[ipnt]),
                        std::ref(vp[ipnt]),
                        std::ref(vs[ipnt])});
                
                    /////////// pattern ///////////          
                    static const double piHalf = numerical::dPi / 2.;
                    double pattern = pow(cos(piHalf * distToOuter(ipnt)), 2.);
                    
                    /////////// gamma ///////////
                    // for theta, change span to arc-length
                    eigen::DColX gamma = eigen::DColX::Zero(std::get<1>(*mWindows[m])(ipnt), 1);

                    if (!mFluid3D[m]) {
                        gamma = abc.getU0Solid(std::abs(span(ipnt)),
                                                  vp[ipnt], vs[ipnt], rho[ipnt]);
                    } else {
                        gamma = abc.getU0Fluid(std::abs(span(ipnt)),
                                                  vp[ipnt], rho[ipnt]);
                    }
                    gamma *= pattern;
                    
                    // release
                    if (gamma.norm() > numerical::dEpsilon) {
                        int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                        GLLPoints[igll].addGamma(std::get<2>(*mWindows[m])(ipnt), gamma);
                    }
                }
            }
        }
    }
    
    // ocean load
    // get data from OceanLoad
    if (mEdgesOnBoundary.at("TOP") != -1 && 
        std::any_of(mOceanLoad.begin(), mOceanLoad.end(), [](const std::unique_ptr<OceanLoad> &b) {return *b;})) {
        // compute normals
        static std::vector<int> ipnts;
        eigen::DMat2P normal1D = computeNormal1D(mEdgesOnBoundary.at("TOP"), sz, J, ipnts);
        // add ocean loads to points
        for (int m = 0; m < getM(); m++) {
            if (!mOceanLoad[m]) continue;
            const eigen::arP_DColX &sumRhoDepth = mOceanLoad[m]->getPointwise();
            for (int ip = 0; ip < spectral::nPED; ip++) {
                int igll = localMesh.mElementGLL(mLocalTag, ipnts[ip]);
                GLLPoints[igll].addOceanLoad(std::get<2>(*mWindows[m])(ipnts[ip]), 
                    mUndulation[m]->computeNormal3D(normal1D.col(ip), sz, ipnts[ip]), 
                    sumRhoDepth[ip]);
            }
        }
    }
    
    // axial
    if (axial()) {
        const std::vector<int> &ipnts =
        vicinity::constants::gEdgeIPnt[mEdgesOnBoundary.at("LEFT")];
        for (int ipnt: ipnts) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].setAxial();
        }
    }
    
    // surface
    if (mEdgesOnBoundary.at("TOP") != -1) {
        const std::vector<int> &ipnts =
        vicinity::constants::gEdgeIPnt[mEdgesOnBoundary.at("TOP")];
        for (int ipnt: ipnts) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].setSurface();
        }
    }
    
    for (int m = 0; m < getM(); m++) {
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            int igll = localMesh.mElementGLL(mLocalTag, ipnt);
            GLLPoints[igll].exchange3DInfoBetweenWindows();
        }
    }
}

// compute dt
double Quad::computeDt(double courant, const ABC &abc) const {
    // 1D coords in reference configuration
    const eigen::DMat2N &szRef = getPointSZ();
    
    eigen::DRowN sint, cost;
    if (!geodesy::isCartesian()) {
        const eigen::DMat2N &rt = geodesy::sz2rtheta(szRef, false);
        sint = rt.row(1).array().sin();
        cost = rt.row(1).array().cos();
    }
    
    double dt = std::numeric_limits<double>::max();
    for (int m = 0; m < getM(); m++) {
        // get vp
        const eigen::DColX &vmaxNr =
        mMaterial[m]->getMaxVelocity().rowwise().maxCoeff();
        
        // compute coords on slices
        std::vector<eigen::DMat2N> szNr;
        const eigen::DMatXN &dZ = mUndulation[m]->getElemental();
        // both 1D and 3D undulation
        if (geodesy::isCartesian()) {
            for (int nr = 0; nr < dZ.rows(); nr++) {
                eigen::DMat2N sz = szRef;
                sz.row(1) += dZ.row(nr);
                szNr.push_back(sz);
            }
        } else {
            for (int nr = 0; nr < dZ.rows(); nr++) {
                eigen::DMat2N sz = szRef;
                sz.row(0) += dZ.row(nr).cwiseProduct(sint);
                sz.row(1) += dZ.row(nr).cwiseProduct(cost);
                szNr.push_back(sz);
            }
        }
        
        // compute hmin and dt slice-wise
        int nrMax = std::max((int)vmaxNr.rows(), (int)szNr.size());
        for (int nr = 0; nr < nrMax; nr++) {
            // vmax
            double vmax = (vmaxNr.rows() == 1) ? vmaxNr(0) : vmaxNr(nr);
            
            // hmin
            const eigen::DMat2N &sz = (szNr.size() == 1) ? szNr[0] : szNr[nr];
            double hmin = std::numeric_limits<double>::max();
            for (int ipnt0 = 0; ipnt0 < spectral::nPEM - 1; ipnt0++) {
                int ipol0 = ipnt0 / spectral::nPED;
                int jpol0 = ipnt0 % spectral::nPED;
                for (int ipnt1 = ipnt0 + 1; ipnt1 < spectral::nPEM; ipnt1++) {
                    int ipol1 = ipnt1 / spectral::nPED;
                    int jpol1 = ipnt1 % spectral::nPED;
                    // only consider neighbouring points
                    if (std::abs(ipol1 - ipol0) <= 1 &&
                        std::abs(jpol1 - jpol0) <= 1) {
                        hmin = std::min(hmin, (sz.col(ipnt1) -
                                               sz.col(ipnt0)).norm());
                    }
                }
            }
            // dt
            dt = std::min(dt, hmin / vmax);
        }
    }
    
    // solid-fluid and clayton BCs are numerically sensitive
    double factorDtForBC = 1.;
    // solid-fluid
    if (mEdgesOnBoundary.at("SOLID_FLUID") != -1) {
        factorDtForBC = .9;
    }
    // clayton ABC
    if (abc.clayton()) {
        for (const std::string &key: abc.getBoundaryKeys()) {
            if (mEdgesOnBoundary.at(key) != -1) {
                factorDtForBC = .5;
                break;
            }
        }
    }
    // decrease DT
    dt *= factorDtForBC;
    
    // courant
    return dt * courant;
}

// release to domain
void Quad::release(const LocalMesh &localMesh,
                   const std::vector<GLLPoint> &GLLPoints,
                   const std::unique_ptr<const AttBuilder> &attBuilder,
                   const std::shared_ptr<WindowInterpolator<numerical::Real>> &interpolator,
                   Domain &domain) {

    static eigen::DMat2N sz;
    static eigen::DMatPP_RM ifPP;
    
    std::array<std::shared_ptr<Point>, spectral::nPEM> points;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        int igll = localMesh.mElementGLL(mLocalTag, ipnt);
        points[ipnt] = GLLPoints[igll].getPoint();
    }
    
    double tol = 2 * numerical::dPi * numerical::dEpsilon / getM();
    
    // prepare overlap comms and 
    // construct splines if needed
    std::vector<eigen::ICol2> nr_ol(getM());
    std::vector<eigen::DCol2> width_interp(getM());
    std::vector<eigen::DCol2> phi1_interp(getM());
    eigen::IMatX2 interpTags = eigen::IMatX2::Constant(getM(), 2, -1);
    if (getM() > 1) {
        for (int m = 0 ; m < getM(); m++) {
            nr_ol[m] = eigen::ICol2::Zero(2, 1);
            int m_prev = window_tools::prevWin(m, getM());
            int m_next = window_tools::nextWin(m, getM());
          
            eigen::DRow4 shape = window_tools::unwrapPhi(std::get<0>(*mWindows[m]), tol);
            
            eigen::DColX winPhi = eigen::DColX::LinSpaced(std::get<1>(*mWindows[m]).maxCoeff(), shape(0), shape(3));
            double dphi_win = winPhi(1) - winPhi(0);
            
            // handle left overlap
            if (std::get<3>(*mWindows[m_prev])) {
                eigen::DRow4 shape_prev = window_tools::unwrapPhi(std::get<0>(*mWindows[m_prev]), tol);
                double dphi_prev_win = (shape_prev(3) - shape_prev(0)) / (std::get<1>(*mWindows[m_prev]).maxCoeff() - 1);
                
                while (winPhi(nr_ol[m](0)) < shape(1) - tol) nr_ol[m](0)++;
                bool aligned = (std::abs(winPhi(nr_ol[m](0)) - shape(1)) < tol) 
                               && (std::abs(dphi_prev_win - dphi_win) < tol);
                
                if (!aligned) {
                    phi1_interp[m](0) = window_tools::setBounds2Pi(winPhi(0));
                    width_interp[m](0) = winPhi(nr_ol[m](0) + interpolation::NrExtend - 1) - winPhi(0);
                    interpTags(m, 0) = interpolator->addSplineFitting(eigen::RColX::LinSpaced(nr_ol[m](0) + interpolation::NrExtend, 0., 1.), 0);
                }
            }
            
            // handle right overlap
            if (std::get<3>(*mWindows[m])) {
                eigen::DRow4 shape_next = window_tools::unwrapPhi(std::get<0>(*mWindows[m_next]), tol);
                double dphi_next_win = (shape_next(3) - shape_next(0)) / (std::get<1>(*mWindows[m_next]).maxCoeff() - 1);
                
                while ((winPhi(winPhi.rows() - nr_ol[m](1) - 1) > shape(2) + tol)) nr_ol[m](1)++;
                bool aligned = (std::abs(winPhi(winPhi.rows() - nr_ol[m](1) - 1) - shape(2)) < tol) 
                               && (std::abs(dphi_next_win - dphi_win) < tol);
                

                if (!aligned) {
                    phi1_interp[m](1) = window_tools::setBounds2Pi(winPhi(winPhi.rows() - nr_ol[m](1) - interpolation::NrExtend - 1));
                    width_interp[m](1) = shape(3) - winPhi(winPhi.rows() - nr_ol[m](1) - interpolation::NrExtend - 1);
                    interpTags(m, 1) = interpolator->addSplineFitting(eigen::RColX::LinSpaced(nr_ol[m](1) + interpolation::NrExtend, 0., 1.), 0);
                }
            }
        }
    }
  
    std::vector<std::unique_ptr<ElementWindow>> eleWins;
    for (int m = 0 ; m < getM(); m++) {
      
        std::unique_ptr<const PRT> prt = mUndulation[m]->createPRT(sz);
        std::array<eigen::RMatX2, 2> overlap;
        
        if (getM() > 1) {
            double dphi = window_tools::setBounds2Pi(std::get<0>(*mWindows[m])(3) - std::get<0>(*mWindows[m])(0));
            dphi /= std::get<1>(*mWindows[m]).maxCoeff();
            
            // create gradient with scaling
            double winFrac = 2 * numerical::dPi / (std::get<1>(*mWindows[m]).maxCoeff() + 2 * interpolation::NrExtend * dphi);
            std::unique_ptr<const GradientQuadrature<numerical::Real>> grad 
            = createGradient<numerical::Real>(sz, ifPP, winFrac);
            
            eigen::DMatX2 ol_left, ol_right;
            int m_prev = window_tools::prevWin(m, getM());
            int m_next = window_tools::nextWin(m, getM());
            
            int dims = mFluid3D[m] ? 3 * spectral::nPEM : 9 * spectral::nPEM;
            if (std::get<3>(*mWindows[m_prev])) {
                eigen::DColX phi_left = eigen::DColX::LinSpaced(nr_ol[m](0), 0., nr_ol[m](0) - 1).array() * dphi + std::get<0>(*mWindows[m])(0);
              
                ol_left = eigen::DMatX2::Zero(nr_ol[m](0), 2);
                ol_left.col(0) = phi_left;
                window_tools::wrapPhi(phi_left);
                ol_left.col(1) = computeWindowFunction(phi_left, m, false);
                
                if (interpTags(m, 0) >= 0) {
                    ol_left.col(0) = ol_left.col(0).array() - phi1_interp[m_prev](1);
                    if (ol_left(0, 0) < 0) ol_left.col(0) = ol_left.col(0).array() + 2 * numerical::dPi;
                    ol_left.col(0) = ol_left.col(0).array() / width_interp[m_prev](1);
                    interpolator->expandWorkspace(interpTags(m_prev, 1), dims);
                }
            }
            
            if (std::get<3>(*mWindows[m])) {
                eigen::DColX phi_right = eigen::DColX::LinSpaced(nr_ol[m](1), -(nr_ol[m](1) - 1), 0.).array() * dphi + std::get<0>(*mWindows[m])(3);
                if (phi_right(0) < 0) phi_right = phi_right.array() + 2 * numerical::dPi;
                ol_right = eigen::DMatX2::Zero(nr_ol[m](1), 2);
                ol_right.col(0) = phi_right;
                window_tools::wrapPhi(phi_right);
                ol_right.col(1) = computeWindowFunction(phi_right, m, false);
                
                if (interpTags(m, 1) >= 0) {
                    ol_right.col(0) = ol_right.col(0).array() - phi1_interp[m_next](0);
                    if (ol_right(0, 0) < 0) ol_right.col(0) = ol_right.col(0).array() + 2 * numerical::dPi;
                    ol_right.col(0) = ol_right.col(0).array() / width_interp[m_next](0);
                    interpolator->expandWorkspace(interpTags(m_next, 0), dims);
                }
            }
            
            std::array<eigen::RMatX2, 2> overlap = {ol_left.cast<numerical::Real>(), ol_right.cast<numerical::Real>()};
            if (mFluid3D[m]) {
                std::unique_ptr<const Acoustic> acoustic = mMaterial[m]->createAcoustic();
                std::array<std::shared_ptr<FluidRCPointWindow<1, numerical::Real>>, spectral::nPEM> wins;
                for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                    int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                    wins[ipnt] = GLLPoints[igll].getFluidPointWindowR(std::get<2>(*mWindows[m])(ipnt));
                }
                eleWins.push_back(std::make_unique<FluidElementWindow<FluidRCPointWindow<1, numerical::Real>>>(grad, prt, acoustic, wins, overlap, interpolator));
            } else { 
                std::unique_ptr<const Elastic> elastic =
                mMaterial[m]->createElastic(attBuilder, computeWeightsCG4(ifPP));
                std::array<std::shared_ptr<SolidRCPointWindow<3, numerical::Real>>, spectral::nPEM> wins;
                for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                    int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                    wins[ipnt] = GLLPoints[igll].getSolidPointWindowR(std::get<2>(*mWindows[m])(ipnt));
                }
                eleWins.push_back(std::make_unique<SolidElementWindow<SolidRCPointWindow<3, numerical::Real>>>(grad, prt, elastic, wins, overlap, interpolator));
            }
            
        } else { // have point windows with storage in Fourier
            std::unique_ptr<const GradientQuadrature<numerical::Real>> grad 
            = createGradient<numerical::Real>(sz, ifPP, 1.);
            
            if (mFluid3D[m]) {
                std::unique_ptr<const Acoustic> acoustic = mMaterial[m]->createAcoustic();
                std::array<std::shared_ptr<FluidRCPointWindow<1, numerical::ComplexR>>, spectral::nPEM> wins;
                for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                    int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                    wins[ipnt] = GLLPoints[igll].getFluidPointWindowC(std::get<2>(*mWindows[m])(ipnt));
                }
                eleWins.push_back(std::make_unique<FluidElementWindow<FluidRCPointWindow<1, numerical::ComplexR>>>(grad, prt, acoustic, wins, overlap, nullptr));
            } else { 
                std::unique_ptr<const Elastic> elastic =
                mMaterial[m]->createElastic(attBuilder, computeWeightsCG4(ifPP));
                std::array<std::shared_ptr<SolidRCPointWindow<3, numerical::ComplexR>>, spectral::nPEM> wins;  
                for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                    int igll = localMesh.mElementGLL(mLocalTag, ipnt);
                    wins[ipnt] = GLLPoints[igll].getSolidPointWindowC(std::get<2>(*mWindows[m])(ipnt));
                }
                eleWins.push_back(std::make_unique<SolidElementWindow<SolidRCPointWindow<3, numerical::ComplexR>>>(grad, prt, elastic, wins, overlap, nullptr));
            }
        }
    }
    
    mElement = std::make_shared<Element>(mGlobalTag, eleWins, points, interpolator, interpTags);
    domain.addElement(mElement);
    
    // free dummy memory
    mMaterial.clear();
    mUndulation.clear();
    mOceanLoad.clear();
}

//////////////////////// integral ////////////////////////
// compute integral factor
eigen::DRowN Quad::
computeIntegralFactor(eigen::DMat2N &sz, std::array<eigen::DMat22,
                      spectral::nPEM> &J) const {
    // xieta and weights
    const eigen::DMat2N &xieta = spectrals::getXiEtaElement(axial());
    const eigen::DRowN &weights = spectrals::getWeightsElement(axial());
    
    // compute integral factor
    eigen::DRowN ifact;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        // element Jacobian
        J[ipnt] = mMapping->jacobian(xieta.col(ipnt));
        // integral weights and Jacobian
        ifact(ipnt) = weights(ipnt) * J[ipnt].determinant();
        // s in integral
        sz.col(ipnt) = mMapping->mapping(xieta.col(ipnt));
        double s = sz(0, ipnt);
        if (axial()) {
            // axial element
            if (ipnt / spectral::nPED == 0) { // ipol == 0
                // axial point, L'Hospital's Rule
                ifact(ipnt) *= J[ipnt](0, 0);
            } else {
                // non-axial point
                ifact(ipnt) *= s / (1. + xieta(0, ipnt));
            }
        } else {
            // non-axial element
            ifact(ipnt) *= s;
        }
    }
    return ifact;
}

// get normal
eigen::DMat2P Quad::computeNormal1D(int edge, const eigen::DMat2N &sz,
                         const std::array<eigen::DMat22, spectral::nPEM> &J,
                         std::vector<int> &ipnts) const {
    // xieta
    const eigen::DMat2N &xieta = spectrals::getXiEtaElement(axial());
    
    // points on edge
    ipnts = vicinity::constants::gEdgeIPnt[edge];
    
    eigen::DMat2P n1D;
    
    // normal
    for (int ip = 0; ip < spectral::nPED; ip++) {
        int ipnt = ipnts[ip];
        // 1D normal
        n1D.col(ip) = mMapping->normal(edge, J[ipnt]);
        // integral weights
        int ipol = ipnt / spectral::nPED;
        int jpol = ipnt % spectral::nPED;
        if (axial()) {
            // edge must be even
            if (edge % 2 != 0) {
                throw std::runtime_error("Quad::computeNormal || Impossible.");
            }
            n1D.col(ip) *= spectrals::gWeightsGLJ(ipol);
        } else {
            n1D.col(ip) *= spectrals::gWeightsGLL(edge % 2 == 0 ? ipol : jpol);
        }
        // s in integral
        double s = sz(0, ipnt);
        if (axial()) {
            // axial element
            if (ipol == 0) {
                // axial point, L'Hospital's Rule
                n1D.col(ip) *= J[ipnt](0, 0);
            } else {
                // non-axial point
                n1D.col(ip) *= s / (1. + xieta(0, ipnt));
            }
        } else {
            // non-axial element
            n1D.col(ip) *= s;
        }
    }
    return n1D;
}

// weights for CG4 attenuation
eigen::DRow4 Quad::computeWeightsCG4(const eigen::DMatPP_RM &ifPP) const {
    if (spectral::nPol != 4) {
        return eigen::DRow4::Zero();
    }
    // weights on CG4 points (marked O)
    // x x x x x
    // x O x O x
    // x x x x x
    // x O x O x
    // x x x x x
    eigen::DRow4 wcg4;
    wcg4(0) = (ifPP(0, 0) + ifPP(0, 1) + ifPP(1, 0) + ifPP(1, 1) +
               0.50 * (ifPP(0, 2) + ifPP(1, 2) + ifPP(2, 0) + ifPP(2, 1)) +
               0.25 * ifPP(2, 2)) / ifPP(1, 1);
    wcg4(1) = (ifPP(0, 3) + ifPP(0, 4) + ifPP(1, 3) + ifPP(1, 4) +
               0.50 * (ifPP(0, 2) + ifPP(1, 2) + ifPP(2, 3) + ifPP(2, 4)) +
               0.25 * ifPP(2, 2)) / ifPP(1, 3);
    wcg4(2) = (ifPP(3, 0) + ifPP(3, 1) + ifPP(4, 0) + ifPP(4, 1) +
               0.50 * (ifPP(2, 0) + ifPP(2, 1) + ifPP(3, 2) + ifPP(4, 2)) +
               0.25 * ifPP(2, 2)) / ifPP(3, 1);
    wcg4(3) = (ifPP(3, 3) + ifPP(3, 4) + ifPP(4, 3) + ifPP(4, 4) +
               0.50 * (ifPP(2, 3) + ifPP(2, 4) + ifPP(3, 2) + ifPP(4, 2)) +
               0.25 * ifPP(2, 2)) / ifPP(3, 3);
    return wcg4;
}

eigen::DColX Quad::computeWindowPhi(int m, int ipnt, bool keep_unwrapped) const {
    eigen::DCol2 phi_lims;
    if (std::get<3>(*mWindows[m])) { // has overlap
        phi_lims(0) = std::get<0>(*mWindows[m])(0);
        phi_lims(1) = std::get<0>(*mWindows[m])(3);
        window_tools::unwrapPhi(phi_lims);
    } else {
        double dphi = std::get<0>(*mWindows[m])(3) - std::get<0>(*mWindows[m])(0);
        if (dphi <= 0) dphi += 2 * numerical::dPi;
        
        phi_lims(0) = std::get<0>(*mWindows[m])(0);
        phi_lims(1) = std::get<0>(*mWindows[m])(0) + (1. - 1. / std::get<1>(*mWindows[m])(ipnt)) * dphi;
    }
    
    eigen::DColX phi = eigen::DColX::LinSpaced(std::get<1>(*mWindows[m])(ipnt), phi_lims(0), phi_lims(1));
    
    if (!keep_unwrapped) window_tools::wrapPhi(phi);
    return phi;
}    

eigen::DColX Quad::computeWindowFunction(eigen::DColX phi, int m, bool relative_phi) const {
    double tol = 2 * numerical::dPi * numerical::dEpsilon;
    eigen::DRow4 winShape = window_tools::unwrapPhi(std::get<0>(*mWindows[m]), tol);
    
    eigen::DColX frac = eigen::DColX::Constant(phi.rows(), 1);
    if (!relative_phi) { // might need unwrapping + check bounds
        window_tools::unwrapPhi(phi);
        if (phi(0) < winShape(0) - tol) phi.array() += 2 * numerical::dPi;
        if (phi(0) < winShape(0) - tol || winShape(3) < phi(phi.rows() - 1) - tol) {
            throw std::runtime_error("Quad::computeWindowFunction || Angle outside window bounds.");
        }
    }
    if (winShape(0) == winShape(1) && winShape(2) == winShape(3)) { // no overlap
        return frac;
    }
    
    if (relative_phi) {
        winShape.array() -= winShape(0);
        winShape.array() *= 2 * numerical::dPi / winShape(3);
    }
    
    int i1 = -1, i2 = -1;
    if (phi(0) < winShape(1) + tol) { // left edge
        for (int i = 0; i < phi.size(); i++) {
            i1 = i;
            if (phi(i) >= winShape(1)) break;
        }
    }
    
    if (phi(phi.size() - 1) > winShape(2) - tol) { // right edge
        for (int i = std::max(i1, 0); i < phi.size(); i++) {
            i2 = i;
            if (phi(i) >= winShape(2)) break;
        }
    }
    
    if (i1 >= 0 && winShape(0) != winShape(1)) {
        frac.topRows(i1) = window_tools::getWindowSlope((phi.topRows(i1).array() - winShape(0)) / (winShape(1) - winShape(0)));
    }
    
    if (i1 >= 0 && i2 >= 0) frac.segment(i1, i2 - i1).array() = 1.;
    
    if (i2 >= 0 && winShape(2) != winShape(3)) {
        frac.bottomRows(phi.size() - i2) = 1. - window_tools::getWindowSlope((winShape(2) - phi.bottomRows(phi.size() - i2).array()) / (winShape(2) - winShape(3))).array();
    }

    return frac;
}

eigen::DMatXX Quad::computeRelativeWindowPhis(const std::vector<double> &phis, std::vector<double> &winFunc) const {
    eigen::DMatXX out = eigen::DMatXX::Constant(phis.size(), mWindows.size(), -1);
    
    for (int m = 0; m < mWindows.size(); m++) {
        double tol = 2 * numerical::dPi * numerical::dEpsilon;
        eigen::DRow4 winShape = window_tools::unwrapPhi(std::get<0>(*mWindows[m]), tol);
        
        for (int i = 0; i < phis.size(); i++) {
            double phi = phis[i];
            if (std::get<0>(*mWindows[m])(0) > std::get<0>(*mWindows[m])(3) && phi < std::get<0>(*mWindows[m])(3)) phi += 2 * numerical::dPi;
              
            if (winShape(0) <= phi && phi <= winShape(3)) {
                out(i, m) = phi - winShape(0);
                out(i, m) *= (winShape(3) - winShape(0));
                
                if (phi <= winShape(1)) {
                    winFunc.push_back(window_tools::getWindowSlope((phi - winShape(0)) / (winShape(1) - winShape(0))));
                } else if (winShape(2) <= phi) {
                    winFunc.push_back(1. - window_tools::getWindowSlope((phi - winShape(2)) / (winShape(3) - winShape(2))));
                } else {
                    winFunc.push_back(1.);
                }
            }
        }
    }
    return out;
}
