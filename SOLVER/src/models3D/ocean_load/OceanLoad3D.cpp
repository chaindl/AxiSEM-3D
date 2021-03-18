//
//  OceanLoad3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models

#include "OceanLoad3D.hpp"
#include "Quad.hpp"
#include "vicinity.hpp"
#include "mpi.hpp"

// apply to Quad
void OceanLoad3D::applyTo(std::vector<Quad> &quads) const {
    if (!isSuperOnly()) {
        for (Quad &quad: quads) {
            // check surface
            int surfEdge = quad.getSurfaceEdge();
            if (surfEdge == -1) {
                continue;
            }
            for (int m = 0; m < quad.getM(); m++) {
                // cardinal coordinates
                const eigen::DMatX3 &spz = computeEdgeSPZ(quad, surfEdge, m);
                // compute values
                eigen::DColX sumRD;
                bool elemInScope = getSumRhoDepth(spz, quad.getNodalSZ(), sumRD);
                // set values to quad
                if (elemInScope) {
                    setSumRhoDepthToQuad(sumRD, quad, m);
                }
            }
        }
    }
    else {
        mpi::enterInfer();
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            // step 1: gather coords on infer and send to super
            std::vector<eigen::DMatX3> spzAll;
            std::vector<eigen::DMat24> szAll;
            if (irank == mpi::rank()) {
                // gather coords
                int nwins = 0;
                for (Quad &quad: quads) {
                    if (quad.getSurfaceEdge() >= 0) {
                        nwins += quad.getM();
                    }
                }
                spzAll.reserve(nwins);
                szAll.reserve(nwins); // not ideal to duplicate sz data, but easiest for comms
                                      // alternatives: create sendVecVecEigen for spzAll or 
                                      // build and send global window tags
                                      
                if (nwins == 0) continue;
                for (Quad &quad: quads) {
                    // check surface
                    int surfEdge = quad.getSurfaceEdge();
                    if (surfEdge == -1) {
                        continue;
                    }
                    for (int m = 0; m < quad.getM(); m++) {
                        spzAll.push_back(computeEdgeSPZ(quad, surfEdge, m));
                        szAll.push_back(quad.getNodalSZ());
                    }
                }
                // send coords to super
                mpi::sendVecEigen(0, spzAll, 0);
                mpi::sendVecEigen(0, szAll, 1);
            }
            
            // step 2: compute values on super and send back to infer
            std::vector<eigen::DColX> sumRD_All;
            std::vector<eigen::IColX> elemInScopeAll;
            if (mpi::root()) {
                // recv coords from infer
                mpi::recvVecEigen(irank, spzAll, 0);
                mpi::recvVecEigen(irank, szAll, 1);
                // allocate values
                int nWin = (int)spzAll.size();
                sumRD_All.reserve(nWin);
                elemInScopeAll.push_back(eigen::IColX::Zero(nWin));
                // compute values
                for (int iw = 0; iw < nWin; iw++) {
                    eigen::DColX sumRD;
                    elemInScopeAll[0](iw) = getSumRhoDepth(spzAll[iw],
                                                           szAll[iw], sumRD);
                    sumRD_All.push_back(sumRD);
                }
                // send values to infer
                mpi::sendVecEigen(irank, sumRD_All, 0);
                mpi::sendVecEigen(irank, elemInScopeAll, 1);
            }
            
            // step 3: set values to quads on infer
            if (irank == mpi::rank()) {
                // recv values from super
                mpi::recvVecEigen(0, sumRD_All, 0);
                mpi::recvVecEigen(0, elemInScopeAll, 1);
                int iw = 0;    
                for (Quad &quad: quads) {
                    // check surface
                    int surfEdge = quad.getSurfaceEdge();
                    if (surfEdge == -1) {
                        continue;
                    }
                    // set values to quads
                    for (int m = 0; m < quad.getM(); m++) {
                        if (elemInScopeAll[0](iw)) {
                            setSumRhoDepthToQuad(sumRD_All[iw], quad, m);
                        }
                        iw++;
                    }
                }
            }
            // do irank one by one
            mpi::barrier();
        }
        mpi::enterWorld();
    }
}

// set sum(rho * depth) to quad
void OceanLoad3D::setSumRhoDepthToQuad(const eigen::DColX &sumRhoDepth,
                                       Quad &quad, const int m) const {
    // edge points
    int surfEdge = quad.getSurfaceEdge();
    const std::vector<int> &ipnts = vicinity::constants::gEdgeIPnt[surfEdge];
    const eigen::IRowN &pointNr = quad.getPointNr(m);
    // flattened to structured
    eigen::arP_DColX sumRD;
    int row = 0;
    for (int ip = 0; ip < spectral::nPED; ip++) {
        int nr = pointNr(ipnts[ip]);
        sumRD[ip] = sumRhoDepth.block(row, 0, nr, 1);
        row += nr;
    }
    // set to Quad
    quad.getOceanLoadPtr(m)->addSumRhoDepth(sumRD);
}


#include "StructuredGridO3D.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const OceanLoad3D> OceanLoad3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             const std::string &modelName, const std::string &keyInparam) {
    // short alias
    const InparamYAML &gm = inparam::gInparamModel;
    const std::string &root = keyInparam;
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "StructuredGridO3D") {
        // file name
        const std::string &fname = gm.get<std::string>(root + ":nc_data_file");
        
        ////////////// coords //////////////
        const std::string &rootc = root + ":coordinates";
        // horizontal
        bool sourceCentered = false, xy = false, ellipticity = false;
        sg_tools::inparamHorizontal(gm, rootc, modelName, className,
                                    sourceCentered, xy, ellipticity);
        // variables
        std::array<std::string, 2> crdVarNames;
        std::array<int, 2> shuffleData;
        sg_tools::inparamVarRank<2>(gm, rootc, modelName, className,
                                    crdVarNames, shuffleData);
        // units
        double lengthUnit = 1., angleUnit = 1.;
        sg_tools::inparamUnits(gm, rootc, xy, lengthUnit, angleUnit);
        
        ////////////// data //////////////
        const std::string &rootd = root + ":data_sum_rho_depth";
        const std::string &dataVarName = gm.get<std::string>(rootd + ":nc_var");
        double factor = gm.get<double>(rootd + ":factor");
        bool superOnly = gm.get<bool>(root + ":store_grid_only_on_leaders");
        
        // construct
        return std::make_shared
        <const StructuredGridO3D>(modelName, fname, crdVarNames, shuffleData,
                                  sourceCentered, xy, ellipticity,
                                  lengthUnit, angleUnit, dataVarName, factor,
                                  superOnly);
    } else {
        // other models
    }
    
    // unknown class
    return nullptr;
}
