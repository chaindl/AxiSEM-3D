//
//  Geometric3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D geometric models

#include "Geometric3D.hpp"
#include "Quad.hpp"
#include "mpi.hpp"

// apply to Quad
void Geometric3D::applyTo(std::vector<Quad> &quads) const {
    if (!isSuperOnly()) {
        for (Quad &quad: quads) {
            for (int m = 0; m < quad.getM(); m++) {
                // cardinal coordinates
                const eigen::DMatX3 &spz = computeElemSPZ(quad, m);
                // compute values
                eigen::DColX und;
                bool elemInScope = getUndulation(spz, quad.getNodalSZ(), und);
                // set values to quad
                if (elemInScope) {
                    setUndulationToQuad(und, quad, m);
                }
            }
        }
    } else {
        mpi::enterInfer();
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            // step 1: gather coords on infer and send to super
            std::vector<eigen::DMatX3> spzAll;
            std::vector<eigen::DMat24> szAll;
            if (irank == mpi::rank()) {
                // gather coords
                int nwins = 0;
                for (Quad &quad: quads) {
                    nwins += quad.getM();
                }
                spzAll.reserve(nwins);
                szAll.reserve(nwins); // not ideal to duplicate sz data, but easiest for comms
                                      // alternatives: create sendVecVecEigen for spzAll or 
                                      // build and send global window tags
                for (Quad &quad: quads) {
                    for (int m = 0; m < quad.getM(); m++) {
                        spzAll.push_back(computeElemSPZ(quad, m));
                        szAll.push_back(quad.getNodalSZ());
                    }
                }
                // send coords to super
                mpi::sendVecEigen(0, spzAll, 0);
                mpi::sendVecEigen(0, szAll, 1);
            }
            
            // step 2: compute values on super and send back to infer
            std::vector<eigen::DColX> undAll;
            std::vector<eigen::IColX> elemInScopeAll;
            if (mpi::root()) {
                // recv coords from infer
                mpi::recvVecEigen(irank, spzAll, 0);
                mpi::recvVecEigen(irank, szAll, 1);
                // allocate values
                int nWin = (int)spzAll.size();
                undAll.reserve(nWin);
                elemInScopeAll.push_back(eigen::IColX::Zero(nWin));
                // compute values
                for (int iw = 0; iw < nWin; iw++) {
                    eigen::DColX und;
                    elemInScopeAll[0](iw) = getUndulation(spzAll[iw], szAll[iw],
                                                          und);
                    undAll.push_back(und);
                }
                // send values to infer
                mpi::sendVecEigen(irank, undAll, 0);
                mpi::sendVecEigen(irank, elemInScopeAll, 1);
            }
            
            // step 3: set values to quads on infer
            if (irank == mpi::rank()) {
                // recv values from super
                mpi::recvVecEigen(0, undAll, 0);
                mpi::recvVecEigen(0, elemInScopeAll, 1);
                // set values to quads
                int iw = 0;
                for (Quad &quad: quads) {
                    for (int m = 0; m < quad.getM(); m++) {
                        if (elemInScopeAll[0](iw)) {
                            setUndulationToQuad(undAll[iw], quads[iw], m);
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

// set undulation to quad
void Geometric3D::setUndulationToQuad(const eigen::DColX &undulation,
                                      Quad &quad, const int m) const {
    // flattened to structured
    const eigen::IRowN &pointNr = quad.getPointNr(m);
    eigen::arN_DColX undArr;
    int row = 0;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        int nr = pointNr(ipnt);
        undArr[ipnt] = undulation.block(row, 0, nr, 1);
        row += nr;
    }
    // set to Quad
    quad.getUndulationPtr(m)->addUndulation(undArr);
}


#include "StructuredGridG3D.hpp"
#include "Ellipticity.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const Geometric3D> Geometric3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             const std::string &modelName, const std::string &keyInparam) {
    // short alias
    const InparamYAML &gm = inparam::gInparamModel;
    const std::string &root = keyInparam;
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "StructuredGridG3D") {
        // file name
        const std::string &fname = gm.get<std::string>(root + ":nc_data_file");
        
        ////////////// coords //////////////
        const std::string &rootc = root + ":coordinates";
        // horizontal
        bool sourceCentered = false, xy = false, ellipticity = false;
        sg_tools::inparamHorizontal(gm, rootc, modelName, className,
                                    sourceCentered, xy, ellipticity);
        // vertical
        bool useDepth = false, depthSolid = false;
        sg_tools::inparamVertical(gm, rootc, modelName, className,
                                  useDepth, depthSolid);
        // variables
        std::array<std::string, 2> crdVarNames;
        std::array<int, 2> shuffleData;
        sg_tools::inparamVarRank<2>(gm, rootc, modelName, className,
                                    crdVarNames, shuffleData);
        // units
        double lengthUnit = 1., angleUnit = 1.;
        sg_tools::inparamUnits(gm, rootc, xy, lengthUnit, angleUnit);
        
        ////////////// undulated range //////////////
        const std::string &rootr = root + ":undulation_range";
        double interface = gm.get<double>(rootr + ":interface");
        double min = gm.get<double>(rootr + ":min_max:[0]");
        double max = gm.get<double>(rootr + ":min_max:[1]");
        if (interface >= max || interface <= min) {
            throw std::runtime_error
            ("Geometric3D::buildInparam || undulation_range:interface "
             "must lie within undulation_range:min_max."
             " || Model name: " + modelName +
             " || Class name: " + className);
        }
        
        ////////////// data //////////////
        const std::string &rootu = root + ":undulation_data";
        const std::string &dataVarName = gm.get<std::string>(rootu + ":nc_var");
        double factor = gm.get<double>(rootu + ":factor");
        bool superOnly = gm.get<bool>(root + ":store_grid_only_on_leaders");
        
        // construct
        return std::make_shared
        <const StructuredGridG3D>(modelName, fname, crdVarNames, shuffleData,
                                  sourceCentered, xy, ellipticity,
                                  useDepth, depthSolid,
                                  interface, min, max, lengthUnit, angleUnit,
                                  dataVarName, factor, superOnly);
    } else if (className == "Ellipticity") {
        // ellipticity can be added only once
        static bool ellipticityAdded = false;
        if (ellipticityAdded) {
            throw
            std::runtime_error("Geometric3D::buildInparam || Ellipticity "
                               "model cannot be added more than once.");
        }
        if (io::gVerboseWarnings) {
            if (geodesy::isCartesian()) {
                io::cout <<
                bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                 "model will be ignored for a Cartesian mesh.");
                
            }
            if (geodesy::getOuterFlattening() < numerical::dEpsilon) {
                io::cout <<
                bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                 "model will be ignored with zero flattening.");
                
            }
        }
        ellipticityAdded = true;
        return std::make_shared<const Ellipticity>(modelName);
    }
    
    // unknown class
    return nullptr;
}
