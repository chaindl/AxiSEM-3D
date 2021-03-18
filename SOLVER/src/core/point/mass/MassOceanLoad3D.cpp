//
//  MassOceanLoad3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  3D mass with ocean load
//  used on solid-fluid boundary when the fluid is modelled as load
//  a = F.n n / (m + m0) + (F - F.n n) / m
//    = [m F.n n + (m + m0) (F - F.n n)] / [m (m + m0)]
//    = [(m + m0) F - m0 F.n n] / [m (m + m0)]
//    = F / m - m0 / [m (m + m0)] F.n n
//    = im F - F.k k, with k = sqrt(m0 / [m (m + m0)]) n
//  a -> acceleration
//  F -> force
//  n -> unit normal of surface
//  m -> mass of solid
//  m0 -> mass of water column above

#include "MassOceanLoad3D.hpp"

// constructor
MassOceanLoad3D::MassOceanLoad3D(const eigen::DColX &mass,
                                 const eigen::DColX &massOcean,
                                 const eigen::DMatX3 &unitNormal):
mIM(mass.cwiseInverse().cast<numerical::Real>()),
mK((massOcean.cwiseQuotient(mass.cwiseProduct(mass + massOcean))
    .cwiseSqrt().asDiagonal() * unitNormal).cast<numerical::Real>()) {
    // nothing
    // im = 1 / m
    // k = sqrt(m0 / [m (m + m0)]) n
}

// check compatibility
void MassOceanLoad3D::checkCompatibility(int nr, bool solid) const {
    // must on solid
    if (!solid) {
        throw std::runtime_error("MassOceanLoad3D::checkCompatibility || "
                                 "Incompatible types: "
                                 "ocean load on fluid point.");
    }
    
    // check size
    if (mIM.rows() != nr) {
        throw std::runtime_error("MassOceanLoad3D::checkCompatibility || "
                                 "Incompatible sizes.");
    }
    
    // expand workspace if needed
    if (sF.rows() < nr) {
        sF.resize(nr, 3);
    }
}

// compute accel in-place for solid
void MassOceanLoad3D::computeAccel(eigen::RMatX3 &stiff3) const {
    // constants
    int nr = (int)mIM.rows();
    
    sF.topRows(nr) = stiff3.topRows(nr);
    
    // a = im F
    stiff3.topRows(nr) = mIM.asDiagonal() * sF.topRows(nr);
    
    // a -= F.k k
    stiff3.topRows(nr) -= (sF.topRows(nr).cwiseProduct(mK)
                       .rowwise().sum().asDiagonal() * mK);
}
