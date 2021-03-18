//
//  Mass3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  3D mass

#include "Mass3D.hpp"

// check compatibility
void Mass3D::checkCompatibility(int nr, bool solid) const {
    // check size
    if (mInvMass.rows() != nr) {
        throw std::runtime_error("Mass3D::checkCompatibility || "
                                 "Incompatible sizes.");
    }
}

// compute accel in-place for fluid
void Mass3D::computeAccel(eigen::RColX &stiffR1) const {
    // constants
    int nr = (int)mInvMass.rows();
    
    // divide by mass in cardinal space
    stiffR1.topRows(nr).array() *= mInvMass.array();
}

// compute accel in-place for solid
void Mass3D::computeAccel(eigen::RMatX3 &stiffR3) const {
    // constants
    int nr = (int)mInvMass.rows();

    // divide by mass in cardinal space
    stiffR3.topRows(nr).applyOnTheLeft(mInvMass.asDiagonal());
}
