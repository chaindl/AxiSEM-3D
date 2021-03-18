//
//  AxialBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  axial boundary condition

#include "AxialBoundary.hpp"

// point
#include "SolidPointWindow.hpp"
#include "FluidPointWindow.hpp"

// domain
#include "Messaging.hpp"

// apply axial masking
void AxialBoundary::apply() const {
    static const numerical::ComplexR czero = 0.;
    static const numerical::ComplexR cJ = {0., 1.};
    static const numerical::Real half = .5;
    
    // solid
    for (const std::shared_ptr<SolidPointWindow> &spw: mSolidPointWindows) {
        eigen::CMatX3 &stiff = spw->getFields().mStiff;
        int nu = spw->getNu_1() - 1;
        
        // alpha = 0
        stiff(0, 0) = czero;
        stiff(0, 1) = czero;
        
        // alpha > 0
        if (nu >= 1) {
            // alpha = 1
            stiff(1, 0) = (stiff(1, 0) - cJ * stiff(1, 1)) * half;
            stiff(1, 1) = cJ * stiff(1, 0);
            stiff(1, 2) = czero;
            
            // alpha >= 2
            stiff.bottomRows(nu - 1).setZero();
        }
    }
    
    // fluid
    for (const std::shared_ptr<FluidPointWindow> &fpw: mFluidPointWindows) {
        fpw->getFields().mStiff.bottomRows(fpw->getNu_1() - 1).setZero();
    }
}

// count info
void AxialBoundary::
countInfo(const Messaging &msg, int &solid, int &fluid) const {
    solid = 0;
    for (const auto &pw: mSolidPointWindows) {
        if (!msg.pointInSmallerRank(pw)) {
            solid++;
        }
    }
    fluid = 0;
    for (const auto &pw: mFluidPointWindows) {
        if (!msg.pointInSmallerRank(pw)) {
            fluid++;
        }
    }
}
