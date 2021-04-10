//
//  FluidSurfaceBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  stress-free boundary condition on fluid

#include "FluidSurfaceBoundary.hpp"

// point
#include "FluidPointWindow.hpp"

// domain
#include "Messaging.hpp"

// apply stress-free boundary condition on fluid
void FluidSurfaceBoundary::apply() const {
    // pressure ≡ 0 or accel ≡ 0
    // so, veloc = disp = everything ≡ 0
    for (const std::shared_ptr<FluidPointWindow> &fpw: mFluidPointWindows) {
        fpw->getFluidFields().mStiff.setZero();
        fpw->getFluidFields().mStiffR.setZero();
    }
}

// count info
int FluidSurfaceBoundary::
countInfo(const Messaging &msg) const {
    int count = 0;
    for (const auto &pw: mFluidPointWindows) {
        if (!msg.pointInSmallerRank(pw)) {
            count++;
        }
    }
    return count;
}
