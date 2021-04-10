//
//  SolidFluidCoupling.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#include "SolidFluidCoupling.hpp"

// point
#include "SolidPointWindow.hpp"
#include "FluidPointWindow.hpp"

class Point;
#include <iostream>
// constructor
SolidFluidCoupling::
SolidFluidCoupling(const std::shared_ptr<SolidPointWindow> &spw,
                   const std::shared_ptr<FluidPointWindow> &fpw):
mSolidPointWindow(spw), mFluidPointWindow(fpw) {
    if (mSolidPointWindow->getMeshTag() != mFluidPointWindow->getMeshTag()) {
        throw std::runtime_error("SolidFluidCoupling::SolidFluidCoupling || "
                                 "The coupled solid and fluid points have "
                                 "different mesh tags (positions).");
    }
}

// compute coupling
void SolidFluidCoupling::apply() const {
    // this order matters!
    coupleSolidToFluid(mSolidPointWindow->getSolidFields().mDispl,
                       mFluidPointWindow->getFluidFields().mStiff);
    coupleFluidToSolid(mFluidPointWindow->getFluidFields().mStiff,
                       mSolidPointWindow->getSolidFields().mStiff);
}


////////////////////////////// virtual //////////////////////////////
// check compatibility
void SolidFluidCoupling::checkCompatibility(int nr) const {
    if (mSolidPointWindow->getNr() != nr || mFluidPointWindow->getNr() != nr) {
        throw std::runtime_error("SolidFluidCoupling::checkCompatibility || "
                                 "Incompatible sizes.");
    }
}
