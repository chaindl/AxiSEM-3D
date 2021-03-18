//
//  FluidPressure.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  pressure source on fluid element

#include "FluidPressure.hpp"
#include "Element.hpp"

// constructor
FluidPressure::FluidPressure(std::unique_ptr<STF> &stf,
                             const std::shared_ptr<const Element> &element,
                             int m, const eigen::CMatXN &pattern):
FluidSource(stf, element, m), mPattern(pattern) {
    // prepare
    mElement->preparePressureSource(mM);
    
    // workspace
    if (sPattern.rows() < mPattern.rows()) {
        sPattern.resize(mPattern.rows(), spectral::nPEM);
    }
}

// apply source at a time step
void FluidPressure::apply(double time) const {
    int nu_1 = (int)mPattern.rows();
    sPattern.topRows(nu_1) = mPattern * mSTF->getValue(time);
    mElement->addPressureSource(mM, sPattern, nu_1);
}
