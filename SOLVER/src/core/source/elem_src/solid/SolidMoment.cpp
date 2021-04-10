//
//  SolidMoment.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  moment source on solid element

#include "SolidMoment.hpp"
#include "Element.hpp"
#include <iostream>
// constructor
SolidMoment::SolidMoment(std::unique_ptr<STF> &stf,
                         const std::shared_ptr<const Element> &element,
                         int m, const eigen::CMatXN6 &pattern):
SolidSource(stf, element, m), mPattern(pattern) {
    // prepare
    mElement->prepareMomentSource(mM);
    
    // workspace
    if (sPattern.rows() < mPattern.rows()) {
        sPattern.resize(mPattern.rows(), spectral::nPEM * 6);
    }
}

// apply source at a time step
void SolidMoment::apply(double time) const {
    int nu_1 = (int)mPattern.rows();
    sPattern.topRows(nu_1) = mPattern * mSTF->getValue(time);
    mElement->addMomentSource(mM, sPattern, nu_1);
}
