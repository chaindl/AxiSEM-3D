//
//  ClaytonFluid1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 1D

#include "ClaytonFluid1D.hpp"
#include "PointWindow.hpp"

// apply ABC
void ClaytonFluid1D::apply() const {
    // get fields
    const eigen::CColX &veloc = mPointWindow->getFluidFields().mVeloc;
    eigen::CColX &stiff = mPointWindow->getFluidFields().mStiff;
    
    // apply
    stiff -= veloc * mAreaOverRhoVp;
}
