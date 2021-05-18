//
//  SolidFluidCoupling1D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition in 1D

#include "SolidFluidCoupling1D.hpp"
#include "PointWindow.hpp"

// constructor
SolidFluidCoupling1D::
SolidFluidCoupling1D(const std::shared_ptr<PointWindow> &spw,
                     const std::shared_ptr<PointWindow> &fpw,
                     double ns_unassmb, double nz_unassmb,
                     double ns_assmb, double nz_assmb,
                     double massFluid):
SolidFluidCoupling(spw, fpw),
mNormalS_UnassembledMPI(ns_unassmb),
mNormalZ_UnassembledMPI(nz_unassmb),
mNormalS_AssembledMPI_InvMassFluid(ns_assmb / massFluid),
mNormalZ_AssembledMPI_InvMassFluid(nz_assmb / massFluid) {
    checkCompatibility(mSolidPointWindow->getNr());
}
