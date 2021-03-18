//
//  GLLPoint.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  GLL point for preloop processing
//  generator of Point and boundary conditions in core

#include "GLLPoint.hpp"
#include "Domain.hpp"

// point
#include "Point.hpp"
#include "SolidPointWindow.hpp"
#include "FluidPointWindow.hpp"
#include "Mass1D.hpp"
#include "Mass3D.hpp"
#include "MassOceanLoad1D.hpp"
#include "MassOceanLoad3D.hpp"
#include "geodesy.hpp"

// boundary
#include "SolidFluidBoundary.hpp"
#include "SolidFluidCoupling1D.hpp"
#include "SolidFluidCoupling3D.hpp"
#include "ABC.hpp"
#include "AbsorbingBoundary.hpp"
#include "ClaytonSolid1D.hpp"
#include "ClaytonSolid3D.hpp"
#include "ClaytonFluid1D.hpp"
#include "ClaytonFluid3D.hpp"
#include "Sponge.hpp"
#include "AxialBoundary.hpp"
#include "FluidSurfaceBoundary.hpp"

#include "WindowSum.hpp"

// release to domain
void GLLPoint::
release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain) {
    // point
    mPoint = std::make_shared<Point>
    (mGlobalTag, mCoords.transpose());
  
    // set up overlap groups 
    // (currently only needed to separate solid and fluid windows
    // in Point but this may be expanded for discontinuous fluids 
    // and/or crack modelling)
    std::vector<int> aligned;
    eigen::DColX windowSumPhi;
    computeWindowSumSampling(windowSumPhi, aligned);
    
    bool hasSolidWins = std::any_of(mMassSolid.begin(), mMassSolid.end(), [](eigen::DColX mass) {return mass.rows() > 0;});
    bool hasFluidWins = std::any_of(mMassFluid.begin(), mMassFluid.end(), [](eigen::DColX mass) {return mass.rows() > 0;});
    std::shared_ptr<SolidWindowSum> sog = nullptr;
    std::shared_ptr<FluidWindowSum> fog = nullptr;
    if (hasSolidWins) {
        sog = std::make_shared<SolidWindowSum>(windowSumPhi.cast<numerical::Real>(), (mWindows.size() == 1));
    }
    if (hasFluidWins) {
        fog = std::make_shared<FluidWindowSum>(windowSumPhi.cast<numerical::Real>(), (mWindows.size() == 1));
    }

    for (int m = 0; m < mWindows.size(); m++) {
        std::shared_ptr<FluidPointWindow> fpw = nullptr;
        std::shared_ptr<SolidPointWindow> spw = nullptr;
        
        //////////////////////////// reduce ////////////////////////////
        // variables to be reduced after setting up all Quads
        op1D_3D::tryReduceTo1D(mMassFluid[m]);
        op1D_3D::tryReduceTo1D(mMassSolid[m]);
        op1D_3D::tryReduceTo1D(mNormalSFU[m]);
        op1D_3D::tryReduceTo1D(mNormalSFA[m]);
        op1D_3D::tryReduceTo1D(mNormalTop[m]);
        op1D_3D::tryReduceTo1D(mSumRhoDepth[m]);
        op1D_3D::tryReduceTo1D(mGamma[m]);
        
        //////////////////////////// point ////////////////////////////
        // solid point
        if (mMassSolid[m].rows() > 0) {
            // mass
            std::unique_ptr<const Mass> mass;
            if (mSumRhoDepth[m].rows() == 0) {
                if (mMassSolid[m].rows() == 1) {
                    // 1D mass
                    mass = std::make_unique<const Mass1D>(mMassSolid[m](0));
                } else {
                    // 3D mass
                    mass = std::make_unique<const Mass3D>(mMassSolid[m]);
                }
            } else {
                if (mMassSolid[m].rows() == 1 && mSumRhoDepth[m].rows() == 1 &&
                    mNormalTop[m].rows() == 1) {
                    // 1D mass with ocean load
                    eigen::DCol2 nsz;
                    nsz << mNormalTop[m](0, 0), mNormalTop[m](0, 2);
                    const eigen::DCol2 &nrt = geodesy::sz2rtheta(nsz, false);
                    mass = std::make_unique<const MassOceanLoad1D>
                    (mMassSolid[m](0), mSumRhoDepth[m](0) * nrt(0), nrt(1));
                } else {
                    // 3D mass with ocean load
                    // mass of ocean
                    eigen::DColX massOcean;
                    const eigen::DColX &area = mNormalTop[m].rowwise().norm();
                    op1D_3D::times(mSumRhoDepth[m], area, massOcean);
                    // unit normal
                    const eigen::DMatX3 &unitNormal =
                    mNormalTop[m].array().colwise() / area.array();
                    // mass
                    mass = std::make_unique<const MassOceanLoad3D>
                    (op1D_3D::to3D(mMassSolid[m], mWindows[m].rows()), op1D_3D::to3D(massOcean, mWindows[m].rows()),
                     op1D_3D::to3D(unitNormal, mWindows[m].rows()));
                }
            }
            // release
            spw = std::make_shared<SolidPointWindow>(mWindows[m].cast<numerical::Real>(), mass, timeScheme, mPoint);
            mPoint->addWindow(spw);
            sog->addSolidWindow(spw, aligned[m]);
        }
        
        // fluid point
        if (mMassFluid[m].rows() > 0) {
            // mass
            std::unique_ptr<const Mass> mass;
            if (mMassFluid[m].rows() == 1) {
                // 1D mass
                mass = std::make_unique<const Mass1D>(mMassFluid[m](0));
            } else {
                // 3D mass
                mass = std::make_unique<const Mass3D>(mMassFluid[m]);
            }
            // release
            fpw = std::make_shared<FluidPointWindow>(mWindows[m].cast<numerical::Real>(), mass, timeScheme, mPoint);
            mPoint->addWindow(fpw);
            fog->addFluidWindow(fpw, aligned[m]);
        }
        
        // check empty
        if (!spw && !fpw) {
            throw std::runtime_error("GLLPoint::release || "
                                     "Window is neither solid nor fluid.");
        }
        
        
        //////////////////////////// boundaries ////////////////////////////
        // solid-fluid
        if (spw && fpw) {
            std::unique_ptr<const SolidFluidCoupling> sfc = nullptr;
            // if a rank contains pure fluid, mNormalSFU is unitialized
            // in such case, mNormalSFU = mNormalSFA
            if (mNormalSFU[m].rows() == 0) {
                mNormalSFU[m] = mNormalSFA[m];
            }
            // 1D or 3D
            if (mNormalSFA[m].rows() == 1 && mMassFluid[m].rows() == 1) {
                // 1D coupling
                sfc = std::make_unique<SolidFluidCoupling1D>
                (spw, fpw,
                 mNormalSFU[m](0, 0), mNormalSFU[m](0, 2),
                 mNormalSFA[m](0, 0), mNormalSFA[m](0, 2),
                 mMassFluid[m](0));
            } else {
                // 3D coupling
                sfc = std::make_unique<SolidFluidCoupling3D>
                (spw, fpw,
                 op1D_3D::to3D(mNormalSFU[m], mWindows[m].rows()), 
                 op1D_3D::to3D(mNormalSFA[m], mWindows[m].rows()),
                 op1D_3D::to3D(mMassFluid[m], mWindows[m].rows()));
            }
            domain.getSolidFluidBoundary()->addSolidFluidCoupling(sfc);
        }
        
        // Clayton ABC
        for (auto itm = mClaytonABC[m].begin(); itm != mClaytonABC[m].end(); itm++) {
            for (auto itv = itm->second.begin(); itv != itm->second.end(); itv++) {
                bool fluid = std::get<0>(*itv);
                const eigen::DMatX3 &nABC = std::get<1>(*itv);
                const eigen::DColX &rhoVp = std::get<2>(*itv);
                const eigen::DColX &rhoVs = std::get<3>(*itv);
                // fluid
                if (fluid) {
                    std::unique_ptr<const ClaytonFluid> clayton = nullptr;
                    if (nABC.rows() == 1 && rhoVp.rows() == 1) {
                        // 1D fluid
                        clayton = std::make_unique<const ClaytonFluid1D>
                        (fpw, rhoVp(0), nABC.row(0).norm());
                    } else {
                        // 3D fluid
                        clayton = std::make_unique<const ClaytonFluid3D>
                        (fpw, op1D_3D::to3D(rhoVp, mWindows[m].rows()),
                         op1D_3D::to3D(nABC, mWindows[m].rows()).rowwise().norm());
                    }
                    domain.getAbsorbingBoundary()->addClaytonFluid(clayton);
                } else {
                    std::unique_ptr<const ClaytonSolid> clayton = nullptr;
                    if (nABC.rows() == 1 && rhoVp.rows() == 1 &&
                        rhoVs.rows() == 1) {
                        // 1D solid
                        eigen::DCol2 nsz;
                        nsz << nABC(0, 0), nABC(0, 2);
                        const eigen::DCol2 &nrt = geodesy::sz2rtheta(nsz, false);
                        clayton = std::make_unique<const ClaytonSolid1D>
                        (spw, rhoVp(0), rhoVs(0), nrt(0), nrt(1));
                    } else {
                        // 3D solid
                        const eigen::DColX &area = nABC.rowwise().norm();
                        const eigen::DMatX3 &unitNormal =
                        nABC.array().colwise() / area.array();
                        clayton = std::make_unique<const ClaytonSolid3D>
                        (spw,
                         op1D_3D::to3D(rhoVp, mWindows[m].rows()), op1D_3D::to3D(rhoVs, mWindows[m].rows()),
                         op1D_3D::to3D(area, mWindows[m].rows()), op1D_3D::to3D(unitNormal, mWindows[m].rows()));
                    }
                    domain.getAbsorbingBoundary()->addClaytonSolid(clayton);
                }
            }
        }
        
        // sponge ABC
        if (abc.sponge() && mCountGammasAdded[m] > 0) {
            if (fpw) {
                std::unique_ptr<const Sponge<FluidPointWindow>> sponge_f = nullptr;
                if (mGamma[m].rows() == 1) {
                    sponge_f = std::make_unique<const
                    Sponge<FluidPointWindow>>(fpw, mGamma[m](0) / mCountGammasAdded[m]);
                } else {
                    sponge_f = std::make_unique<const
                    Sponge<FluidPointWindow>>(fpw, mGamma[m] / mCountGammasAdded[m]);
                }
                domain.getAbsorbingBoundary()->addSpongeFluid(sponge_f);
            }
            if (spw) {
                std::unique_ptr<const Sponge<SolidPointWindow>> sponge_s = nullptr;
                if (mGamma[m].rows() == 1) {
                    sponge_s = std::make_unique<const
                    Sponge<SolidPointWindow>>(spw, mGamma[m](0) / mCountGammasAdded[m]);
                } else {
                    sponge_s = std::make_unique<const
                    Sponge<SolidPointWindow>>(spw, mGamma[m] / mCountGammasAdded[m]);
                }
                domain.getAbsorbingBoundary()->addSpongeSolid(sponge_s);
            }
        }
        
        // fluid surface without ABC
        if (fpw && mSurface &&
            std::find(abc.getBoundaryKeys().begin(), abc.getBoundaryKeys().end(),
                      "TOP") == abc.getBoundaryKeys().end()) {
            domain.getFluidSurfaceBoundary()->addPointWindow(fpw);
        }
        
        // axial boundary
        if (mAxial) {
            if (spw) {
                domain.getAxialBoundary()->addPointWindow(spw);
            }
            if (fpw) {
                domain.getAxialBoundary()->addPointWindow(fpw);
            }
        }
        
        mSolidPointWindows.push_back(spw);
        mFluidPointWindows.push_back(fpw);
    }
    
    if (sog) mPoint->addWindowSum(sog);
    if (fog) mPoint->addWindowSum(fog);
    
    // free dummy memory
    mMassFluid.clear();
    mMassSolid.clear();
    mNormalSFU.clear();
    mNormalSFA.clear();
    mClaytonABC.clear();
    mNormalTop.clear();
    mSumRhoDepth.clear();
    mGamma.clear();
    
    // release
    domain.addPoint(mPoint);
}
