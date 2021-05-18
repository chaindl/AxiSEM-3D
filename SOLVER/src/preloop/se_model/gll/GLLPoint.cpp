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
#include "PointWindow.hpp"
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
release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain,
    const std::shared_ptr<WindowInterpolator<numerical::Real>> &interpolator) {
    // point
    mPoint = std::make_shared<Point>
    (mGlobalTag, mCoords.transpose());
  
    // set up overlap groups 
    // (currently only needed to separate solid and fluid windows
    // in Point but this may be expanded for discontinuous fluids 
    // and/or crack modelling)
    eigen::DColX knotsWhole;
    std::vector<eigen::DColX> relPhis(mWindows.size());
    std::vector<eigen::IColX> posIndices(mWindows.size());
    int nrWS;
    computeWindowSumSampling(knotsWhole, relPhis, posIndices, nrWS);
    
    bool hasSolidWins = std::any_of(mMassSolid.begin(), mMassSolid.end(), 
                        [](eigen::DColX mass) {return mass.rows() > 0;});
    bool hasFluidWins = std::any_of(mMassFluid.begin(), mMassFluid.end(), 
                        [](eigen::DColX mass) {return mass.rows() > 0;});
    
    int interpTag = -1;
    if (knotsWhole.size() > 0) {
        int dim = hasSolidWins ? 3 : 1;
        int nr = (*std::max_element(mWindows.begin(), mWindows.end(), 
            [](const eigen::DMatX2 &win1, const eigen::DMatX2 &win2) {return (win1.rows() < win2.rows());})).rows();
        interpTag = interpolator->addSplineFitting(knotsWhole.cast<numerical::Real>(), dim);
    }
    
    std::shared_ptr<SFWindowSum<3>> sog = nullptr;
    std::shared_ptr<SFWindowSum<1>> fog = nullptr;
    if (hasSolidWins) {
        sog = std::make_shared<SFWindowSum<3>>(interpolator, interpTag, nrWS);
    }
    if (hasFluidWins) {
        fog = std::make_shared<SFWindowSum<1>>(interpolator, interpTag, nrWS);
    }
    
    for (int m = 0; m < mWindows.size(); m++) {
        std::shared_ptr<PointWindow> fpw = nullptr;
        std::shared_ptr<PointWindow> spw = nullptr;
        
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
            // interpolator
            interpTag = -1;
            if (relPhis[m].size() > 0) interpTag = interpolator->addSplineFitting(
                                                   eigen::RColX::LinSpaced(mWindows[m].rows(), 0, 1), 3);
            
            // release
            if (mGlobalWin[m]) {
                spw = std::make_shared<SolidRCPointWindow<3, numerical::ComplexR>>(mass, mWindows[m].cast<numerical::Real>(), mPoint->getMeshTag());
                spw->checkCompatibility(timeScheme);
            } else {
                spw = std::make_shared<SolidRCPointWindow<3, numerical::Real>>(mass, mWindows[m].cast<numerical::Real>(), mPoint->getMeshTag());
                spw->checkCompatibility(timeScheme);
            }
            mPoint->addWindow(spw);
            sog->addWindow(spw, interpTag, posIndices[m], relPhis[m].cast<numerical::Real>());
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
            // interpolator
            interpTag = -1;
            if (relPhis[m].size() > 0) interpTag = interpolator->addSplineFitting(
                                                   eigen::RColX::LinSpaced(mWindows[m].rows(), 0, 1), 1);
            
            // release
            if (mGlobalWin[m]) {
                fpw = std::make_shared<FluidRCPointWindow<1, numerical::ComplexR>>(mass, mWindows[m].cast<numerical::Real>(), mPoint->getMeshTag());
                fpw->checkCompatibility(timeScheme);
            } else {
                fpw = std::make_shared<FluidRCPointWindow<1, numerical::Real>>(mass, mWindows[m].cast<numerical::Real>(), mPoint->getMeshTag());
                fpw->checkCompatibility(timeScheme);
            }
            mPoint->addWindow(fpw);
            fog->addWindow(fpw, interpTag, posIndices[m], relPhis[m].cast<numerical::Real>());
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
                        if (mGlobalWin[m]) {
                            clayton = std::make_unique<const ClaytonFluid1D<numerical::ComplexR>>
                                      (fpw, rhoVp(0), nABC.row(0).norm());
                        } else {
                            clayton = std::make_unique<const ClaytonFluid1D<numerical::Real>>
                                     (fpw, rhoVp(0), nABC.row(0).norm());  
                        }
                    } else {
                        // 3D fluid
                        if (mGlobalWin[m]) {
                            clayton = std::make_unique<const ClaytonFluid3D_C>
                                      (fpw, op1D_3D::to3D(rhoVp, mWindows[m].rows()),
                                      op1D_3D::to3D(nABC, mWindows[m].rows()).rowwise().norm());
                        } else {
                            clayton = std::make_unique<const ClaytonFluid3D_R>
                                      (fpw, op1D_3D::to3D(rhoVp, mWindows[m].rows()),
                                      op1D_3D::to3D(nABC, mWindows[m].rows()).rowwise().norm());
                        }
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
                        if (mGlobalWin[m]) {
                            clayton = std::make_unique<const ClaytonSolid1D<numerical::ComplexR>>
                                      (spw, rhoVp(0), rhoVs(0), nrt(0), nrt(1));
                        } else {
                            clayton = std::make_unique<const ClaytonSolid1D<numerical::Real>>
                                      (spw, rhoVp(0), rhoVs(0), nrt(0), nrt(1));
                        }
                    } else {
                        // 3D solid
                        const eigen::DColX &area = nABC.rowwise().norm();
                        const eigen::DMatX3 &unitNormal =
                        nABC.array().colwise() / area.array();
                        if (mGlobalWin[m]) {
                            clayton = std::make_unique<const ClaytonSolid3D_C>
                                      (spw, op1D_3D::to3D(rhoVp, mWindows[m].rows()), op1D_3D::to3D(rhoVs, mWindows[m].rows()),
                                      op1D_3D::to3D(area, mWindows[m].rows()), op1D_3D::to3D(unitNormal, mWindows[m].rows()));
                        } else {
                            clayton = std::make_unique<const ClaytonSolid3D_R>
                                      (spw, op1D_3D::to3D(rhoVp, mWindows[m].rows()), op1D_3D::to3D(rhoVs, mWindows[m].rows()),
                                      op1D_3D::to3D(area, mWindows[m].rows()), op1D_3D::to3D(unitNormal, mWindows[m].rows()));
                        }
                    }
                    domain.getAbsorbingBoundary()->addClaytonSolid(clayton);
                }
            }
        }
        
        // sponge ABC
        if (abc.sponge() && mCountGammasAdded[m] > 0) {
            if (fpw) {
                std::unique_ptr<const Sponge<1>> sponge_f = nullptr;
                if (mGamma[m].rows() == 1) {
                    if (mGlobalWin[m]) {
                        sponge_f = std::make_unique<const
                        Sponge_C<1>>(fpw, mGamma[m](0) / mCountGammasAdded[m]);
                    } else {
                        sponge_f = std::make_unique<const
                        Sponge_R<1>>(fpw, mGamma[m](0) / mCountGammasAdded[m]);
                    }
                } else {
                    if (mGlobalWin[m]) {
                        sponge_f = std::make_unique<const
                        Sponge_C<1>>(fpw, mGamma[m] / mCountGammasAdded[m]);
                    } else {
                        sponge_f = std::make_unique<const
                        Sponge_R<1>>(fpw, mGamma[m] / mCountGammasAdded[m]);
                    }
                }
                domain.getAbsorbingBoundary()->addSpongeFluid(sponge_f);
            }
            if (spw) {
                std::unique_ptr<const Sponge<3>> sponge_s = nullptr;
                if (mGamma[m].rows() == 1) {
                    if (mGlobalWin[m]) {
                        sponge_s = std::make_unique<const
                        Sponge_C<3>>(spw, mGamma[m](0) / mCountGammasAdded[m]);
                    } else {
                        sponge_s = std::make_unique<const
                        Sponge_R<3>>(spw, mGamma[m](0) / mCountGammasAdded[m]);
                    }
                } else {
                    if (mGlobalWin[m]) {
                        sponge_s = std::make_unique<const
                        Sponge_C<3>>(spw, mGamma[m] / mCountGammasAdded[m]);
                    } else {
                        sponge_s = std::make_unique<const
                        Sponge_R<3>>(spw, mGamma[m] / mCountGammasAdded[m]);
                    }
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
