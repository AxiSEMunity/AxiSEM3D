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
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
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

// release to domain
void GLLPoint::
release(const ABC &abc, const TimeScheme &timeScheme, Domain &domain) {
    //////////////////////////// reduce ////////////////////////////
    // variables to be reduced after setting up all Quads
    op1D_3D::tryReduceTo1D(mMassFluid);
    op1D_3D::tryReduceTo1D(mMassSolid);
    op1D_3D::tryReduceTo1D(mNormalSFU);
    op1D_3D::tryReduceTo1D(mNormalSFA);
    op1D_3D::tryReduceTo1D(mNormalTop);
    op1D_3D::tryReduceTo1D(mSumRhoDepth);
    op1D_3D::tryReduceTo1D(mGamma);
    
    
    //////////////////////////// point ////////////////////////////
    // solid point
    if (mMassSolid.rows() > 0) {
        // mass
        std::unique_ptr<const Mass> mass;
        if (mSumRhoDepth.rows() == 0) {
            if (mMassSolid.rows() == 1) {
                // 1D mass
                mass = std::make_unique<const Mass1D>(mMassSolid(0));
            } else {
                // 3D mass
                mass = std::make_unique<const Mass3D>(mMassSolid);
            }
        } else {
            if (mMassSolid.rows() == 1 && mSumRhoDepth.rows() == 1 &&
                mNormalTop.rows() == 1) {
                // 1D mass with ocean load
                eigen::DCol2 nsz;
                nsz << mNormalTop(0, 0), mNormalTop(0, 2);
                const eigen::DCol2 &nrt = geodesy::sz2rtheta(nsz, false);
                mass = std::make_unique<const MassOceanLoad1D>
                (mMassSolid(0), mSumRhoDepth(0) * nrt(0), nrt(1));
            } else {
                // 3D mass with ocean load
                // mass of ocean
                eigen::DColX massOcean;
                const eigen::DColX &area = mNormalTop.rowwise().norm();
                op1D_3D::times(mSumRhoDepth, area, massOcean);
                // unit normal
                const eigen::DMatX3 &unitNormal =
                mNormalTop.array().colwise() / area.array();
                // mass
                mass = std::make_unique<const MassOceanLoad3D>
                (op1D_3D::to3D(mMassSolid, mNr), op1D_3D::to3D(massOcean, mNr),
                 op1D_3D::to3D(unitNormal, mNr));
            }
        }
        // point
        mSolidPoint = std::make_shared<SolidPoint>
        (mNr, mCoords.transpose(), mGlobalTag, mass, timeScheme);
        // release
        domain.addSolidPoint(mSolidPoint);
    }
    
    // fluid point
    if (mMassFluid.rows() > 0) {
        // mass
        std::unique_ptr<const Mass> mass;
        if (mMassFluid.rows() == 1) {
            // 1D mass
            mass = std::make_unique<const Mass1D>(mMassFluid(0));
        } else {
            // 3D mass
            mass = std::make_unique<const Mass3D>(mMassFluid);
        }
        // point
        mFluidPoint = std::make_shared<FluidPoint>
        (mNr, mCoords.transpose(), mGlobalTag, mass, timeScheme);
        // release
        domain.addFluidPoint(mFluidPoint);
    }
    
    // check empty
    if (!mSolidPoint && !mFluidPoint) {
        throw std::runtime_error("GLLPoint::release || "
                                 "Point is neither solid nor fluid.");
    }
    
    
    //////////////////////////// boundaries ////////////////////////////
    // solid-fluid
    if (mSolidPoint && mFluidPoint) {
        std::unique_ptr<const SolidFluidCoupling> sfc = nullptr;
        // if a rank contains pure fluid, mNormalSFU is unitialized
        // in such case, mNormalSFU = mNormalSFA
        if (mNormalSFU.rows() == 0) {
            mNormalSFU = mNormalSFA;
        }
        // 1D or 3D
        if (mNormalSFA.rows() == 1 && mMassFluid.rows() == 1) {
            // 1D coupling
            sfc = std::make_unique<SolidFluidCoupling1D>
            (mSolidPoint, mFluidPoint,
             mNormalSFU(0, 0), mNormalSFU(0, 2),
             mNormalSFA(0, 0), mNormalSFA(0, 2),
             mMassFluid(0));
        } else {
            // 3D coupling
            sfc = std::make_unique<SolidFluidCoupling3D>
            (mSolidPoint, mFluidPoint,
             op1D_3D::to3D(mNormalSFU, mNr), op1D_3D::to3D(mNormalSFA, mNr),
             op1D_3D::to3D(mMassFluid, mNr));
        }
        domain.getSolidFluidBoundary()->addSolidFluidCoupling(sfc);
    }
    
    // Clayton ABC
    for (auto itm = mClaytonABC.begin(); itm != mClaytonABC.end(); itm++) {
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
                    (mFluidPoint, rhoVp(0), nABC.row(0).norm());
                } else {
                    // 3D fluid
                    clayton = std::make_unique<const ClaytonFluid3D>
                    (mFluidPoint, op1D_3D::to3D(rhoVp, mNr),
                     op1D_3D::to3D(nABC, mNr).rowwise().norm());
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
                    (mSolidPoint, rhoVp(0), rhoVs(0), nrt(0), nrt(1));
                } else {
                    // 3D solid
                    const eigen::DColX &area = nABC.rowwise().norm();
                    const eigen::DMatX3 &unitNormal =
                    nABC.array().colwise() / area.array();
                    clayton = std::make_unique<const ClaytonSolid3D>
                    (mSolidPoint,
                     op1D_3D::to3D(rhoVp, mNr), op1D_3D::to3D(rhoVs, mNr),
                     op1D_3D::to3D(area, mNr), op1D_3D::to3D(unitNormal, mNr));
                }
                domain.getAbsorbingBoundary()->addClaytonSolid(clayton);
            }
        }
    }
    
    // sponge ABC
    if (abc.sponge() && mCountGammasAdded > 0) {
        if (mFluidPoint) {
            std::unique_ptr<const Sponge<FluidPoint>> sponge_f = nullptr;
            if (mGamma.rows() == 1) {
                sponge_f = std::make_unique<const
                Sponge<FluidPoint>>(mFluidPoint, mGamma(0) / mCountGammasAdded);
            } else {
                sponge_f = std::make_unique<const
                Sponge<FluidPoint>>(mFluidPoint, mGamma / mCountGammasAdded);
            }
            domain.getAbsorbingBoundary()->addSpongeFluid(sponge_f);
        }
        if (mSolidPoint) {
            std::unique_ptr<const Sponge<SolidPoint>> sponge_s = nullptr;
            if (mGamma.rows() == 1) {
                sponge_s = std::make_unique<const
                Sponge<SolidPoint>>(mSolidPoint, mGamma(0) / mCountGammasAdded);
            } else {
                sponge_s = std::make_unique<const
                Sponge<SolidPoint>>(mSolidPoint, mGamma / mCountGammasAdded);
            }
            domain.getAbsorbingBoundary()->addSpongeSolid(sponge_s);
        }
    }
    
    // axial boundary
    if (mAxial) {
        if (mSolidPoint) {
            domain.getAxialBoundary()->addPoint(mSolidPoint);
        }
        if (mFluidPoint) {
            domain.getAxialBoundary()->addPoint(mFluidPoint);
        }
    }
    
    // fluid surface without ABC
    if (mFluidPoint && mSurface &&
        std::find(abc.getBoundaryKeys().begin(), abc.getBoundaryKeys().end(),
                  "TOP") == abc.getBoundaryKeys().end()) {
        domain.getFluidSurfaceBoundary()->addPoint(mFluidPoint);
    }
    
    // free dummy memory
    mMassFluid.resize(0);
    mMassSolid.resize(0);
    mNormalSFU.resize(0, 3);
    mNormalSFA.resize(0, 3);
    mClaytonABC.clear();
    mNormalTop.resize(0, 3);
    mSumRhoDepth.resize(0);
    mGamma.resize(0);
}
