//
//  SolidFluidCoupling1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition in 1D

#ifndef SolidFluidCoupling1D_hpp
#define SolidFluidCoupling1D_hpp

#include "SolidFluidCoupling.hpp"

class SolidFluidCoupling1D: public SolidFluidCoupling {
public:
    // constructor
    SolidFluidCoupling1D(const std::shared_ptr<SolidPoint> &sp,
                         const std::shared_ptr<FluidPoint> &fp,
                         double ns_unassmb, double nz_unassmb,
                         double ns_assmb, double nz_assmb,
                         double massFluid);
    
    // solid => fluid
    void coupleSolidToFluid(const eigen::CMatX3 &solidDispl,
                            eigen::CColX &fluidStiff) const {
        fluidStiff += mNormalS_UnassembledMPI * solidDispl.col(0);
        fluidStiff += mNormalZ_UnassembledMPI * solidDispl.col(2);
    }
    
    // fluid => solid
    void coupleFluidToSolid(const eigen::CColX &fluidStiff,
                            eigen::CMatX3 &solidStiff) const {
        solidStiff.col(0) -= mNormalS_AssembledMPI_InvMassFluid * fluidStiff;
        solidStiff.col(2) -= mNormalZ_AssembledMPI_InvMassFluid * fluidStiff;
    }
    
private:
    // These two normal vectors enable isochronous MPI communication for solid
    // and fluid domains. Though it is bad practice to mix MPI and physics,
    // but this trick can lead to significant performance boost.
    const numerical::Real mNormalS_UnassembledMPI;
    const numerical::Real mNormalZ_UnassembledMPI;
    const numerical::Real mNormalS_AssembledMPI_InvMassFluid;
    const numerical::Real mNormalZ_AssembledMPI_InvMassFluid;
};

#endif /* SolidFluidCoupling1D_hpp */
