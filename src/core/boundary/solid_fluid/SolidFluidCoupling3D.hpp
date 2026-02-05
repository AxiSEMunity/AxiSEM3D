//
//  SolidFluidCoupling3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition in 3D

#ifndef SolidFluidCoupling3D_hpp
#define SolidFluidCoupling3D_hpp

#include "SolidFluidCoupling.hpp"

class SolidFluidCoupling3D: public SolidFluidCoupling {
public:
    // constructor
    SolidFluidCoupling3D(const std::shared_ptr<SolidPoint> &sp,
                         const std::shared_ptr<FluidPoint> &fp,
                         const eigen::DMatX3 &n_unassmb,
                         const eigen::DMatX3 &n_assmb,
                         const eigen::DColX &massFluid);
    
private:
    // check compatibility
    void checkCompatibility(int nr) const;
    
public:
    // solid => fluid
    void coupleSolidToFluid(const eigen::CMatX3 &solidDispl,
                            eigen::CColX &fluidStiff) const;
    
    // fluid => solid
    void coupleFluidToSolid(const eigen::CColX &fluidStiff,
                            eigen::CMatX3 &solidStiff) const;
    
private:
    // These two normal vectors enable isochronous MPI communication for solid
    // and fluid domains. Though it is bad practice to mix MPI and physics,
    // but this trick can lead to significant performance boost.
    const eigen::RMatX3 mNormal_UnassembledMPI;
    const eigen::RMatX3 mNormal_AssembledMPI_InvMassFluid;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static eigen::RMatX3 sSolidR = eigen::RMatX3(0, 3);
    inline static eigen::CMatX3 sSolidC = eigen::CMatX3(0, 3);
    inline static eigen::RColX sFluidR = eigen::RColX(0);
    inline static eigen::CColX sFluidC = eigen::CColX(0);
};

#endif /* SolidFluidCoupling3D_hpp */
