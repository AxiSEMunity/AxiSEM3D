//
//  ClaytonSolid3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 3D
//  f = v.n n * rpa + (v - v.n n) * rsa
//    = rsa v + v.k k, with k = sqrt(rpa - rsa) n
//  f -> stifness force
//  v -> velocity
//  n -> unit normal of surface
//  rsa -> rho * vs * area
//  rpa -> rho * vp * area

#ifndef ClaytonSolid3D_hpp
#define ClaytonSolid3D_hpp

#include "ClaytonSolid.hpp"
#include "eigen_point.hpp"

class ClaytonSolid3D: public ClaytonSolid {
public:
    // constructor
    ClaytonSolid3D(const std::shared_ptr<SolidPoint> &sp,
                   const eigen::DColX &rhoVp, const eigen::DColX &rhoVs,
                   const eigen::DColX &area,
                   const eigen::DMatX3 &unitNormal):
    ClaytonSolid(sp),
    mRSA(rhoVs.cwiseProduct(area).cast<numerical::Real>()),
    mK(((rhoVp - rhoVs).cwiseProduct(area).cwiseSqrt().asDiagonal()
        * unitNormal).cast<numerical::Real>()) {
        // check compatibility
        checkCompatibility();
    }
    
private:
    // check compatibility
    void checkCompatibility();
    
public:
    // apply ABC
    void apply() const;
    
private:
    // rsa = rho * vs * area
    const eigen::RColX mRSA;
    // k = sqrt(rpa - rsa) n
    const eigen::RMatX3 mK;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // workspace
    // V = FFT(velocity)
    inline static eigen::RMatX3 sVR = eigen::RMatX3(0, 3);
    // a = rsa V + V.k k
    inline static eigen::RMatX3 sAR = eigen::RMatX3(0, 3);
    inline static eigen::CMatX3 sAC = eigen::CMatX3(0, 3);
};

#endif /* ClaytonSolid3D_hpp */
