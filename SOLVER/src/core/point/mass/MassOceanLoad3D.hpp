//
//  MassOceanLoad3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  3D mass with ocean load
//  used on solid-fluid boundary when the fluid is modelled as load
//  a = F.n n / (m + m0) + (F - F.n n) / m
//    = [m F.n n + (m + m0) (F - F.n n)] / [m (m + m0)]
//    = [(m + m0) F - m0 F.n n] / [m (m + m0)]
//    = F / m - m0 / [m (m + m0)] F.n n
//    = im F - F.k k, with k = sqrt(m0 / [m (m + m0)]) n
//  a -> acceleration
//  F -> force
//  n -> unit normal of surface
//  m -> mass of solid
//  m0 -> mass of water column above

#ifndef MassOceanLoad3D_hpp
#define MassOceanLoad3D_hpp

#include "Mass.hpp"

class MassOceanLoad3D: public Mass {
public:
    // constructor
    MassOceanLoad3D(const eigen::DColX &mass,
                    const eigen::DColX &massOcean,
                    const eigen::DMatX3 &unitNormal);
    
    // check compatibility
    void checkCompatibility(int nr, bool solid) const;
    
    // compute accel in-place for fluid
    void computeAccel(eigen::CColX &stiff1) const {
        throw std::runtime_error("MassOceanLoad3D::computeAccel || "
                                 "Incompatible types: "
                                 "ocean load on fluid point.");
    }
    
    // compute accel in-place for solid
    void computeAccel(eigen::CMatX3 &stiff3) const;
    
private:
    // im = 1 / m
    const eigen::RColX mIM;
    // k = sqrt(m0 / [m (m + m0)]) n
    const eigen::RMatX3 mK;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    // F = FFT(stiff)
    inline static eigen::RMatX3 sF = eigen::RMatX3(0, 3);
    // a = im F - F.k k
    inline static eigen::RMatX3 sA = eigen::RMatX3(0, 3);
};

#endif /* MassOceanLoad3D_hpp */
