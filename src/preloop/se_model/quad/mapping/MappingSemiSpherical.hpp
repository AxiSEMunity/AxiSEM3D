//
//  MappingSemiSpherical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/22/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

#ifndef MappingSemiSpherical_hpp
#define MappingSemiSpherical_hpp

#include "Mapping.hpp"

class MappingSemiSpherical: public Mapping {
public:
    // constructor
    MappingSemiSpherical(const eigen::DMat24 &nodalSZ): Mapping(nodalSZ) {
        // find curved outer
        const auto &r = nodalSZ.colwise().norm();
        double distTol = mMinEdgeLength * 1e-3;
        if (std::abs(r(0) - r(1)) < distTol && r(0) > r(2)) {
            mCurvedOuter = 0;
        } else if (std::abs(r(1) - r(2)) < distTol && r(1) > r(3)) {
            mCurvedOuter = 1;
        } else if (std::abs(r(2) - r(3)) < distTol && r(2) > r(0)) {
            mCurvedOuter = 2;
        } else if (std::abs(r(3) - r(0)) < distTol && r(3) > r(1)) {
            mCurvedOuter = 3;
        } else {
            throw
            std::runtime_error("MappingSemiSpherical::MappingSemiSpherical || "
                               "Error identifying curved outer edge.");
        }
    }
    
    // forward mapping: (ξ,η) -> (s,z)
    eigen::DCol2 mapping(const eigen::DCol2 &xieta) const {
        // compute coords rotated by Q2
        double s0, z0, s1, z1, r2, t2, t3, xi, eta;
        computeCoordsQ2(xieta, s0, z0, s1, z1, r2, t2, t3, xi, eta);
        
        // intermediate
        double s01 = ((1. - xi) * s0 + (1. + xi) * s1) * .5;
        double z01 = ((1. - xi) * z0 + (1. + xi) * z1) * .5;
        double t32 = ((1. - xi) * t3 + (1. + xi) * t2) * .5;
        double etam = (1. - eta) * .5;
        double r2p = (1. + eta) * r2 * .5;
        
        // compute in new system
        eigen::DCol2 sz;
        sz(0) = etam * s01 + r2p * sin(t32);
        sz(1) = etam * z01 + r2p * cos(t32);
        
        // rotate back
        return sOrthogQ2[mCurvedOuter].transpose() * sz;
    }
    
    // Jacobian: ∂(s,z) / ∂(ξ,η)
    eigen::DMat22 jacobian(const eigen::DCol2 &xieta) const {
        // compute coords rotated by Q2
        double s0, z0, s1, z1, r2, t2, t3, xi, eta;
        computeCoordsQ2(xieta, s0, z0, s1, z1, r2, t2, t3, xi, eta);
        
        // intermediate
        double s01 = ((1. - xi) * s0 + (1. + xi) * s1) * .5;
        double z01 = ((1. - xi) * z0 + (1. + xi) * z1) * .5;
        double t32 = ((1. - xi) * t3 + (1. + xi) * t2) * .5;
        double etam = (1. - eta) * .5;
        double r2p = (1. + eta) * r2 * .5;
        // derivative
        double s01_xi = (s1 - s0) * .5;
        double z01_xi = (z1 - z0) * .5;
        double t32_xi = (t2 - t3) * .5;
        double etam_eta = -.5;
        double r2p_eta = r2 * .5;
        
        // compute in new system
        eigen::DMat22 J;
        J(0, 0) = etam * s01_xi + r2p * cos(t32) * t32_xi;
        J(1, 0) = etam * z01_xi - r2p * sin(t32) * t32_xi;
        J(0, 1) = etam_eta * s01 + r2p_eta * sin(t32);
        J(1, 1) = etam_eta * z01 + r2p_eta * cos(t32);
        
        // rotate back
        const eigen::DMat22 &Q2 = sOrthogQ2[mCurvedOuter];
        return Q2.transpose() * J * Q2;
    }
    
private:
    // compute coords rotated by Q2
    void computeCoordsQ2(const eigen::DCol2 &xieta,
                         double &s0, double &z0, double &s1, double &z1,
                         double &r2, double &t2, double &t3,
                         double &xi, double &eta)
    const {
        // rotate CS such that mCurvedOuter lies on edge 2
        static eigen::DMat24 szQ2, rtQ2;
        static eigen::DCol2 xietaQ2;
        rotateQ2(xieta, szQ2, rtQ2, xietaQ2);
        
        // the coords are not permutated, so edge 2 should be
        // accessed by Mapping::cycle4(mCurvedOuter - 2)
        s0 = szQ2(0, Mapping::cycle4(mCurvedOuter - 2));
        z0 = szQ2(1, Mapping::cycle4(mCurvedOuter - 2));
        s1 = szQ2(0, Mapping::cycle4(mCurvedOuter - 1));
        z1 = szQ2(1, Mapping::cycle4(mCurvedOuter - 1));
        r2 = rtQ2(0, Mapping::cycle4(mCurvedOuter + 0));
        t2 = rtQ2(1, Mapping::cycle4(mCurvedOuter + 0));
        t3 = rtQ2(1, Mapping::cycle4(mCurvedOuter + 1));
        xi = xietaQ2(0);
        eta = xietaQ2(1);
        
        // handling angle jumps in atan2
        makeAnglesClose(t3, t2);
    }
};

#endif /* MappingSemiSpherical_hpp */
