//
//  MappingSpherical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/21/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  spherical mapping

#ifndef MappingSpherical_hpp
#define MappingSpherical_hpp

#include "Mapping.hpp"

class MappingSpherical: public Mapping {
public:
    // constructor
    MappingSpherical(const eigen::DMat24 &nodalSZ): Mapping(nodalSZ) {
        // find curved outer
        const auto &r = nodalSZ.colwise().norm();
        double distTol = mMinEdgeLength * 1e-3;
        if (std::abs(r(0) - r(1)) < distTol &&
            std::abs(r(2) - r(3)) < distTol) {
            mCurvedOuter = r(2) > r(0) ? 2 : 0;
        } else if (std::abs(r(1) - r(2)) < distTol &&
                   std::abs(r(3) - r(0)) < distTol) {
            mCurvedOuter = r(3) > r(1) ? 3 : 1;
        } else {
            throw std::runtime_error("MappingSpherical::MappingSpherical || "
                                     "Error identifying curved outer edge.");
        }
    }
    
    // forward mapping: (ξ,η) -> (s,z)
    eigen::DCol2 mapping(const eigen::DCol2 &xieta) const {
        // compute coords rotated by Q2
        double r0, r2, t0, t1, t2, t3, xi, eta;
        computeCoordsQ2(xieta, r0, r2, t0, t1, t2, t3, xi, eta);
        
        // intermediate
        double t01 = ((1. - xi) * t0 + (1. + xi) * t1) * .5;
        double t32 = ((1. - xi) * t3 + (1. + xi) * t2) * .5;
        double r0m = (1. - eta) * r0 * .5;
        double r2p = (1. + eta) * r2 * .5;
        
        // compute in new system
        eigen::DCol2 sz;
        sz(0) = r0m * sin(t01) + r2p * sin(t32);
        sz(1) = r0m * cos(t01) + r2p * cos(t32);
        
        // rotate back
        return sOrthogQ2[mCurvedOuter].transpose() * sz;
    }
    
    // Jacobian: ∂(s,z) / ∂(ξ,η)
    eigen::DMat22 jacobian(const eigen::DCol2 &xieta) const {
        // compute coords rotated by Q2
        double r0, r2, t0, t1, t2, t3, xi, eta;
        computeCoordsQ2(xieta, r0, r2, t0, t1, t2, t3, xi, eta);
        
        // intermediate
        double t01 = ((1. - xi) * t0 + (1. + xi) * t1) * .5;
        double t32 = ((1. - xi) * t3 + (1. + xi) * t2) * .5;
        double r0m = (1. - eta) * r0 * .5;
        double r2p = (1. + eta) * r2 * .5;
        // derivative
        double t01_xi = (t1 - t0) * .5;
        double t32_xi = (t2 - t3) * .5;
        double r0m_eta = -r0 * .5;
        double r2p_eta = r2 * .5;
        
        // compute in new system
        eigen::DMat22 J;
        J(0, 0) = (r0m * cos(t01) * t01_xi + r2p * cos(t32) * t32_xi);
        J(1, 0) = (r0m * sin(t01) * t01_xi + r2p * sin(t32) * t32_xi) * (-1.);
        J(0, 1) = r0m_eta * sin(t01) + r2p_eta * sin(t32);
        J(1, 1) = r0m_eta * cos(t01) + r2p_eta * cos(t32);
        
        // rotate back
        const eigen::DMat22 &Q2 = sOrthogQ2[mCurvedOuter];
        return Q2.transpose() * J * Q2;
    }
    
private:
    // compute coords rotated by Q2
    void computeCoordsQ2(const eigen::DCol2 &xieta,
                         double &r0, double &r2, double &t0, double &t1,
                         double &t2, double &t3, double &xi, double &eta)
    const {
        // rotate CS such that mCurvedOuter lies on edge 2
        static eigen::DMat24 szQ2, rtQ2;
        static eigen::DCol2 xietaQ2;
        rotateQ2(xieta, szQ2, rtQ2, xietaQ2);
        
        // the coords are not permutated, so edge 2 should be
        // accessed by Mapping::cycle4(mCurvedOuter - 2)
        r0 = rtQ2(0, Mapping::cycle4(mCurvedOuter - 2));
        r2 = rtQ2(0, Mapping::cycle4(mCurvedOuter + 0));
        t0 = rtQ2(1, Mapping::cycle4(mCurvedOuter - 2));
        t1 = rtQ2(1, Mapping::cycle4(mCurvedOuter - 1));
        t2 = rtQ2(1, Mapping::cycle4(mCurvedOuter + 0));
        t3 = rtQ2(1, Mapping::cycle4(mCurvedOuter + 1));
        xi = xietaQ2(0);
        eta = xietaQ2(1);
        
        // handling angle jumps in atan2
        makeAnglesClose(t3, t2);
        makeAnglesClose(t0, t1);
    }
};

#endif /* MappingSpherical_hpp */
