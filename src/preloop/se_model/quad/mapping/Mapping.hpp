//
//  Mapping.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/21/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  geometric mapping on quadrilateral

#ifndef Mapping_hpp
#define Mapping_hpp

#include "geodesy.hpp"
#include "eigen_sem.hpp"

class Mapping {
public:
    // constructor
    Mapping(const eigen::DMat24 &nodalSZ): mNodalSZ(nodalSZ) {
        // minimum edge length
        eigen::DRow4 edgeLength;
        edgeLength(0) = (nodalSZ.col(0) - nodalSZ.col(1)).norm();
        edgeLength(1) = (nodalSZ.col(1) - nodalSZ.col(2)).norm();
        edgeLength(2) = (nodalSZ.col(2) - nodalSZ.col(3)).norm();
        edgeLength(3) = (nodalSZ.col(3) - nodalSZ.col(0)).norm();
        mMinEdgeLength = edgeLength.minCoeff();
    }
    
    // destructor
    virtual ~Mapping() = default;
    
    // forward mapping: (ξ,η) -> (s,z)
    virtual eigen::DCol2 mapping(const eigen::DCol2 &xieta) const = 0;
    
    // Jacobian: ∂(s,z) / ∂(ξ,η)
    virtual eigen::DMat22 jacobian(const eigen::DCol2 &xieta) const = 0;
    
    // inverse mapping: (s,z) -> (ξ,η)
    // return true if (s,z) is inside this element
    bool inverseMapping(const eigen::DCol2 &sz, eigen::DCol2 &xieta,
                        double maxIter = 10, double tolerance = 1e-9) const {
        // Newton
        xieta.setZero();
        double absSZ = tolerance * mMinEdgeLength;
        int iter = 0;
        for (; iter < maxIter; iter++) {
            const eigen::DCol2 &dsz = sz - mapping(xieta);
            if (dsz.norm() < absSZ) {
                break;
            }
            xieta += jacobian(xieta).inverse() * dsz;
        }
        
        // Newton failed
        if (iter == maxIter) {
            return false;
        }
        
        // check inside/outside with a larger tolerance
        double boundXiEta = 1. + tolerance * 20.;
        return xieta.minCoeff() > -boundXiEta && xieta.maxCoeff() < boundXiEta;
    }
    
    // normal
    eigen::DCol2 normal(int edge, const eigen::DMat22 &J) const {
        eigen::DCol2 normal;
        if (edge == 0) {
            // xi increases from -1 to 1, eta = -1
            // n = (0, 1, 0) x (ds/dxi, 0, dz/dxi)
            normal(0) = J(1, 0);
            normal(1) = -J(0, 0);
        } else if (edge == 1) {
            // eta increases from -1 to 1, xi = 1
            // n = (0, 1, 0) x (ds/deta, 0, dz/deta)
            normal(0) = J(1, 1);
            normal(1) = -J(0, 1);
        } else if (edge == 2) {
            // xi decreases from 1 to -1, eta = 1
            // n = -(0, 1, 0) x (ds/dxi, 0, dz/dxi)
            normal(0) = -J(1, 0);
            normal(1) = J(0, 0);
        } else {
            // eta decreases from 1 to -1, xi = -1
            // n = -(0, 1, 0) x (ds/deta, 0, dz/deta)
            normal(0) = -J(1, 1);
            normal(1) = J(0, 1);
        }
        return normal;
    }
    
    // get nodal coordinates
    const eigen::DMat24 &getNodalSZ() const {
        return mNodalSZ;
    }
    
    // get minimum edge length
    double getMinEdgeLength() const {
        return mMinEdgeLength;
    }
    
protected:
    // rotate CS such that mCurvedOuter lies on edge 2
    void rotateQ2(const eigen::DCol2 &xieta,
                  eigen::DMat24 &szQ2, eigen::DMat24 &rtQ2,
                  eigen::DCol2 &xietaQ2) const {
        // direct multiplication for sz and xieta
        szQ2 = sOrthogQ2[mCurvedOuter] * mNodalSZ;
        xietaQ2 = sOrthogQ2[mCurvedOuter] * xieta;
        
        // get r and theta
        // NOTE: 1) must use atan2 (instead of acos) because the rotation
        //          spans all four quadrants
        //       2) no need to check undefined outputs because spherical
        //          elements won't apprear near the planet center (sz=0)
        rtQ2 = geodesy::xy2sphi(szQ2, false, 0.);
        // need the angle from z-axis
        rtQ2.row(1) = numerical::dPi * .5 - rtQ2.row(1).array();
    }
    
    
    ////////////////////////////// data //////////////////////////////
protected:
    // nodal coordinates
    const eigen::DMat24 mNodalSZ;
    
    // minimum edge length
    double mMinEdgeLength = 0.;
    
    // curved outer
    int mCurvedOuter = -1;
    
    
    /////////////////////// static ///////////////////////
public:
    // cycle 0->1->2->3->0->1->2->3
    static int cycle4(int p) {
        if (p > 3) {
            return p - 4;
        }
        if (p < 0) {
            return p + 4;
        }
        return p;
    }
    
protected:
    // making two angles "numerically" close to each other
    // NOTE: this is essential because atan2 includes jumps while
    //       angles will be interpolated directly in mapping
    static void makeAnglesClose(double &a, double &b) {
        if (a - b > numerical::dPi) {
            b += 2. * numerical::dPi;
        }
        if (b - a > numerical::dPi) {
            a += 2. * numerical::dPi;
        }
    }
    
    // coordinate rotation matrix Q2
    inline static std::array<eigen::DMat22, 4> sOrthogQ2 = {
        (eigen::DMat22(2, 2) << -1., 0., 0., -1.).finished(),
        (eigen::DMat22(2, 2) << 0., -1., 1., 0.).finished(),
        (eigen::DMat22(2, 2) << 1., 0., 0., 1.).finished(),
        (eigen::DMat22(2, 2) << 0., 1., -1., 0.).finished()};
};

#endif /* Mapping_hpp */
