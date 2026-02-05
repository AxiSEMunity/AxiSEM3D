//
//  Undulation.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  vertical undulation
//  generator of PRT in core

#ifndef Undulation_hpp
#define Undulation_hpp

// data
#include "PhysicalProperty.hpp"

// gradient
#include "SolverFFTW.hpp"
class Quad;
namespace eigen {
    using numerical::ComplexD;
    using spectral::nPED;
    using spectral::nPEM;
    typedef Eigen::Matrix<ComplexD, nPED, nPED, Eigen::RowMajor> ZMatPP_RM;
    typedef std::vector<std::array<ZMatPP_RM, 1>> vec_ar1_ZMatPP_RM;
    typedef std::vector<std::array<ZMatPP_RM, 3>> vec_ar3_ZMatPP_RM;
    typedef Eigen::Matrix<double, Eigen::Dynamic, nPEM * 3> DMatXN3;
}

// release
#include "PRT.hpp"

class Undulation {
public:
    // add undulation
    void addUndulation(const eigen::arN_DColX &und) {
        mDeltaZ.addGLL(und);
    }
    
    // get elemental
    eigen::DMatXN getElemental() const {
        return mDeltaZ.getElemental();
    }
    
    // get pointwise
    eigen::arN_DColX getPointwise() const {
        return mDeltaZ.getPointwise();
    }
    
    // finishing 3D properties
    void finishing3D() const;
    
    // finished 3D properties
    void finished3D(const Quad &myQuad);
    
    // get Jacobian for mass
    eigen::arN_DColX getMassJacobian(const eigen::DMat2N &sz) const;
    
    // create PRT
    std::unique_ptr<const PRT> createPRT(const eigen::DMat2N &sz) const;
    
    // compute 3D normal at a point
    eigen::DMatX3 computeNormal3D(const eigen::DCol2 &n1D,
                                  const eigen::DMat2N &sz, int ipnt) const;
    
    
    ///////////////////////// data /////////////////////////
private:
    // delta Z
    PhysicalProperty<spectral::nPEM> mDeltaZ;
    
    // gradient of delta Z
    eigen::DMatXN3 mDeltaZ_RTZ;
    
    
    ///////////////////////// static /////////////////////////
public:
    // finished 3D properties
    static void finished3D();
    
private:
    // static fft variables
    inline static SolverFFTW<double, spectral::nPEM> sFFT_N1;
    inline static SolverFFTW<double, spectral::nPEM * 3> sFFT_N3;
    inline static eigen::vec_ar1_ZMatPP_RM sDeltaZ_Fourier;
    inline static eigen::vec_ar3_ZMatPP_RM sDeltaZ_SPZ_Fourier;
};

#endif /* Undulation_hpp */
