//
//  eigen_element.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/20/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  eigen typedef for element

#ifndef eigen_element_hpp
#define eigen_element_hpp

#include "eigen.hpp"
#include "numerical.hpp"
#include "spectral.hpp"

namespace eigen {
    using numerical::ComplexR;
    using numerical::Real;
    using spectral::nPEM;
    using spectral::nPED;
    
    /////////////// Fourier space ///////////////
    // wavefields
    typedef Eigen::Matrix<ComplexR, nPED, nPED, Eigen::RowMajor> CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 1>> vec_ar1_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 3>> vec_ar3_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 6>> vec_ar6_CMatPP_RM;
    typedef std::vector<std::array<CMatPP_RM, 9>> vec_ar9_CMatPP_RM;
    
    // properties
    typedef Eigen::Matrix<double, nPED, nPED, Eigen::RowMajor> DMatPP_RM;
    typedef Eigen::Matrix<Real, nPED, nPED, Eigen::RowMajor> RMatPP_RM;
    
    /////////////// Cardinal space ///////////////
    // complex wavefields
    typedef Eigen::Matrix<ComplexR, Eigen::Dynamic, nPEM> CMatXN;
    typedef Eigen::Matrix<ComplexR, Eigen::Dynamic, nPEM * 3> CMatXN3;
    typedef Eigen::Matrix<ComplexR, Eigen::Dynamic, nPEM * 6> CMatXN6;
    typedef Eigen::Matrix<ComplexR, Eigen::Dynamic, nPEM * 9> CMatXN9;
    
    // real wavefields
    typedef Eigen::Matrix<Real, Eigen::Dynamic, nPEM> RMatXN;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, nPEM * 3> RMatXN3;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, nPEM * 6> RMatXN6;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, nPEM * 9> RMatXN9;
    
    // properties
    typedef Eigen::Matrix<double, Eigen::Dynamic, nPEM> DMatXN;
}

#endif /* eigen_element_hpp */
