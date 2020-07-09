//
//  eigen_station.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  eigen for station

#ifndef eigen_station_hpp
#define eigen_station_hpp

#include "eigen_tensor.hpp"
#include "eigen_element.hpp"

namespace eigen {
    // tensor
    typedef Eigen::Tensor<numerical::Real, 3, Eigen::RowMajor> RTensor3;
    typedef Eigen::Tensor<numerical::Real, 1, Eigen::RowMajor> RTensor1;
    typedef Eigen::array<Eigen::DenseIndex, 3> IArray3;
    typedef Eigen::array<Eigen::DenseIndex, 2> IArray2;
    typedef Eigen::array<Eigen::DenseIndex, 1> IArray1;
    
    // weights
    typedef Eigen::Matrix<numerical::Real, 1, Eigen::Dynamic> RRowX;
    typedef Eigen::Matrix<double, 1, nPEM> DRowN;
    
    // inplane interpolation
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CMatX1;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> CMatX3;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 6> CMatX6;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 9> CMatX9;
    
    // Fourier interpolation
    typedef Eigen::Matrix<numerical::Real, 1, 1> RRow1;
    typedef Eigen::Matrix<numerical::Real, 1, 3> RRow3;
    typedef Eigen::Matrix<numerical::Real, 1, 6> RRow6;
    typedef Eigen::Matrix<numerical::Real, 1, 9> RRow9;
    
    // buffers
    // alawys use double for time output
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CColX;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> RMatX1_RM;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, 3, Eigen::RowMajor> RMatX3_RM;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, 6, Eigen::RowMajor> RMatX6_RM;
    typedef Eigen::Matrix<Real, Eigen::Dynamic, 9, Eigen::RowMajor> RMatX9_RM;
}

#endif /* eigen_station_hpp */
