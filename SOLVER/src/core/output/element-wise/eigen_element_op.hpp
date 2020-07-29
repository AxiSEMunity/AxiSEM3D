//
//  eigen_element_op.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 28/7/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  eigen for element output

#ifndef eigen_element_op_hpp
#define eigen_element_op_hpp

#include "eigen_station.hpp"
#include "eigen_generic.hpp"

namespace eigen {
    // tensor
    typedef Eigen::Tensor<numerical::Real, 5, Eigen::RowMajor> RTensor5;
    typedef Eigen::Tensor<numerical::Real, 4, Eigen::RowMajor> RTensor4;
    typedef Eigen::array<Eigen::DenseIndex, 5> IArray5;
    typedef Eigen::array<Eigen::DenseIndex, 4> IArray4;
    
    // element-na info
    typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> IMatX4_RM;
}

#endif /* eigen_element_op_hpp */
