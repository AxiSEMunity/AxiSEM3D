//
//  eigen_point.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 2/21/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  eigen typedef for point

#ifndef eigen_point_hpp
#define eigen_point_hpp

#include "eigen.hpp"
#include "numerical.hpp"

namespace eigen {
    // coords
    typedef Eigen::Matrix<double, 1, 2> DRow2;
    
    // fields on a fluid point
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DColX;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RColX;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CColX;
    
    // fields on a solid point
    typedef Eigen::Matrix<double, Eigen::Dynamic, 3> DMatX3;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 3> RMatX3;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> CMatX3;
}

#endif /* eigen_point_hpp */
