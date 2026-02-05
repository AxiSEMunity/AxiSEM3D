//
//  eigen_generic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/7/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  generic eigen matrix and vector types

#ifndef eigen_generic_hpp
#define eigen_generic_hpp

#include "eigen.hpp"
#include "numerical.hpp"

namespace eigen {
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using numerical::Real;
    using numerical::ComplexR;
    using numerical::ComplexD;
    
    // Matrix Nomenclature
    // DMatXN_RM => D           Mat      X          N          _RM
    //              ^           ^        ^          ^          ^
    //              datatype    shape    1st dim    2nd dim    storage order
    // datatype options:
    //     I => int
    //     R => Real
    //     D => double
    //     C => ComplexR
    //     Z => ComplexD
    // shape options:
    //     Mat => matrix
    //     Row => row vector
    //     Col => column vector
    // dimension options:
    //     P => nPntEdge
    //     N => nPntElem
    //     X => Eigen::Dynamic
    //     3 => 3 (constant numbers)
    // storage order options:
    //     _RM => Eigen::RowMajor
    //     _CM => Eigen::ColMajor
    //     default => Eigen::ColMajor
    
    //////////////////////////// int ////////////////////////////
    typedef Eigen::Matrix<int, Dynamic, Dynamic> IMatXX;
    typedef Eigen::Matrix<int, Dynamic, Dynamic, RowMajor> IMatXX_RM;
    typedef Eigen::Matrix<int, Dynamic, 1> IColX;
    typedef Eigen::Matrix<int, 1, Dynamic> IRowX;
    
    //////////////////////////// Real ////////////////////////////
    typedef Eigen::Matrix<Real, Dynamic, Dynamic> RMatXX;
    typedef Eigen::Matrix<Real, Dynamic, Dynamic, RowMajor> RMatXX_RM;
    typedef Eigen::Matrix<Real, Dynamic, 1> RColX;
    typedef Eigen::Matrix<Real, 1, Dynamic> RRowX;
    
    //////////////////////////// double ////////////////////////////
    typedef Eigen::Matrix<double, Dynamic, Dynamic> DMatXX;
    typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> DMatXX_RM;
    typedef Eigen::Matrix<double, Dynamic, 1> DColX;
    typedef Eigen::Matrix<double, 1, Dynamic> DRowX;
    
    //////////////////////////// Real complex ////////////////////////////
    typedef Eigen::Matrix<ComplexR, Dynamic, Dynamic> CMatXX;
    typedef Eigen::Matrix<ComplexR, Dynamic, Dynamic, RowMajor> CMatXX_RM;
    typedef Eigen::Matrix<ComplexR, Dynamic, 1> CColX;
    typedef Eigen::Matrix<ComplexR, 1, Dynamic> CRowX;
    
    //////////////////////////// double complex ////////////////////////////
    typedef Eigen::Matrix<ComplexD, Dynamic, Dynamic> ZMatXX;
    typedef Eigen::Matrix<ComplexD, Dynamic, Dynamic, RowMajor> ZMatXX_RM;
    typedef Eigen::Matrix<ComplexD, Dynamic, 1> ZColX;
    typedef Eigen::Matrix<ComplexD, 1, Dynamic> ZRowX;
}

#endif /* eigen_generic_hpp */
