//
//  eigen_mesh.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  eigen typedef for mesh

#ifndef eigen_mesh_hpp
#define eigen_mesh_hpp

#include "eigen_generic.hpp"
#include "spectral.hpp"

namespace eigen {
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    
    ///////// nodal-level /////////
    // connectivity
    typedef Eigen::Matrix<int, Dynamic, 4, RowMajor> IMatX4_RM;
    // nodal coords
    typedef Eigen::Matrix<double, Dynamic, 2, RowMajor> DMatX2_RM;
    
    ///////// GLL-level /////////
    // GLL tag on elements
    typedef Eigen::Matrix<int, Dynamic, spectral::nPEM, RowMajor> IMatXN_RM;
}

#endif /* eigen_mesh_hpp */
