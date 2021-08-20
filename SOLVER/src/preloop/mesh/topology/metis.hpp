//
//  metis.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  dual graph partioning by Metis
//  compatible with both 32-bit and 64-bit builds of Metis

#ifndef metis_hpp
#define metis_hpp

#include "eigen_mesh.hpp"

namespace metis {
    // form neighbourhood of connectivity
    void formNeighbourhood(const eigen::IMatX4_RM &connectivity, int ncommon,
                           std::vector<eigen::IColX> &neighbours);
    
    // domain decomposition
    double decompose(const eigen::IMatX4_RM &connectivity,
                     const eigen::IColX &solid_fluid,
                     const eigen::DColX &weights, int npart, int rseed,
                     eigen::IColX &elemRank);
}

#endif /* metis_hpp */
