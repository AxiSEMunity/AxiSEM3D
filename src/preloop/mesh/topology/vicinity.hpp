//
//  vicinity.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/19/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element-GLL vicinity

#ifndef vicinity_hpp
#define vicinity_hpp

#include "eigen_mesh.hpp"
#include <vector>

namespace vicinity {
    namespace constants {
#if _NPOL == 1
        // nPol = 1
        const std::vector<int> gNodeIPnt = {0, 2, 3, 1};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 2},
            {2, 3},
            {3, 1},
            {1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3};
#elif _NPOL == 2
        // nPol = 2
        const std::vector<int> gNodeIPnt = {0, 6, 8, 2};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 3, 6},
            {6, 7, 8},
            {8, 5, 2},
            {2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 5, 6, 7, 8};
#elif _NPOL == 3
        // nPol = 3
        const std::vector<int> gNodeIPnt = {0, 12, 15, 3};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 4, 8, 12},
            {12, 13, 14, 15},
            {15, 11, 7, 3},
            {3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
#elif _NPOL == 4
        // nPol = 4
        const std::vector<int> gNodeIPnt = {0, 20, 24, 4};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 5, 10, 15, 20},
            {20, 21, 22, 23, 24},
            {24, 19, 14, 9, 4},
            {4, 3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24};
#elif _NPOL == 5
        // nPol = 5
        const std::vector<int> gNodeIPnt = {0, 30, 35, 5};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 6, 12, 18, 24, 30},
            {30, 31, 32, 33, 34, 35},
            {35, 29, 23, 17, 11, 5},
            {5, 4, 3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 31, 32,
            33, 34, 35};
#elif _NPOL == 6
        // nPol = 6
        const std::vector<int> gNodeIPnt = {0, 42, 48, 6};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 7, 14, 21, 28, 35, 42},
            {42, 43, 44, 45, 46, 47, 48},
            {48, 41, 34, 27, 20, 13, 6},
            {6, 5, 4, 3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 5, 6, 7, 13, 14, 20, 21, 27, 28, 34, 35, 41, 42,
            43, 44, 45, 46, 47, 48};
#elif _NPOL == 7
        // nPol = 7
        const std::vector<int> gNodeIPnt = {0, 56, 63, 7};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 8, 16, 24, 32, 40, 48, 56},
            {56, 57, 58, 59, 60, 61, 62, 63},
            {63, 55, 47, 39, 31, 23, 15, 7},
            {7, 6, 5, 4, 3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 23, 24, 31, 32, 39, 40, 47,
            48, 55, 56, 57, 58, 59, 60, 61, 62, 63};
#elif _NPOL == 8
        // nPol = 8
        const std::vector<int> gNodeIPnt = {0, 72, 80, 8};
        const std::vector<std::vector<int>> gEdgeIPnt = {
            {0, 9, 18, 27, 36, 45, 54, 63, 72},
            {72, 73, 74, 75, 76, 77, 78, 79, 80},
            {80, 71, 62, 53, 44, 35, 26, 17, 8},
            {8, 7, 6, 5, 4, 3, 2, 1, 0}};
        const std::vector<int> gEdgeIPntAll = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17, 18, 26, 27, 35, 36, 44, 45,
            53, 54, 62, 63, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80};
#endif
    }
    
    // connectivity to element-GLL vicinity
    int connectivityToElementGLL(const eigen::IMatX4_RM &connectivity,
                                 eigen::IMatXN_RM &elementGLL,
                                 std::vector<eigen::IColX> &neighbours);
    
    // build local to global GLL
    eigen::IColX buildL2G_GLL(int nLocalGLL,
                              const eigen::IColX &myL2G_Element,
                              const eigen::IMatXN_RM &myElementGLL,
                              const eigen::IMatXN_RM &globalElementGLL);
    
    // build mpi communication
    // std::vector as the output to use mpi::scatter
    // rank0 gll0_0 gll0_1 ... -1 rank1 gll1_0 gll1_1 ... -1 rank2 ...
    std::vector<int>
    buildCommMPI(const eigen::IColX &myL2G_Element,
                 const std::vector<eigen::IColX> &globalNeighbours,
                 const eigen::IColX &elemRank, int myRank,
                 const eigen::IMatX4_RM &globalConnectivity,
                 const eigen::IMatXN_RM &myElementGLL,
                 const eigen::IColX &myL2G_GLL);
}

#endif /* vicinity_hpp */
