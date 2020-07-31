//
//  ElementIO.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  IO for station output

#include "ElementIO.hpp"
#include "io.hpp"
#include "mpi.hpp"

// initialize
void ElementIO::initialize(const std::string &groupName,
                           int numRecordSteps,
                           const std::vector<std::string> &channels,
                           int npnts, const std::vector<int> &naGrid,
                           const eigen::IMatX4_RM &elemNaInfo,
                           const eigen::DMatXX_RM &elemCoords) {
    // mpi globals
    int nelem = (int)elemNaInfo.rows();
    std::vector<int> nelemg;
    mpi::gather(nelem, nelemg, MPI_ALL);
    mNumElementsGlobal = std::accumulate(nelemg.begin(), nelemg.end(), 0);
    if (mNumElementsGlobal == 0) {
        // no station at all
        mRankWithMaxNumElements = -1;
    } else {
        mRankWithMaxNumElements =
        (int)(std::max_element(nelemg.begin(), nelemg.end()) - nelemg.begin());
    }
    
    // create group dir
    if (mpi::root()) {
        io::mkdir(io::gOutputDirectory + "/elements/" + groupName);
    }
    // must wait for mkdir
    mpi::barrier();
}
