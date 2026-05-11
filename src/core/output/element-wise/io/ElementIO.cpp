// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  IO for station output

#include "ElementIO.hpp"
#include "io.hpp"
#include "mpi.hpp"

// initialize
void
ElementIO::initialize(const std::string& groupName,
    int numRecordSteps,
    const std::vector<std::string>& channels,
    int npnts,
    const std::vector<int>& naGrid,
    const axisem3d::eigen::IMatX4_RM& elemNaInfo,
    const axisem3d::eigen::DMatXX_RM& elemCoords) {
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
