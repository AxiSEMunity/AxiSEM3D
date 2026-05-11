// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  IO for element output

#ifndef ElementIO_hpp
#define ElementIO_hpp

#include "eigen_element_op.hpp"

class ElementIO {
  public:
  // destructor
  virtual ~ElementIO() = default;

  // initialize
  virtual void
  initialize(const std::string& groupName,
      int numRecordSteps,
      const std::vector<std::string>& channels,
      int npnts,
      const std::vector<int>& naGrid,
      const axisem3d::eigen::IMatX4_RM& elemNaInfo,
      const axisem3d::eigen::DMatXX_RM& elemCoords);

  // finalize
  virtual void
  finalize() = 0;

  // dump to file
  virtual void
  dumpToFile(const axisem3d::eigen::DColX& bufferTime,
      const std::vector<axisem3d::eigen::RTensor5>& bufferFields,
      int bufferLine) = 0;

  // set flush
  void
  setFlush(bool flush) {
    mFlush = flush;
  }

  //////////////// mpi global properties ////////////////
  protected:
  // total number of elements
  int mNumElementsGlobal = 0;

  // mpi rank with maximum number of elements
  int mRankWithMaxNumElements = -1;

  // flush
  bool mFlush = false;
};

#endif /* ElementIO_hpp */
