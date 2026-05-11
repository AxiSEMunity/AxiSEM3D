// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Ascii IO for station output
//  only for a small number of stations

#ifndef StationIO_Ascii_hpp
#define StationIO_Ascii_hpp

#include "StationIO.hpp"
#include <memory>

class StationIO_Ascii : public StationIO {
  public:
  // constructor
  StationIO_Ascii(bool stationCentric) : mStationCentric(stationCentric) {
    // nothing
  }

  // initialize
  void
  initialize(const std::string& groupName,
      int numRecordSteps,
      const std::vector<std::string>& channels,
      const std::vector<std::string>& stKeys);

  // finalize
  void
  finalize();

  // dump to file
  void
  dumpToFile(const axisem3d::eigen::DColX& bufferTime,
      const axisem3d::eigen::RTensor3& bufferFields,
      int bufferLine);

  private:
  // create file stream
  void
  createFileStream(const std::string& varKey, int precision);

  // station centric
  const bool mStationCentric;

  // file streams
  std::vector<std::unique_ptr<std::ofstream>> mFileStreams;
};

#endif /* StationIO_Ascii_hpp */
