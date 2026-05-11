// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  local mesh: skeleton + element-GLL vicinity + mpi communication

#ifndef LocalMesh_hpp
#define LocalMesh_hpp

#include "eigen_mesh.hpp"
#include <vector>

class ExodusMesh;

class LocalMesh {
  // friends
  friend class SE_Model;
  friend class Quad;

  public:
  /////////////////////////// build //////////////////////////
  // constructor
  LocalMesh(const ExodusMesh& exodusMesh, const axisem3d::eigen::DColX& weights);

  // free memory after building SE_Model
  void
  freeMemorySE_ModelBuilt() {
    mConnectivity.resize(0, 4);
    mNodalCoords.resize(0, 2);
    mGeometryType.resize(0);
    mIsElementFluid.resize(0);
    mNodalNr.resize(0);
    mL2G_GLL.resize(0);
  }

  // free memory after SE_Model released to domain
  void
  freeMemorySE_ModelReleased() {
    mCommProc.clear();
    mCommMyGLL.clear();
  }

  private:
  // domain decomposition
  static void
  decomposeExodusMesh(const ExodusMesh& exodusMesh,
      const axisem3d::eigen::DColX& weights,
      axisem3d::eigen::IColX& elemRank);

  // build local skeleton
  void
  buildLocalSkeleton(const ExodusMesh& exodusMesh, const axisem3d::eigen::IColX& elemRank);

  // build element-GLL vicinity and MPI communication
  void
  buildElementGLL_CommMPI(const ExodusMesh& exodusMesh, const axisem3d::eigen::IColX& elemRank);

  /////////////////////////// info //////////////////////////
  public:
  // verbose
  std::string
  verbose(const std::string& meshTitle,
      const std::string& weightsKey,
      const axisem3d::eigen::DColX& weights) const;

  // plot domain decomposition
  void
  plotDD(const std::string& fname, const axisem3d::eigen::DColX& weights) const;

  /////////////////////////// data //////////////////////////
  private:
  // local skeleton (super-only in ExodusMesh)
  axisem3d::eigen::IColX mL2G_Element;
  axisem3d::eigen::IMatX4_RM mConnectivity;
  axisem3d::eigen::DMatX2_RM mNodalCoords;
  axisem3d::eigen::IColX mGeometryType;
  axisem3d::eigen::IColX mIsElementFluid;
  axisem3d::eigen::IColX mNodalNr;

  // element-GLL vicinity
  axisem3d::eigen::IColX mL2G_GLL;
  axisem3d::eigen::IMatXN_RM mElementGLL;

  // mpi communication
  std::vector<int> mCommProc;
  std::vector<std::vector<int>> mCommMyGLL;
};

#endif /* LocalMesh_hpp */
