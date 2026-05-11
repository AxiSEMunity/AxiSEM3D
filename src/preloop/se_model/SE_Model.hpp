// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  spectral element model
//  Exodus mesh ---
//                | -> mesh --------
//  local mesh ----                | -> SE model (this) -> domain
//                     3D models ---

#ifndef SE_Model_hpp
#define SE_Model_hpp

#include "eigen_sem.hpp"

// constructor
class ExodusMesh;
class ABC;
class LocalMesh;
class Model3D;

// components
#include "Quad.hpp"
#include "GLLPoint.hpp"

// release
class AttBuilder;
class Domain;

// source-receiver
class Geometric3D;
#include "RTreeND.hpp"

class SE_Model {
  public:
  // step 1: constructor
  SE_Model(const ExodusMesh& exodusMesh,
      const ABC& abc,
      const LocalMesh& localMesh,
      const std::vector<std::shared_ptr<const Model3D>>& models3D,
      bool useLuckyNumbers);

  // step 2: get dt for attenuation
  double
  computeDt(double courant, const ABC& abc, axisem3d::eigen::DCol2& sz) const;

  // step 3: release to domain
  void
  release(const ABC& abc,
      const LocalMesh& localMesh,
      const std::unique_ptr<const AttBuilder>& attBuilder,
      const TimeScheme& timeScheme,
      Domain& domain);

  // step 4: measure
  axisem3d::eigen::DColX
  measureCost(const ExodusMesh& exodusMesh,
      const LocalMesh& localMesh,
      const TimeScheme& timeScheme,
      bool measurePoint) const;

  // verbose
  std::string
  verbose(const std::string& titleSEM) const;

  ////////////////// source-receiver //////////////////
  // get total undulation
  double
  getTotalUndulation(const axisem3d::eigen::DRow3& spzRef) const;

  // undulated to reference
  axisem3d::eigen::DRow3
  undulatedToReference(const axisem3d::eigen::DRow3& spzUnd) const;

  // form inplane RTree
  void
  formInplaneRTree();

  // locate inplane
  int
  locateInplane(const axisem3d::eigen::DCol2& sz, bool inFluid) const;

  // compute inplane factor
  axisem3d::eigen::DRowN
  computeInplaneFactor(const axisem3d::eigen::DCol2& sz, int iquad) const;

  // get quads
  const std::vector<Quad>&
  getQuads() const {
    return mQuads;
  }

  ////////////////// wavefield scanning //////////////////
  // initialize scanning on points
  void
  initScanningOnPoints(bool vertexOnly) const;

  private:
  template <class Vec>
  static Vec
  interpLagrange(double target, const Vec& bases) {
    Vec weights = Vec::Ones(bases.rows(), bases.cols());
    for (int idgr = 0; idgr < bases.size(); idgr++) {
      for (int ibase = 0; ibase < bases.size(); ibase++) {
        if (ibase != idgr) {
          weights[idgr] *= target - bases[ibase];
          weights[idgr] /= bases[idgr] - bases[ibase];
        }
      }
    }
    return weights;
  }

  ////////////////////// data //////////////////////
  private:
  ///////// components /////////
  // quads
  std::vector<Quad> mQuads;
  // GLL points
  std::vector<GLLPoint> mGLLPoints;

  ///////// source-receiver /////////
  // distance tolerance of mesh
  const double mDistToleranceMesh;

  // 3D geometric models
  std::vector<std::shared_ptr<const Geometric3D>> mModelsG3D;

  // local mesh range
  axisem3d::eigen::DCol2 mRangeS = axisem3d::eigen::DCol2::Zero();
  axisem3d::eigen::DCol2 mRangeZ = axisem3d::eigen::DCol2::Zero();
  axisem3d::eigen::DCol2 mRangeR = axisem3d::eigen::DCol2::Zero();
  axisem3d::eigen::DCol2 mRangeT = axisem3d::eigen::DCol2::Zero();

  // inplane
  std::unique_ptr<RTreeND<2, 1, int>> mRTreeFluid = nullptr;
  std::unique_ptr<RTreeND<2, 1, int>> mRTreeSolid = nullptr;
};

#endif /* SE_Model_hpp */
