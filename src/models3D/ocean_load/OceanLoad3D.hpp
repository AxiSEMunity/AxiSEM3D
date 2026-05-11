// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D ocean-load models

#ifndef OceanLoad3D_hpp
#define OceanLoad3D_hpp

#include "Model3D.hpp"

class OceanLoad3D : public Model3D {
  public:
  // constructor
  OceanLoad3D(const std::string& modelName) : Model3D(modelName) {
    // nothing
  }

  // destructor
  virtual ~OceanLoad3D() = default;

  // apply to Quad
  virtual void
  applyTo(std::vector<Quad>& quads) const;

  protected:
  // get sum(rho * depth)
  virtual bool
  getSumRhoDepth(const axisem3d::eigen::DMatX3& spz,
      const axisem3d::eigen::DMat24& nodalSZ,
      axisem3d::eigen::DColX& sumRhoDepth) const = 0;

  // set sum(rho * depth) to quad
  virtual void
  setSumRhoDepthToQuad(const axisem3d::eigen::DColX& sumRhoDepth, Quad& quad) const;

  ////////////////////////////// static //////////////////////////////
  public:
  // build from inparam
  static std::shared_ptr<const OceanLoad3D>
  buildInparam(const ExodusMesh& exodusMesh,
      const LocalMesh& localMesh,
      const std::string& modelName,
      const std::string& keyInparam);
};

#endif /* OceanLoad3D_hpp */
