// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D geometric models

#ifndef Geometric3D_hpp
#define Geometric3D_hpp

#include "Model3D.hpp"

class Geometric3D : public Model3D {
  public:
  // constructor
  Geometric3D(const std::string& modelName) : Model3D(modelName) {
    // nothing
  }

  // destructor
  virtual ~Geometric3D() = default;

  // apply to Quad
  virtual void
  applyTo(std::vector<Quad>& quads) const;

  protected:
  // get undulation on an element
  virtual bool
  getUndulation(const axisem3d::eigen::DMatX3& spz,
      const axisem3d::eigen::DMat24& nodalSZ,
      axisem3d::eigen::DColX& undulation) const {
    // no element check by virtual
    return getUndulation(spz, undulation);
  }

  // set undulation to quad
  virtual void
  setUndulationToQuad(const axisem3d::eigen::DColX& undulation, Quad& quad) const;

  public:
  // get undulation on points
  virtual bool
  getUndulation(const axisem3d::eigen::DMatX3& spz, axisem3d::eigen::DColX& undulation) const = 0;

  ////////////////////////////// static //////////////////////////////
  public:
  // build from inparam
  static std::shared_ptr<const Geometric3D>
  buildInparam(const ExodusMesh& exodusMesh,
      const LocalMesh& localMesh,
      const std::string& modelName,
      const std::string& keyInparam);
};

#endif /* Geometric3D_hpp */
