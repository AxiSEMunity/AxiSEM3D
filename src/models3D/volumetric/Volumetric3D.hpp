// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D volumetric models
//  density, velocity, elasticity, attenuation

#ifndef Volumetric3D_hpp
#define Volumetric3D_hpp

#include "Model3D.hpp"
#include <map>

namespace axisem3d::eigen {
  // anisotropy
  typedef Eigen::Matrix<int, 6, 6> IMat66;
  typedef Eigen::Matrix<double, 6, 6> DMat66;
  typedef Eigen::Matrix<double, 3, 3> DMat33;
} // namespace axisem3d::eigen

class Volumetric3D : public Model3D {
  public:
  // reference kind
  enum class ReferenceKind { ABS, REF1D, REF3D, REF_PERTURB };
  inline static const std::map<ReferenceKind, std::string> sReferenceKindStr = {
      {ReferenceKind::ABS, "ABS"},
      {ReferenceKind::REF1D, "REF1D"},
      {ReferenceKind::REF3D, "REF3D"},
      {ReferenceKind::REF_PERTURB, "REF_PERTURB"}};

  // constructor
  Volumetric3D(const std::string& modelName) : Model3D(modelName) {
    // nothing
  }

  // destructor
  virtual ~Volumetric3D() = default;

  // apply to Quad
  virtual void
  applyTo(std::vector<Quad>& quads) const;

  protected:
  // using reference or undulated geometry
  virtual bool
  usingUndulatedGeometry() const = 0;

  // get property info
  virtual void
  getPropertyInfo(
      std::vector<std::string>& propKeys, std::vector<ReferenceKind>& refKinds) const = 0;

  // get properties
  virtual bool
  getProperties(const axisem3d::eigen::DMatX3& spz,
      const axisem3d::eigen::DMat24& nodalSZ,
      axisem3d::eigen::IMatXX& inScopes,
      axisem3d::eigen::DMatXX& propValues) const = 0;

  // set properties to quad
  virtual void
  setPropertiesToQuad(const std::vector<std::string>& propKeys,
      const std::vector<ReferenceKind>& refKinds,
      const axisem3d::eigen::IMatXX& inScopes,
      const axisem3d::eigen::DMatXX& propValues,
      Quad& quad) const;

  ////////////////////////////// static //////////////////////////////
  // Bond transformation for rotating Cijkl
  static void
  bondTransformation(const axisem3d::eigen::DMat66& inCijkl,
      double alpha,
      double beta,
      double gamma,
      axisem3d::eigen::DMat66& outCijkl);

  public:
  // build from inparam
  static std::shared_ptr<const Volumetric3D>
  buildInparam(const ExodusMesh& exodusMesh,
      const LocalMesh& localMesh,
      const std::string& modelName,
      const std::string& keyInparam);
};

#endif /* Volumetric3D_hpp */
