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

// configuration for clamping and Gaussian smoothing of geometric models
struct ClampSmoothConfig {
  bool enabled = false;
  double clampMin = -1e10;
  double clampMax = 1e10;
  double sigmaDegreeOrMeter = 1.;
};

class Geometric3D : public Model3D {
  public:
  // constructor
  Geometric3D(
      const std::string& modelName, const ClampSmoothConfig& clampSmooth = ClampSmoothConfig()) :
      Model3D(modelName), mClampSmooth(clampSmooth) {
    // nothing
  }

  // destructor
  virtual ~Geometric3D() = default;

  // apply to Quad
  virtual void
  applyTo(std::vector<Quad>& quads) const;

  protected:
  // units of the two axes on the model's native horizontal grid
  enum class SmoothGridUnit { METER, RADIAN, DEGREE };

  // clamp in metres, then smooth on the model's native horizontal grid
  void
  applyClampSmooth(eigen::DMatXX& data,
      const std::array<double, 2>& gridSpacing,
      SmoothGridUnit gridUnit,
      const std::array<bool, 2>& periodicAxes = {false, false}) const;

  // verbose for clamp and smooth
  std::string
  verboseClampSmooth(int indent, int width) const;

  // get undulation on an element
  virtual bool
  getUndulation(
      const eigen::DMatX3& spz, const eigen::DMat24& nodalSZ, eigen::DColX& undulation) const {
    // no element check by virtual
    return getUndulation(spz, undulation);
  }

  // set undulation to quad
  virtual void
  setUndulationToQuad(const eigen::DColX& undulation, Quad& quad) const;

  public:
  // get undulation on points
  virtual bool
  getUndulation(const eigen::DMatX3& spz, eigen::DColX& undulation) const = 0;

  ////////////////////////////// static //////////////////////////////
  public:
  // build from inparam
  static std::shared_ptr<const Geometric3D>
  buildInparam(const ExodusMesh& exodusMesh,
      const LocalMesh& localMesh,
      const std::string& modelName,
      const std::string& keyInparam);

  protected:
  // common clamp-and-smooth configuration
  const ClampSmoothConfig mClampSmooth;
};

#endif /* Geometric3D_hpp */
