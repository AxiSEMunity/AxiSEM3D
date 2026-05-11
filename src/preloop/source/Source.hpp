// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  source: generator of ElementSource in core

#ifndef Source_hpp
#define Source_hpp

class SE_Model;
class Domain;
class STF;
class Mechanism;
#include <memory>

#include "eigen.hpp"
namespace axisem3d::eigen {
  typedef Eigen::Matrix<double, 1, 3> DRow3;
  typedef Eigen::Matrix<double, 3, 3> DMat33;
} // namespace axisem3d::eigen

class Source {
  public:
  // constructor
  Source(const std::string& name,
      bool axial,
      bool sourceCentered,
      bool ellipticity,
      bool useDepth,
      bool depthSolid,
      bool undulatedGeometry,
      const axisem3d::eigen::DRow3& crdIn) :
      mName(name), mAxial(axial), mSourceCentered(sourceCentered), mEllipticity(ellipticity),
      mUseDepth(useDepth), mDepthSolid(depthSolid), mUndulatedGeometry(undulatedGeometry),
      mCrdIn(crdIn) {
    // nothing
  }

  private:
  // build from inparam
  static std::shared_ptr<const Source>
  buildInparam(int sindex);

  public:
  // compute spz
  static axisem3d::eigen::DRow3
  computeSPZ(const SE_Model& sem,
      const axisem3d::eigen::DRow3& crdIn,
      bool sourceCentered,
      bool xy,
      bool ellipticity,
      bool useDepth,
      bool depthSolid,
      bool undulatedGeometry,
      const std::string& errInfo,
      bool enforceOnAxis);

  private:
  // compute rotation matrix Q from input to (z, s, phi)
  const axisem3d::eigen::DMat33&
  computeQzsp(const axisem3d::eigen::DRow3& spz, bool ellipticity) const;

  // verbose
  std::string
  verbose(int sindex, const STF& stf, const Mechanism& mechanism) const;

  public:
  // release sources to domain
  static void
  release(const SE_Model& sem, Domain& domain, double dt, double& minT0);

  private:
  // data from inparam
  const std::string mName;
  const bool mAxial;
  const bool mSourceCentered;
  const bool mEllipticity;
  const bool mUseDepth;
  const bool mDepthSolid;
  const bool mUndulatedGeometry;
  const axisem3d::eigen::DRow3 mCrdIn;
};

#endif /* Source_hpp */
