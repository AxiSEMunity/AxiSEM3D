// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  material properties
//  generator of Acoustic and Elastic in core

#ifndef Material_hpp
#define Material_hpp

#include "PhysicalProperty.hpp"
#include <map>

// constructor
class ExodusMesh;

// 3D models
#include "Volumetric3D.hpp"

// release
#include "Acoustic.hpp"
#include "Elastic.hpp"
class AttBuilder;

class Material {
  public:
  ////////////////////////// input //////////////////////////
  // constructor
  Material(const ExodusMesh& exodusMesh, const axisem3d::eigen::DMat24& nodalSZ, bool axial);

  // add a 3D property
  void
  addProperty3D(const std::string& propKey,
      const Volumetric3D::ReferenceKind& refKind,
      const axisem3d::eigen::arN_IColX& inScope,
      const axisem3d::eigen::arN_DColX& propValue);

  // finished 3D properties
  void
  finished3D();

  ////////////////////////// output //////////////////////////
  // get maximum velocity for dt
  axisem3d::eigen::DMatXN
  getMaxVelocity() const;

  // get rho, vp, vs for Clayton
  void
  getPointwiseRhoVpVs(axisem3d::eigen::arN_DColX& rho,
      axisem3d::eigen::arN_DColX& vp,
      axisem3d::eigen::arN_DColX& vs) const;

  // get mass for GLL-point setup
  axisem3d::eigen::arN_DColX
  getMass(const axisem3d::eigen::DRowN& integralFactor,
      const axisem3d::eigen::arN_DColX& jacobianPRT,
      bool fluid) const;

  // create Acoustic
  std::unique_ptr<Acoustic>
  createAcoustic() const;

  // create Elastic
  std::unique_ptr<Elastic>
  createElastic(const std::unique_ptr<const AttBuilder>& attBuilder,
      const axisem3d::eigen::DRow4& weightsCG4) const;

  private:
  // anisotropic
  std::unique_ptr<Elastic>
  createAnisotropic(const std::unique_ptr<const AttBuilder>& attBuilder,
      const axisem3d::eigen::DRow4& weightsCG4) const;

  // transversely isotropic
  std::unique_ptr<Elastic>
  createTISO(const std::unique_ptr<const AttBuilder>& attBuilder,
      const axisem3d::eigen::DRow4& weightsCG4) const;

  // isotropic
  std::unique_ptr<Elastic>
  createIsotropic(const std::unique_ptr<const AttBuilder>& attBuilder,
      const axisem3d::eigen::DRow4& weightsCG4) const;

  ////////////////////// property //////////////////////
  private:
  // get property pointwise
  axisem3d::eigen::arN_DColX
  getPointwise(const std::string& key) const {
    return getProperty(key).getPointwise();
  }

  // get property elemental
  axisem3d::eigen::DMatXN
  getElemental(const std::string& key) const {
    return getProperty(key).getElemental();
  }

  // get property for get
  const NodalPhysicalProperty&
  getProperty(const std::string& key) const {
    try {
      return mProperties.at(key);
    } catch (...) {
      throw std::runtime_error("Material::getProperty || "
                               "Unacceptable property key: " +
          key +
          "||"
          "Rheology type: " +
          RheologyTypeStr.at(currentRheology()));
    }
  }

  // get property for set
  NodalPhysicalProperty&
  getProperty(const std::string& key) {
    try {
      return mProperties.at(key);
    } catch (...) {
      throw std::runtime_error("Material::getProperty || "
                               "Unacceptable property key: " +
          key +
          "||"
          "Rheology type: " +
          RheologyTypeStr.at(currentRheology()));
    }
  }

  ////////////////////// rheology //////////////////////
  // rheology type
  enum class RheologyType { FLUID, ISO, TISO, ANISO };
  inline static const std::map<RheologyType, std::string> RheologyTypeStr = {
      {RheologyType::FLUID, "Fluid"},
      {RheologyType::ISO, "Isotropic solid"},
      {RheologyType::TISO, "Transversely isotropic solid"},
      {RheologyType::ANISO, "Full anisotropic solid"}};

  // current rheology
  RheologyType
  currentRheology() const {
    if (mProperties.find("C11") != mProperties.end()) {
      return RheologyType::ANISO;
    } else if (mProperties.find("VPV") != mProperties.end()) {
      return RheologyType::TISO;
    } else if (mProperties.find("VS") != mProperties.end()) {
      return RheologyType::ISO;
    } else {
      return RheologyType::FLUID;
    }
  }

  // evolve from ISO to TISO
  void
  evolveISO_TISO() {
    // copy vp, vs
    mProperties.insert({"VPV", mProperties.at("VP")});
    mProperties.insert({"VPH", mProperties.at("VP")});
    mProperties.insert({"VSV", mProperties.at("VS")});
    mProperties.insert({"VSH", mProperties.at("VS")});
    // eta
    mProperties.insert({"ETA",
        NodalPhysicalProperty(axisem3d::eigen::DRow4::Ones(), mProperties.at("RHO").axial())});
    // erase vp, vs
    mProperties.erase("VP");
    mProperties.erase("VS");
  }

  // evolve from TISO to ANISO
  void
  evolveTISO_ANISO() {
    // density and velocity
    const NodalPhysicalProperty& rho = mProperties.at("RHO");
    const NodalPhysicalProperty& vpv = mProperties.at("VPV");
    const NodalPhysicalProperty& vph = mProperties.at("VPH");
    const NodalPhysicalProperty& vsv = mProperties.at("VSV");
    const NodalPhysicalProperty& vsh = mProperties.at("VSH");
    const NodalPhysicalProperty& eta = mProperties.at("ETA");
    // A, C, F, L, N
    const NodalPhysicalProperty& A = rho * vph.pow(2.);
    const NodalPhysicalProperty& C = rho * vpv.pow(2.);
    const NodalPhysicalProperty& L = rho * vsv.pow(2.);
    const NodalPhysicalProperty& N = rho * vsh.pow(2.);
    const NodalPhysicalProperty& F = eta * (A - L * 2.);
    const NodalPhysicalProperty& zero =
        NodalPhysicalProperty(axisem3d::eigen::DRow4::Zero(), mProperties.at("RHO").axial());
    // Cijkl
    // non-zero
    mProperties.insert({"C11", A});
    mProperties.insert({"C22", A});
    mProperties.insert({"C33", C});
    mProperties.insert({"C44", L});
    mProperties.insert({"C55", L});
    mProperties.insert({"C66", N});
    mProperties.insert({"C12", A - N * 2.});
    mProperties.insert({"C13", F});
    mProperties.insert({"C23", F});
    // zero
    mProperties.insert({"C14", zero});
    mProperties.insert({"C15", zero});
    mProperties.insert({"C16", zero});
    mProperties.insert({"C24", zero});
    mProperties.insert({"C25", zero});
    mProperties.insert({"C26", zero});
    mProperties.insert({"C34", zero});
    mProperties.insert({"C35", zero});
    mProperties.insert({"C36", zero});
    mProperties.insert({"C45", zero});
    mProperties.insert({"C46", zero});
    mProperties.insert({"C56", zero});
    // erase
    mProperties.erase("VPV");
    mProperties.erase("VPH");
    mProperties.erase("VSV");
    mProperties.erase("VSH");
    mProperties.erase("ETA");
  }

  /////////////////////////// data ///////////////////////////
  private:
  // properties
  std::map<std::string, NodalPhysicalProperty> mProperties;
};

#endif /* Material_hpp */
