// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  Exodus mesh created by salvus mesher

#ifndef ExodusMesh_hpp
#define ExodusMesh_hpp

#include "eigen_mesh.hpp"
#include <map>
#include <vector>

class NetCDF_Reader;
class NrField;

class ExodusMesh {
  public:
  /////////////////////////////////////////////
  /////////////////// build ///////////////////
  /////////////////////////////////////////////

  ///////////////// construction /////////////////
  // constructor
  ExodusMesh(const std::string& meshFile);

  // verbose
  std::string
  verbose() const;

  private:
  // global variables
  void
  readBcastGlobal(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // connectivity
  void
  readBcastConnectivity(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // coordinates
  void
  readBcastCoordinates(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // side sets
  void
  readBcastSideSets(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // elemental variables
  void
  readBcastElemental(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // radial variables
  void
  readBcastRadial(const NetCDF_Reader& reader, double& memSup, double& memAll);
  // ellipticity
  void
  readEllipticity(const NetCDF_Reader& reader, double& memSup, double& memAll);

  public:
  ///////////////// nr field /////////////////
  // Nr at nodes
  void
  formNrAtNodes(const NrField& nrField, bool boundByInplane, bool useLuckyNumbers);

  // verbose Nr
  std::string
  verboseNr(bool boundByInplane, bool useLuckyNumbers) const;

  ////////////////////////////////////////////
  /////////////////// data ///////////////////
  ////////////////////////////////////////////
  private:
  // file name
  const std::string mFileName;

  // global variables and records
  std::map<std::string, double> mGlobalVariables;
  std::map<std::string, std::string> mGlobalRecords;
  // mesh generation cmdline
  std::string mCmdMeshGen = "";

  // side sets
  std::string mKeyLeftSS, mKeyRightSS, mKeyBottomSS, mKeyTopSS;
  std::map<std::string, std::map<int, int>> mSideSets;

  // storage type
  bool mElementNodesStorage = false;

  // radial variables
  // NOTE: these are elemental variables depending ONLY on radius (depth)
  //       all material properties in Exodus are assumed to be radial
  std::vector<double> mRadialCoords;
  std::map<std::string, axisem3d::eigen::DColX> mRadialVariables;

  // discontinuities
  axisem3d::eigen::DColX mDiscontinuities = axisem3d::eigen::DColX(0);

  // ellipticity (root-only)
  axisem3d::eigen::DMatXX_RM mEllipticityCurve = axisem3d::eigen::DMatXX_RM(0, 0);

  // super-only
  // wrap over super-only variables for thread safety
  struct ExodusSuperOnly {
    // connectivity
    axisem3d::eigen::IMatX4_RM mConnectivity = axisem3d::eigen::IMatX4_RM(0, 4);

    // coords
    axisem3d::eigen::DMatX2_RM mNodalCoords = axisem3d::eigen::DMatX2_RM(0, 2);

    // element-wise variables
    axisem3d::eigen::IColX mGeometryType = axisem3d::eigen::IColX(0);
    axisem3d::eigen::IColX mIsElementFluid = axisem3d::eigen::IColX(0);

    // Nr at nodes
    axisem3d::eigen::IColX mNodalNr = axisem3d::eigen::IColX(0);
  } mSuperOnly__AccessOnlyBy__mySuperOnly;

  // set access to super-only
  ExodusSuperOnly&
  mySuperOnly();

  // get access to super-only
  const ExodusSuperOnly&
  mySuperOnly() const;

  ///////////////////////////////////////////////////////////
  /////////////////////////// get ///////////////////////////
  ///////////////////////////////////////////////////////////

  public:
  ///////////////// super-only properties /////////////////
  // number of nodes
  int
  getNumNodes() const {
    return (int)mySuperOnly().mNodalCoords.rows();
  }

  // number of quads
  int
  getNumQuads() const {
    return (int)mySuperOnly().mConnectivity.rows();
  }

  // coords
  const axisem3d::eigen::DMatX2_RM&
  getNodalCoords() const {
    return mySuperOnly().mNodalCoords;
  }

  // connectivity
  const axisem3d::eigen::IMatX4_RM&
  getConnectivity() const {
    return mySuperOnly().mConnectivity;
  }

  // geometry
  const axisem3d::eigen::IColX&
  getGeometryType() const {
    return mySuperOnly().mGeometryType;
  }

  // fluid
  const axisem3d::eigen::IColX&
  getIsElementFluid() const {
    return mySuperOnly().mIsElementFluid;
  }

  // Nr at nodes
  const axisem3d::eigen::IColX&
  getNrAtNodes() const {
    return mySuperOnly().mNodalNr;
  }

  ///////////////// non-super-only properties /////////////////
  // Cartesian
  bool
  isCartesian() const {
    return mGlobalRecords.at("crdsys") != "spherical";
  }

  // attenuation
  bool
  hasAttenuation() const {
    return mGlobalVariables.find("nr_lin_solids") != mGlobalVariables.end();
  }

  // isotropic
  bool
  isIsotropic() const {
    return mRadialVariables.find("VP") != mRadialVariables.end();
  }

  // global variables
  double
  getGlobalVariable(const std::string& key) const {
    return mGlobalVariables.at(key);
  }

  // discontinuities
  const axisem3d::eigen::DColX&
  getDiscontinuities() const {
    return mDiscontinuities;
  }

  // ellipticity
  const axisem3d::eigen::DMatXX_RM&
  getEllipticityCurve() const {
    return mEllipticityCurve;
  }

  // radial coords
  const std::vector<double>&
  getRadialCoords() const {
    return mRadialCoords;
  }

  // radial variables
  const std::map<std::string, axisem3d::eigen::DColX>&
  getRadialVariables() const {
    return mRadialVariables;
  }

  ///////////////// side methods /////////////////
  // get side
  int
  getSide(const std::string& ssName, int iquad) const {
    try {
      return mSideSets.at(ssName).at(iquad);
    } catch (...) {
      return -1;
    }
  }

  // left
  int
  getLeftSide(int iquad) const {
    return getSide(mKeyLeftSS, iquad);
  }

  // right
  int
  getRightSide(int iquad) const {
    return getSide(mKeyRightSS, iquad);
  }

  // bottom
  int
  getBottomSide(int iquad) const {
    return getSide(mKeyBottomSS, iquad);
  }

  // top
  int
  getTopSide(int iquad) const {
    return getSide(mKeyTopSS, iquad);
  }

  // axial
  int
  getAxialSide(int iquad) const {
    return getLeftSide(iquad);
  }

  // get boundary info for ABC
  // input: "RIGHT", "BOTTOM" or "TOP"
  // output: 1) whether mesh contains this boundary (return)
  //         2) location of this boundary (outer)
  //         3) mesh span in the normal direction of this boundary (meshSpan)
  bool
  boundaryInfoABC(const std::string& boundaryKey, double& outer, double& meshSpan) const;

  ///////////////// surface /////////////////
  // mesh surface
  double
  getMeshSurface() const {
    return mRadialCoords.back() + getGlobalVariable("dist_tolerance");
  }

  // solid surface
  double
  getSolidSurface() const {
    for (int irad = (int)mRadialCoords.size() - 1; irad > 0; irad--) {
      // get vs
      double vs =
          (isIsotropic() ? mRadialVariables.at("VS")(irad) : mRadialVariables.at("VSV")(irad));
      // solid
      if (vs > numerical::dEpsilon) {
        return (mRadialCoords[irad] + getGlobalVariable("dist_tolerance"));
      }
    }
    // pure fluid, return a level below mesh bottom
    return getMeshBottom() - getGlobalVariable("dist_tolerance") * 1000.;
  }

  // mesh bottom
  double
  getMeshBottom() const {
    return mRadialCoords.front() - getGlobalVariable("dist_tolerance");
  }
};

#endif /* ExodusMesh_hpp */
