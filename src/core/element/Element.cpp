// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  spectral element

#include "Element.hpp"

// point
#include "Point.hpp"

// transform
#include "CoordTransformCartesian.hpp"
#include "CoordTransformSpherical.hpp"
#include "geodesy.hpp"

// side
#include "vicinity.hpp"

using spectral::nPEM;

/////////////////////////// point ///////////////////////////
// point set
void
Element::pointSet(bool elemInFourier) {
  // order
  mNr = 0;
  for (int ipnt = 0; ipnt < nPEM; ipnt++) {
    mNr = std::max(getPoint(ipnt).getNr(), mNr);
  }
  mNu_1 = mNr / 2 + 1;

  // check compatibility
  mGradQuad->checkCompatibility(mNr);
  if (mPRT) {
    mPRT->checkCompatibility(mNr, elemInFourier);
    // set to RTZ for PRT
    setToRTZ_ByMaterialOrPRT();
  }
}

// find boundary points by tag
std::vector<int>
Element::findBoundaryPointsByTag(const std::vector<int>& boundaryMeshTags) const {
  // point indices on an element boundary
  const std::vector<int>& allSidePoints = vicinity::constants::gEdgeIPntAll;

  // result
  std::vector<int> pointsFound;
  pointsFound.reserve(allSidePoints.size());

  // check mesh tag
  for (int ipnt : allSidePoints) {
    if (vector_tools::findSortedUnique(boundaryMeshTags, getPoint(ipnt).getMeshTag())) {
      // found a point on the boundary
      pointsFound.push_back(ipnt);
    }
  }

  // already sorted and unique
  return pointsFound;
}

// find boundary points by crds
std::vector<int>
Element::findBoundaryPointsByCrds(const std::vector<double>& boundaryCrdsRorZ,
    const std::vector<double>& boundaryCrdsTorS,
    double distTol) const {
  // point indices on an element boundary
  const std::vector<int>& allSidePoints = vicinity::constants::gEdgeIPntAll;

  // collect coords
  axisem3d::eigen::DColX crdRorZ(nPEM);
  axisem3d::eigen::DColX crdTorS(nPEM);
  for (int ipnt = 0; ipnt < nPEM; ipnt++) {
    const axisem3d::eigen::DRow2& sz = getPoint(ipnt).getCoords();
    if (geodesy::isCartesian()) {
      crdRorZ(ipnt) = sz(1);
      crdTorS(ipnt) = sz(0);
    } else {
      crdRorZ(ipnt) = sz.norm();
      crdTorS(ipnt) = (crdRorZ(ipnt) < numerical::dEpsilon) ? 0. : acos(sz(1) / crdRorZ(ipnt));
    }
  }

  // result
  std::vector<int> pointsFound;

  // check r or z
  for (double bRorZ : boundaryCrdsRorZ) {
    for (int ipnt : allSidePoints) {
      if (std::abs(crdRorZ[ipnt] - bRorZ) < distTol) {
        pointsFound.push_back(ipnt);
      }
    }
  }

  // check theta or s
  for (double bTorS : boundaryCrdsTorS) {
    for (int ipnt : allSidePoints) {
      double dist = std::abs(crdTorS[ipnt] - bTorS);
      if (!geodesy::isCartesian()) {
        // change angular distance to meter
        dist *= crdRorZ(ipnt);
      }
      if (dist < distTol) {
        pointsFound.push_back(ipnt);
      }
    }
  }

  // sort and unique points
  vector_tools::sortUnique(pointsFound);
  return pointsFound;
}

/////////////////////////// crd transform ///////////////////////////
// set to RTZ by material or PRT
void
Element::setToRTZ_ByMaterialOrPRT() {
  // coordinate transform
  createCoordTransform();
}

// create coordinate transform
// internally called by setToRTZ_ByMaterialOrPRT
// externally called by prepareWavefieldOutput
void
Element::createCoordTransform() {
  if (!mTransform) {
    // difference between Cartesian and spherical
    if (geodesy::isCartesian()) {
      mTransform = std::make_unique<const CoordTransformCartesian>();
    } else {
      // compute theta
      axisem3d::eigen::DMatPP_RM theta;
      for (int ipol = 0; ipol < spectral::nPED; ipol++) {
        for (int jpol = 0; jpol < spectral::nPED; jpol++) {
          int ipnt = ipol * spectral::nPED + jpol;
          const axisem3d::eigen::DRow2& sz = getPoint(ipnt).getCoords();
          const axisem3d::eigen::DRow2& rt = geodesy::sz2rtheta(sz, true);
          theta(ipol, jpol) = rt(1);
        }
      }
      mTransform = std::make_unique<const CoordTransformSpherical>(theta);
    }
  }
}
