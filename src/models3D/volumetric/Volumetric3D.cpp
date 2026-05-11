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

#include "Volumetric3D.hpp"
#include "Quad.hpp"
#include "mpi.hpp"

// apply to Quad
void
Volumetric3D::applyTo(std::vector<Quad>& quads) const {
  // info
  std::vector<std::string> propKeys;
  std::vector<ReferenceKind> refKinds;
  getPropertyInfo(propKeys, refKinds);

  // data
  if (!isSuperOnly()) {
    for (Quad& quad : quads) {
      // cardinal coordinates
      const axisem3d::eigen::DMatX3& spz = computeElemSPZ(quad, usingUndulatedGeometry());
      // compute values
      axisem3d::eigen::IMatXX inScopes;
      axisem3d::eigen::DMatXX propValues;
      bool elemInScope = getProperties(spz, quad.getNodalSZ(), inScopes, propValues);
      // set values to quad
      if (elemInScope) {
        setPropertiesToQuad(propKeys, refKinds, inScopes, propValues, quad);
      }
    }
  } else {
    mpi::enterInfer();
    for (int irank = 0; irank < mpi::nproc(); irank++) {
      // step 1: gather coords on infer and send to super
      std::vector<axisem3d::eigen::DMatX3> spzAll;
      std::vector<axisem3d::eigen::DMat24> szAll;
      if (irank == mpi::rank()) {
        // gather coords
        spzAll.reserve(quads.size());
        szAll.reserve(quads.size());
        for (Quad& quad : quads) {
          spzAll.push_back(computeElemSPZ(quad, usingUndulatedGeometry()));
          szAll.push_back(quad.getNodalSZ());
        }
        // send coords to super
        mpi::sendVecEigen(0, spzAll, 0);
        mpi::sendVecEigen(0, szAll, 1);
      }

      // step 2: compute values on super and send back to infer
      std::vector<axisem3d::eigen::IMatXX> inScopesAll;
      std::vector<axisem3d::eigen::DMatXX> propValuesAll;
      std::vector<axisem3d::eigen::IColX> elemInScopeAll;
      if (mpi::root()) {
        // recv coords from infer
        mpi::recvVecEigen(irank, spzAll, 0);
        mpi::recvVecEigen(irank, szAll, 1);
        // allocate values
        int nQuad = (int)spzAll.size();
        inScopesAll.reserve(nQuad);
        propValuesAll.reserve(nQuad);
        elemInScopeAll.push_back(axisem3d::eigen::IColX::Zero(nQuad));
        // compute values
        for (int iq = 0; iq < nQuad; iq++) {
          axisem3d::eigen::IMatXX inScopes;
          axisem3d::eigen::DMatXX propValues;
          elemInScopeAll[0](iq) = getProperties(spzAll[iq], szAll[iq], inScopes, propValues);
          inScopesAll.push_back(inScopes);
          propValuesAll.push_back(propValues);
        }
        // send values to infer
        mpi::sendVecEigen(irank, inScopesAll, 0);
        mpi::sendVecEigen(irank, propValuesAll, 1);
        mpi::sendVecEigen(irank, elemInScopeAll, 2);
      }

      // step 3: set values to quads on infer
      if (irank == mpi::rank()) {
        // recv values from super
        mpi::recvVecEigen(0, inScopesAll, 0);
        mpi::recvVecEigen(0, propValuesAll, 1);
        mpi::recvVecEigen(0, elemInScopeAll, 2);
        // set values to quads
        for (int iquad = 0; iquad < spzAll.size(); iquad++) {
          if (elemInScopeAll[0](iquad)) {
            setPropertiesToQuad(
                propKeys, refKinds, inScopesAll[iquad], propValuesAll[iquad], quads[iquad]);
          }
        }
      }
      // do irank one by one
      mpi::barrier();
    }
    mpi::enterWorld();
  }
}

// set properties to quad
void
Volumetric3D::setPropertiesToQuad(const std::vector<std::string>& propKeys,
    const std::vector<ReferenceKind>& refKinds,
    const axisem3d::eigen::IMatXX& inScopes,
    const axisem3d::eigen::DMatXX& propValues,
    Quad& quad) const {
  // property loop
  const axisem3d::eigen::IRowN& pointNr = quad.getPointNr();
  int nprop = (int)propKeys.size();
  for (int iprop = 0; iprop < nprop; iprop++) {
    // flattened to structured
    axisem3d::eigen::arN_IColX inScope;
    axisem3d::eigen::arN_DColX propValue;
    int row = 0;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
      int nr = pointNr(ipnt);
      inScope[ipnt] = inScopes.block(row, iprop, nr, 1);
      propValue[ipnt] = propValues.block(row, iprop, nr, 1);
      row += nr;
    }
    // set to material
    quad.getMaterialPtr()->addProperty3D(propKeys[iprop], refKinds[iprop], inScope, propValue);
  }
}

// Bond transformation for rotating Cijkl
void
Volumetric3D::bondTransformation(const axisem3d::eigen::DMat66& inCijkl,
    double alpha,
    double beta,
    double gamma,
    axisem3d::eigen::DMat66& outCijkl) {
  // R
  static axisem3d::eigen::DMat33 R1, R2, R3, R;
  R1 << 1., 0., 0., 0., cos(alpha), sin(alpha), 0., -sin(alpha), cos(alpha);
  R2 << cos(beta), 0., sin(beta), 0., 1., 0., -sin(beta), 0, cos(beta);
  R3 << cos(gamma), sin(gamma), 0., -sin(gamma), cos(gamma), 0., 0., 0., 1.;
  R = R1 * R2 * R3;

  // K
  static axisem3d::eigen::DMat33 K1, K2, K3, K4;
  K1.array() = R.array().pow(2.);
  K2 << R(0, 1) * R(0, 2), R(0, 2) * R(0, 0), R(0, 0) * R(0, 1), R(1, 1) * R(1, 2),
      R(1, 2) * R(1, 0), R(1, 0) * R(1, 1), R(2, 1) * R(2, 2), R(2, 2) * R(2, 0), R(2, 0) * R(2, 1);
  K3 << R(1, 0) * R(2, 0), R(1, 1) * R(2, 1), R(1, 2) * R(2, 2), R(2, 0) * R(0, 0),
      R(2, 1) * R(0, 1), R(2, 2) * R(0, 2), R(0, 0) * R(1, 0), R(0, 1) * R(1, 1), R(0, 2) * R(1, 2);
  K4 << R(1, 1) * R(2, 2) + R(1, 2) * R(2, 1), R(1, 2) * R(2, 0) + R(1, 0) * R(2, 2),
      R(1, 0) * R(2, 1) + R(1, 1) * R(2, 0), R(2, 1) * R(0, 2) + R(2, 2) * R(0, 1),
      R(2, 2) * R(0, 0) + R(2, 0) * R(0, 2), R(2, 0) * R(0, 1) + R(2, 1) * R(0, 0),
      R(0, 1) * R(1, 2) + R(0, 2) * R(1, 1), R(0, 2) * R(1, 0) + R(0, 0) * R(1, 2),
      R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0);
  static axisem3d::eigen::DMat66 K;
  K.block(0, 0, 3, 3) = K1;
  K.block(0, 3, 3, 3) = 2. * K2;
  K.block(3, 0, 3, 3) = K3;
  K.block(3, 3, 3, 3) = K4;

  // rotation
  outCijkl = K * inCijkl * K.transpose();
}

#include "StructuredGridV3D.hpp"
#include "sg_tools.hpp"

// build from inparam
std::shared_ptr<const Volumetric3D>
Volumetric3D::buildInparam(const ExodusMesh& exodusMesh,
    const LocalMesh& localMesh,
    const std::string& modelName,
    const std::string& keyInparam) {
  // short alias
  const InparamYAML& gm = inparam::gInparamModel;
  const std::string& root = keyInparam;

  // class name
  const std::string& className = gm.get<std::string>(root + ":class_name");

  // init class
  if (className == "StructuredGridV3D") {
    // file name
    const std::string& fname = gm.get<std::string>(root + ":nc_data_file");

    ////////////// coords //////////////
    const std::string& rootc = root + ":coordinates";
    // horizontal
    bool sourceCentered = false, xy = false, ellipticity = false;
    sg_tools::inparamHorizontal(gm, rootc, modelName, className, sourceCentered, xy, ellipticity);
    // vertical
    bool useDepth = false, depthSolid = false;
    sg_tools::inparamVertical(gm, rootc, modelName, className, useDepth, depthSolid);
    // variables
    std::array<std::string, 3> crdVarNames;
    std::array<int, 3> shuffleData;
    sg_tools::inparamVarRank<3>(gm, rootc, modelName, className, crdVarNames, shuffleData);
    // units
    double lengthUnit = 1., angleUnit = 1.;
    sg_tools::inparamUnits(gm, rootc, xy, lengthUnit, angleUnit);
    // other options
    bool undulated = gm.get<bool>(rootc + ":undulated_geometry");
    bool center = gm.get<bool>(rootc + ":whole_element_inplane");

    ////////////// properties //////////////
    // size
    const std::string& rootp = root + ":properties";
    int nprop = gm.get<int>(rootp + ":[?]");
    // info
    std::vector<std::tuple<std::string, std::string, double, ReferenceKind>> propertyInfo;
    for (int iprop = 0; iprop < nprop; iprop++) {
      // get key
      const std::string& key = gm.get<std::string>(rootp + ":{" + bstring::toString(iprop) + "}");
      const std::string& rootpi = rootp + ":[" + bstring::toString(iprop) + "]:" + key;
      // var, factor, ref
      const std::string& vname = gm.get<std::string>(rootpi + ":nc_var");
      double factor = gm.get<double>(rootpi + ":factor");
      ReferenceKind ref = gm.getWithLimits<ReferenceKind>(rootpi + ":reference_kind",
          {{"ABS", ReferenceKind::ABS},
              {"REF1D", ReferenceKind::REF1D},
              {"REF3D", ReferenceKind::REF3D},
              {"REF_PERTURB", ReferenceKind::REF_PERTURB}});
      propertyInfo.push_back({key, vname, factor, ref});
    }

    // super-only
    bool superOnly = gm.get<bool>(root + ":store_grid_only_on_leaders");

    // construct
    return std::make_shared<const StructuredGridV3D>(modelName,
        fname,
        crdVarNames,
        shuffleData,
        sourceCentered,
        xy,
        ellipticity,
        useDepth,
        depthSolid,
        undulated,
        lengthUnit,
        angleUnit,
        center,
        propertyInfo,
        superOnly);
  } else {
    // other models
  }

  // unknown class
  return nullptr;
}
