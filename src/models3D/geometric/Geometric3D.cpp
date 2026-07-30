// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: MIT
// Copyright (C) 2019 - 2026 by the AxiSEM3D authors
//
// This file is part of the AxiSEM3D library. See the LICENSE file for details.
//
// -----------------------------------------------------------------------------

//  3D geometric models

#include "Geometric3D.hpp"
#include "Quad.hpp"
#include "mpi.hpp"
#include "bstring.hpp"
#include "inparam.hpp"
#include <cmath>

namespace {
  // skimage/scipy default; intentionally not exposed in YAML
  constexpr double sGaussianTruncate = 4.;

  ClampSmoothConfig
  readClampSmooth(const InparamYAML& gm,
      const std::string& root,
      const std::string& modelName,
      const std::string& className) {
    ClampSmoothConfig config;
    // Crust1 keeps smoothing on by default; other geometric models opt in
    config.enabled = className == "Crust1G3D";

    const std::string roots = root + ":clamp_smooth";
    if (gm.contains(roots)) {
      if (gm.contains(roots + ":enabled")) {
        config.enabled = gm.get<bool>(roots + ":enabled");
      }
      if (gm.contains(roots + ":clamp_range")) {
        const std::vector<double>& range = gm.getVector<double>(roots + ":clamp_range");
        if (range.size() != 2) {
          throw std::runtime_error("Geometric3D::buildInparam || "
                                   "clamp_smooth:clamp_range must contain two values. || "
                                   "Model name: " +
              modelName + " || Class name: " + className);
        }
        config.clampMin = range[0];
        config.clampMax = range[1];
      }
      if (gm.contains(roots + ":sigma_degree_or_meter")) {
        config.sigmaDegreeOrMeter = gm.get<double>(roots + ":sigma_degree_or_meter");
      }
    }

    if (config.clampMin > config.clampMax) {
      throw std::runtime_error("Geometric3D::buildInparam || "
                               "clamp_smooth:clamp_range must be ascending. || "
                               "Model name: " +
          modelName + " || Class name: " + className);
    }
    if (config.enabled &&
        (!(config.sigmaDegreeOrMeter > 0.) || !std::isfinite(config.sigmaDegreeOrMeter))) {
      throw std::runtime_error("Geometric3D::buildInparam || "
                               "clamp_smooth:sigma_degree_or_meter must be positive and finite. || "
                               "Model name: " +
          modelName + " || Class name: " + className);
    }
    return config;
  }
} // namespace

// clamp in metres, then smooth on the model's native horizontal grid
void
Geometric3D::applyClampSmooth(eigen::DMatXX& data,
    const std::array<double, 2>& gridSpacing,
    SmoothGridUnit gridUnit,
    const std::array<bool, 2>& periodicAxes) const {
  if (!mClampSmooth.enabled || data.size() == 0) {
    return;
  }

  // clamp always acts on the final physical undulation in metres
  data = data.array().max(mClampSmooth.clampMin).min(mClampSmooth.clampMax);

  // user sigma is in degrees on spherical grids and metres on Cartesian grids;
  // source-centred spherical coordinates are stored internally in radians
  double sigmaGridUnit = mClampSmooth.sigmaDegreeOrMeter;
  if (gridUnit == SmoothGridUnit::RADIAN) {
    sigmaGridUnit *= numerical::dDegree;
  }

  // separable Gaussian with model-specific boundary topology
  for (int axis = 0; axis < 2; axis++) {
    if (!(gridSpacing[axis] > 0.) || !std::isfinite(gridSpacing[axis])) {
      throw std::runtime_error("Geometric3D::applyClampSmooth || "
                               "Invalid spacing on the native model grid. || Model name: " +
          mModelName);
    }
    const double sigmaIndex = sigmaGridUnit / gridSpacing[axis];
    const int radius = (int)std::round(sGaussianTruncate * sigmaIndex);
    if (radius == 0) {
      continue;
    }

    eigen::DColX kernel(radius * 2 + 1);
    for (int iker = -radius; iker <= radius; iker++) {
      kernel(iker + radius) = std::exp(-.5 * iker * iker / (sigmaIndex * sigmaIndex));
    }
    kernel /= kernel.sum();

    const int nAxis = axis == 0 ? (int)data.rows() : (int)data.cols();
    eigen::DMatXX result = eigen::DMatXX::Zero(data.rows(), data.cols());
    for (int irow = 0; irow < data.rows(); irow++) {
      for (int icol = 0; icol < data.cols(); icol++) {
        for (int iker = -radius; iker <= radius; iker++) {
          int index = (axis == 0 ? irow : icol) + iker;
          if (periodicAxes[axis]) {
            index %= nAxis;
            if (index < 0) {
              index += nAxis;
            }
          } else {
            index = std::max(0, std::min(index, nAxis - 1));
          }
          const int jrow = axis == 0 ? index : irow;
          const int jcol = axis == 1 ? index : icol;
          result(irow, icol) += kernel(iker + radius) * data(jrow, jcol);
        }
      }
    }
    data = result;
  }
}

// verbose for clamp and smooth
std::string
Geometric3D::verboseClampSmooth(int indent, int width) const {
  using namespace bstring;
  std::stringstream ss;
  ss << boxSubTitle(indent, "Clamp and Gaussian smoothing");
  ss << boxEquals(indent + 2, width, "enabled", mClampSmooth.enabled);
  ss << boxEquals(
      indent + 2, width, "clamp range (m)", range(mClampSmooth.clampMin, mClampSmooth.clampMax));
  ss << boxEquals(indent + 2, width, "sigma (degree or m)", mClampSmooth.sigmaDegreeOrMeter);
  return ss.str();
}

// apply to Quad
void
Geometric3D::applyTo(std::vector<Quad>& quads) const {
  if (!isSuperOnly()) {
    for (Quad& quad : quads) {
      // cardinal coordinates
      const eigen::DMatX3& spz = computeElemSPZ(quad);
      // compute values
      eigen::DColX und;
      bool elemInScope = getUndulation(spz, quad.getNodalSZ(), und);
      // set values to quad
      if (elemInScope) {
        setUndulationToQuad(und, quad);
      }
    }
  } else {
    mpi::enterInfer();
    for (int irank = 0; irank < mpi::nproc(); irank++) {
      // step 1: gather coords on infer and send to super
      std::vector<eigen::DMatX3> spzAll;
      std::vector<eigen::DMat24> szAll;
      if (irank == mpi::rank()) {
        // gather coords
        spzAll.reserve(quads.size());
        szAll.reserve(quads.size());
        for (Quad& quad : quads) {
          spzAll.push_back(computeElemSPZ(quad));
          szAll.push_back(quad.getNodalSZ());
        }
        // send coords to super
        mpi::sendVecEigen(0, spzAll, 0);
        mpi::sendVecEigen(0, szAll, 1);
      }

      // step 2: compute values on super and send back to infer
      std::vector<eigen::DColX> undAll;
      std::vector<eigen::IColX> elemInScopeAll;
      if (mpi::root()) {
        // recv coords from infer
        mpi::recvVecEigen(irank, spzAll, 0);
        mpi::recvVecEigen(irank, szAll, 1);
        // allocate values
        int nQuad = (int)spzAll.size();
        undAll.reserve(nQuad);
        elemInScopeAll.push_back(eigen::IColX::Zero(nQuad));
        // compute values
        for (int iq = 0; iq < nQuad; iq++) {
          eigen::DColX und;
          elemInScopeAll[0](iq) = getUndulation(spzAll[iq], szAll[iq], und);
          undAll.push_back(und);
        }
        // send values to infer
        mpi::sendVecEigen(irank, undAll, 0);
        mpi::sendVecEigen(irank, elemInScopeAll, 1);
      }

      // step 3: set values to quads on infer
      if (irank == mpi::rank()) {
        // recv values from super
        mpi::recvVecEigen(0, undAll, 0);
        mpi::recvVecEigen(0, elemInScopeAll, 1);
        // set values to quads
        for (int iq = 0; iq < spzAll.size(); iq++) {
          if (elemInScopeAll[0](iq)) {
            setUndulationToQuad(undAll[iq], quads[iq]);
          }
        }
      }
      // do irank one by one
      mpi::barrier();
    }
    mpi::enterWorld();
  }
}

// set undulation to quad
void
Geometric3D::setUndulationToQuad(const eigen::DColX& undulation, Quad& quad) const {
  // flattened to structured
  const eigen::IRowN& pointNr = quad.getPointNr();
  eigen::arN_DColX undArr;
  int row = 0;
  for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
    int nr = pointNr(ipnt);
    undArr[ipnt] = undulation.block(row, 0, nr, 1);
    row += nr;
  }
  // set to Quad
  quad.getUndulationPtr()->addUndulation(undArr);
}

#include "StructuredGridG3D.hpp"
#include "Ellipticity.hpp"
#include "sg_tools.hpp"
#include "Crust1G3D.hpp"

// build from inparam
std::shared_ptr<const Geometric3D>
Geometric3D::buildInparam(const ExodusMesh& exodusMesh,
    const LocalMesh& localMesh,
    const std::string& modelName,
    const std::string& keyInparam) {
  // short alias
  const InparamYAML& gm = inparam::gInparamModel;
  const std::string& root = keyInparam;

  // class name
  const std::string& className = gm.get<std::string>(root + ":class_name");

  // clamp and Gaussian smoothing shared by geometric models
  const ClampSmoothConfig& clampSmooth = readClampSmooth(gm, root, modelName, className);

  // init class
  if (className == "StructuredGridG3D") {
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
    std::array<std::string, 2> crdVarNames;
    std::array<int, 2> shuffleData;
    sg_tools::inparamVarRank<2>(gm, rootc, modelName, className, crdVarNames, shuffleData);
    // units
    double lengthUnit = 1., angleUnit = 1.;
    sg_tools::inparamUnits(gm, rootc, xy, lengthUnit, angleUnit);

    ////////////// undulated range //////////////
    const std::string& rootr = root + ":undulation_range";
    double interface = gm.get<double>(rootr + ":interface");
    double min = gm.get<double>(rootr + ":min_max:[0]");
    double max = gm.get<double>(rootr + ":min_max:[1]");
    if (interface >= max || interface <= min) {
      throw std::runtime_error("Geometric3D::buildInparam || undulation_range:interface "
                               "must lie within undulation_range:min_max."
                               " || Model name: " +
          modelName + " || Class name: " + className);
    }

    ////////////// data //////////////
    const std::string& rootu = root + ":undulation_data";
    const std::string& dataVarName = gm.get<std::string>(rootu + ":nc_var");
    double factor = gm.get<double>(rootu + ":factor");
    bool superOnly = gm.get<bool>(root + ":store_grid_only_on_leaders");

    // construct
    return std::make_shared<const StructuredGridG3D>(modelName,
        fname,
        crdVarNames,
        shuffleData,
        sourceCentered,
        xy,
        ellipticity,
        useDepth,
        depthSolid,
        interface,
        min,
        max,
        lengthUnit,
        angleUnit,
        dataVarName,
        factor,
        superOnly,
        clampSmooth);
  } else if (className == "Ellipticity") {
    if (clampSmooth.enabled) {
      throw std::runtime_error("Geometric3D::buildInparam || "
                               "clamp_smooth is not supported for the analytical Ellipticity "
                               "model. || Model name: " +
          modelName);
    }
    // ellipticity can be added only once
    static bool ellipticityAdded = false;
    if (ellipticityAdded) {
      throw std::runtime_error("Geometric3D::buildInparam || Ellipticity "
                               "model cannot be added more than once.");
    }
    if (io::gVerboseWarnings) {
      if (geodesy::isCartesian()) {
        io::cout << bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                     "model will be ignored for a Cartesian mesh.");
      }
      if (geodesy::getOuterFlattening() < numerical::dEpsilon) {
        io::cout << bstring::warning("Geometric3D::buildInparam || Ellipticity "
                                     "model will be ignored with zero flattening.");
      }
    }
    ellipticityAdded = true;
    return std::make_shared<const Ellipticity>(modelName);
  } else if (className == "Crust1G3D") {
    double rSurf = gm.get<double>(root + ":surface_radius");
    double rMoho = gm.get<double>(root + ":moho_radius");
    double rBase = gm.get<double>(root + ":base_radius");
    bool includeSediment = gm.get<bool>(root + ":include_sediment");
    bool includeIce = gm.get<bool>(root + ":include_ice");
    bool ellipticity = gm.get<bool>(root + ":ellipticity");
    double surfaceFactor = gm.get<double>(root + ":surface_factor");
    double mohoFactor = gm.get<double>(root + ":moho_factor");
    return std::make_shared<const Crust1G3D>(modelName,
        rSurf,
        rMoho,
        rBase,
        includeSediment,
        includeIce,
        ellipticity,
        surfaceFactor,
        mohoFactor,
        clampSmooth);
  }

  // unknown class
  return nullptr;
}
