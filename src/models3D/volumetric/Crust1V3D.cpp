//
//  Crust1V3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 10/14/24.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

//  3D volumetric model based on Crust 1.0

#include "Crust1V3D.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include "ExodusMesh.hpp"
#include "spectrals.hpp"
#include "bstring.hpp"
#include "c1_tools.hpp"
#include <cmath>

// constructor
Crust1V3D::Crust1V3D(const std::string& modelName,
    double rSurf,
    double rMoho,
    bool includeSediment,
    bool includeIce,
    bool ellipticity,
    const ExodusMesh& exodusMesh) :
    Volumetric3D(modelName), mRSurf(rSurf), mRMoho(rMoho), mIncludeSediment(includeSediment),
    mIncludeIce(includeIce), mEllipticity(ellipticity) {
  // read file
  int nrow = sNLat * sNLon;
  eigen::DMatXX bnd, v_p, v_s, rho;
  bnd = v_p = v_s = rho = eigen::DMatXX::Zero(nrow, sNLayer);
  if (mpi::root()) {
    std::string path = io::gProjectDirectory + "/src/models3D/crust1_data";
    std::fstream fsbnd, fsv_p, fsv_s, fsrho;
    fsbnd.open(path + "/crust1.bnds", std::fstream::in);
    fsv_p.open(path + "/crust1.vp", std::fstream::in);
    fsv_s.open(path + "/crust1.vs", std::fstream::in);
    fsrho.open(path + "/crust1.rho", std::fstream::in);
    if (!fsbnd || !fsv_p || !fsv_s || !fsrho) {
      throw std::runtime_error("Crust1V3D::Crust1V3D || "
                               "Error opening crust1.0 data files "
                               "at directory: ||" +
          path);
    }
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < sNLayer; j++) {
        fsbnd >> bnd(i, j);
        fsv_p >> v_p(i, j);
        fsv_s >> v_s(i, j);
        fsrho >> rho(i, j);
      }
    }
    fsbnd.close();
    fsv_p.close();
    fsv_s.close();
    fsrho.close();
  }
  // broadcast
  mpi::bcastEigen(bnd);
  mpi::bcastEigen(v_p);
  mpi::bcastEigen(v_s);
  mpi::bcastEigen(rho);

  // cast to integer theta
  mRl = mVp = mVs = mRh = eigen::DMatXX::Zero(nrow + sNLon, sNLayer);
  for (int col = 0; col < sNLayer; col++) {
    mRl.block(0, col, sNLon, 1).fill(bnd.block(0, col, sNLon, 1).sum() / sNLon);
    mVp.block(0, col, sNLon, 1).fill(v_p.block(0, col, sNLon, 1).sum() / sNLon);
    mVs.block(0, col, sNLon, 1).fill(v_s.block(0, col, sNLon, 1).sum() / sNLon);
    mRh.block(0, col, sNLon, 1).fill(rho.block(0, col, sNLon, 1).sum() / sNLon);
    mRl.block(nrow, col, sNLon, 1).fill(bnd.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
    mVp.block(nrow, col, sNLon, 1).fill(v_p.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
    mVs.block(nrow, col, sNLon, 1).fill(v_s.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
    mRh.block(nrow, col, sNLon, 1).fill(rho.block(nrow - sNLon, col, sNLon, 1).sum() / sNLon);
  }
  for (int i = 1; i < sNLat; i++) {
    mRl.block(i * sNLon, 0, sNLon, sNLayer) =
        (bnd.block(i * sNLon, 0, sNLon, sNLayer) + bnd.block((i - 1) * sNLon, 0, sNLon, sNLayer)) *
        .5;
    mVp.block(i * sNLon, 0, sNLon, sNLayer) =
        (v_p.block(i * sNLon, 0, sNLon, sNLayer) + v_p.block((i - 1) * sNLon, 0, sNLon, sNLayer)) *
        .5;
    mVs.block(i * sNLon, 0, sNLon, sNLayer) =
        (v_s.block(i * sNLon, 0, sNLon, sNLayer) + v_s.block((i - 1) * sNLon, 0, sNLon, sNLayer)) *
        .5;
    mRh.block(i * sNLon, 0, sNLon, sNLayer) =
        (rho.block(i * sNLon, 0, sNLon, sNLayer) + rho.block((i - 1) * sNLon, 0, sNLon, sNLayer)) *
        .5;
  }
  // reverse south to north
  mRl = mRl.colwise().reverse().eval();
  mVp = mVp.colwise().reverse().eval();
  mVs = mVs.colwise().reverse().eval();
  mRh = mRh.colwise().reverse().eval();
  for (int i = 0; i <= sNLat; i++) {
    mRl.block(i * sNLon, 0, sNLon, sNLayer) =
        mRl.block(i * sNLon, 0, sNLon, sNLayer).colwise().reverse().eval();
    mVp.block(i * sNLon, 0, sNLon, sNLayer) =
        mVp.block(i * sNLon, 0, sNLon, sNLayer).colwise().reverse().eval();
    mVs.block(i * sNLon, 0, sNLon, sNLayer) =
        mVs.block(i * sNLon, 0, sNLon, sNLayer).colwise().reverse().eval();
    mRh.block(i * sNLon, 0, sNLon, sNLayer) =
        mRh.block(i * sNLon, 0, sNLon, sNLayer).colwise().reverse().eval();
  }

  // convert to SI
  mVp *= 1e3;
  mVs *= 1e3;
  mRh *= 1e3;
  mRl *= 1e3;

  // layers
  if (mIncludeIce) {
    mIncludeSediment = true;
  }
  int colSurf = c1_tools::columnSurf(mIncludeSediment, mIncludeIce);
  int colMoho = 8;

  // linear mapping to sphere
  const eigen::DColX& rmoho = eigen::DColX::Constant(nrow + sNLon, mRMoho);
  const eigen::DColX& rdiff = (eigen::DColX::Constant(nrow + sNLon, mRSurf - mRMoho).array() /
      (mRl.col(colSurf) - mRl.col(colMoho)).array())
                                  .matrix();
  const eigen::DMatXX copyRl = mRl;
  for (int i = 0; i < sNLayer; i++) {
    mRl.col(i).array() =
        rdiff.array() * (copyRl.col(i) - copyRl.col(colMoho)).array() + rmoho.array();
  }

  // above: physical layers
  ///////////////////////////////////////////////////////
  // blow: elemental layers

  // find all z-coordinates on the axis
  if (mpi::root()) {
    const auto& sz = exodusMesh.getNodalCoords();
    for (int i = 0; i < sz.rows(); i++) {
      double s = sz(i, 0);
      double z = sz(i, 1);
      if (s <= exodusMesh.getGlobalVariable("dist_tolerance") && z > 0)
        mElementBoundaries.push_back(z);
    }
    std::sort(mElementBoundaries.begin(), mElementBoundaries.end(), std::greater<double>());
  }
  mpi::bcast(mElementBoundaries);
  for (int i = 0; i < mElementBoundaries.size(); i++) {
    if (mElementBoundaries[i] < mRMoho - 1.) {
      mElementBoundaries.erase(mElementBoundaries.begin() + i, mElementBoundaries.end());
      break;
    }
  }
  for (int i = (int)mElementBoundaries.size() - 1; i >= 0; i--) {
    if (mElementBoundaries[i] > mRSurf + 1.) {
      mElementBoundaries.erase(mElementBoundaries.begin(), mElementBoundaries.begin() + i + 1);
      break;
    }
  }
  if (mElementBoundaries.size() < 2) {
    throw std::runtime_error("Crust1V3D::Crust1V3D || No element layer found in crust.");
  }
  if (std::abs(mRSurf - mElementBoundaries[0]) > 1.) {
    throw std::runtime_error("Crust1V3D::Crust1V3D || Conflict in surface radius.");
  }
  if (std::abs(mRMoho - mElementBoundaries[mElementBoundaries.size() - 1]) > 1.) {
    throw std::runtime_error("Crust1V3D::Crust1V3D || Conflict in Moho radius.");
  }

  // form gll boundaries
  int nEleCrust = (int)mElementBoundaries.size() - 1;
  const auto& eta = spectrals::gPositionsGLL;
  int numGll = 1 + spectral::nPol * nEleCrust;
  mRlGLL = eigen::DColX(numGll);
  mRlGLL(0) = mRSurf; // surface
  int index = 1;
  for (int iele = 0; iele < nEleCrust; iele++) {
    double eBoundTop = mElementBoundaries[iele];
    double eBoundBot = mElementBoundaries[iele + 1];
    double eHeight = eBoundTop - eBoundBot;
    for (int jpol = 1; jpol <= spectral::nPol; jpol++) {
      double z = eBoundTop - (eta(jpol) - eta(0)) / (eta(spectral::nPol) - eta(0)) * eHeight;
      mRlGLL(index++) = z;
    }
  }

  // zero properties
  mVpGLL = mVsGLL = mRhGLL = eigen::DMatXX::Zero(nrow + sNLon, numGll);

  // form values at GLL boundaries
  for (int igll = 0; igll < numGll; igll++) {
    int itop = igll >= 1 ? igll - 1 : 0;
    int imid = igll;
    int ibot = igll + 1;
    if (ibot > numGll - 1) {
      ibot = numGll - 1;
    }
    double gll_top = mRlGLL(itop);
    double gll_mid = mRlGLL(imid);
    double gll_bot = mRlGLL(ibot);
    double gll_mid_value = 2. / (gll_top - gll_bot);
    for (int row = 0; row < mRl.rows(); row++) {
      // integrate
      for (int ilayer = colSurf; ilayer < colMoho; ilayer++) {
        double phy_top = mRl(row, ilayer);
        double phy_bot = mRl(row, ilayer + 1);
        // top to mid
        double top = std::min(phy_top, gll_top);
        double bot = std::max(phy_bot, gll_mid);
        if (top > bot) {
          double gllt = gll_mid_value / (gll_top - gll_mid) * (gll_top - top);
          double gllb = gll_mid_value / (gll_top - gll_mid) * (gll_top - bot);
          double area = .5 * (gllt + gllb) * (top - bot);
          mVpGLL(row, igll) += mVp(row, ilayer) * area;
          mVsGLL(row, igll) += mVs(row, ilayer) * area;
          mRhGLL(row, igll) += mRh(row, ilayer) * area;
        }
        // mid to bot
        top = std::min(phy_top, gll_mid);
        bot = std::max(phy_bot, gll_bot);
        if (top > bot) {
          double gllt = gll_mid_value / (gll_mid - gll_bot) * (top - gll_bot);
          double gllb = gll_mid_value / (gll_mid - gll_bot) * (bot - gll_bot);
          double area = .5 * (gllt + gllb) * (top - bot);
          mVpGLL(row, igll) += mVp(row, ilayer) * area;
          mVsGLL(row, igll) += mVs(row, ilayer) * area;
          mRhGLL(row, igll) += mRh(row, ilayer) * area;
        }
      }
    }
  }

  // grid lat and lon
  mGridLat = eigen::DColX(sNLat + 1);
  mGridLon = eigen::DColX(sNLon + 1); // one bigger than data
  for (int i = 0; i < sNLat + 1; i++) {
    mGridLat[i] = i * 1. - 90.;
  }
  for (int i = 0; i < sNLon + 1; i++) {
    mGridLon[i] = i * 1. - 179.5;
  }
}

// get property info
void
Crust1V3D::getPropertyInfo(
    std::vector<std::string>& propKeys, std::vector<ReferenceKind>& refKinds) const {
  propKeys = {"VP", "VS", "RHO"};
  refKinds = {ReferenceKind::ABS, ReferenceKind::ABS, ReferenceKind::ABS};
}

// get properties
bool
Crust1V3D::getProperties(const eigen::DMatX3& spz,
    const eigen::DMat24& nodalSZ,
    eigen::IMatXX& inScopes,
    eigen::DMatXX& propValues) const {
  //////////////////////// coords ////////////////////////
  // check center
  if (!inplaneScope<eigen::DCol2>(
          nodalSZ.rowwise().mean(), false, 0., 0., true, mRMoho, mRSurf, false, false)) {
    return false;
  }

  // compute grid coords
  eigen::DMatX3 crdGrid =
      coordsFromMeshToModel(spz, false, false, mEllipticity, false, false, false, mModelName);

  //////////////////////// values ////////////////////////
  // allocate
  int nProperties = 3;
  int nCardinals = (int)spz.rows();
  inScopes = eigen::IMatXX::Ones(nCardinals, nProperties);
  propValues = eigen::DMatXX::Zero(nCardinals, nProperties);

  // point loop
  for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
    double lat = crdGrid(ipnt, 0);
    double lon = crdGrid(ipnt, 1);
    double r = crdGrid(ipnt, 2);
    const auto& val = getPropertiesPoint(r, lat, lon);
    propValues(ipnt, 0) = val[0];
    propValues(ipnt, 1) = val[1];
    propValues(ipnt, 2) = val[2];
  }
  return true;
}

// verbose
std::string
Crust1V3D::verbose() const {
  using namespace bstring;
  std::stringstream ss;
  ss << boxSubTitle(0, mModelName + " ", '~');
  ss << boxEquals(2, 22, "class name", "Crust1V3D");
  ss << boxEquals(2, 22, "surface radius", mRSurf);
  ss << boxEquals(2, 22, "moho radius", mRMoho);
  ss << boxEquals(2, 22, "include sediment", mIncludeSediment);
  ss << boxEquals(2, 22, "include ice", mIncludeIce);
  ss << boxEquals(2, 22, "ellipticity correction", mEllipticity);
  return ss.str();
}

// get properties at a point
std::vector<double>
Crust1V3D::getPropertiesPoint(double r, double lat, double lon) const {
  // regularise
  if (lon > 180.) {
    lon -= 360.;
  }
  if (lon < -179.5) {
    lon += 360.;
  }

  // interpolation on sphere
  int llat[2], llon[2];
  double wlat[2], wlon[2];
  c1_tools::interpLinear(lat, mGridLat, llat[0], wlat[0]);
  c1_tools::interpLinear(lon, mGridLon, llon[0], wlon[0]);
  llat[1] = llat[0] + 1;
  llon[1] = llon[0] + 1;
  wlat[1] = 1. - wlat[0];
  wlon[1] = 1. - wlon[0];
  if (llon[1] == sNLon) {
    llon[1] = 0;
  }

  // element inside but point slightly outside
  r = std::max(std::min(r, mRSurf * 0.999999), mRMoho * 1.000001);
  double v_p = 0.;
  double v_s = 0.;
  double rho = 0.;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      double weight = wlat[i] * wlon[j];
      int rowdata = llat[i] * sNLon + llon[j];
      bool found = false;
      for (int iLayer = 1; iLayer < mRlGLL.rows(); iLayer++) {
        if (r > mRlGLL(iLayer)) {
          double v_ptop = mVpGLL(rowdata, iLayer - 1);
          double v_stop = mVsGLL(rowdata, iLayer - 1);
          double rhotop = mRhGLL(rowdata, iLayer - 1);
          double v_pbot = mVpGLL(rowdata, iLayer);
          double v_sbot = mVsGLL(rowdata, iLayer);
          double rhobot = mRhGLL(rowdata, iLayer);
          double top = mRlGLL(iLayer - 1);
          double bot = mRlGLL(iLayer);
          double vp = (v_ptop - v_pbot) / (top - bot) * (r - bot) + v_pbot;
          double vs = (v_stop - v_sbot) / (top - bot) * (r - bot) + v_sbot;
          double rh = (rhotop - rhobot) / (top - bot) * (r - bot) + rhobot;
          v_p += vp * weight;
          v_s += vs * weight;
          rho += rh * weight;
          found = true;
          break;
        }
      }
      if (!found) {
        throw std::runtime_error("Crust1V3D::getPropertiesPoint || Please report this bug.");
      }
    }
  }

  // when sediment is considered, Vs at some locations may be too small
  v_s = std::max(v_s, 500.);
  v_p = std::max(v_p, sqrt(2) * v_s);
  return std::vector<double>{v_p, v_s, rho};
}
