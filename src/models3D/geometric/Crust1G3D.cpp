//
//  Crust1G3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 16/10/2024.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

//  3D geometric model based on Crust 1.0

#include "Crust1G3D.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "bstring.hpp"
#include "c1_tools.hpp"
#include <cmath>

// constructor
Crust1G3D::
Crust1G3D(const std::string &modelName,
          double rSurf, double rMoho, double rBase,
          bool includeSediment, bool includeIce, bool ellipticity,
          double surfaceFactor, double mohoFactor,
          int gaussianOrder, double gaussianDev):
Geometric3D(modelName),
mRSurf(rSurf), mRMoho(rMoho), mRBase(rBase),
mIncludeSediment(includeSediment), mIncludeIce(includeIce),
mEllipticity(ellipticity),
mSurfFactor(surfaceFactor), mMohoFactor(mohoFactor),
mGaussianOrder(gaussianOrder), mGaussianDev(gaussianDev) {
    // read raw data
    int nrow = sNLat * sNLon;
    eigen::DMatXX elevation = eigen::DMatXX::Zero(nrow, sNLayer);
    if (mpi::root()) {
        std::string path =
        io::gProjectDirectory + "/src/models3D/crust1_data";
        std::fstream fs(path + "/crust1.bnds", std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Crust1G3D::Crust1G3D || "
                                     "Error opening crust1.0 data file: ||" +
                                     path + "/crust1.bnds");
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < sNLayer; j++) {
                fs >> elevation(i, j);
            }
        }
        fs.close();
    }
    // broadcast
    mpi::bcastEigen(elevation);
    
    // surface and moho undulations
    // NOTE: ellipticity should not be considered here because all geometric
    //       models are defined independently w.r.t. reference sphere
    // layers
    if (mIncludeIce) {
        mIncludeSediment = true;
    }
    int colSurf = c1_tools::columnSurf(mIncludeSediment, mIncludeIce);
    int colMoho = 8;
    eigen::DColX deltaRSurfVec = elevation.col(colSurf) * 1e3;
    eigen::DColX deltaRMohoVec = elevation.col(colMoho) * 1e3;
    deltaRMohoVec += eigen::DColX::Constant(nrow, mRSurf - mRMoho);
    // cast to matrix
    eigen::DMatXX deltaRSurf(sNLat, sNLon);
    eigen::DMatXX deltaRMoho(sNLat, sNLon);
    for (int i = 0; i < sNLat; i++) {
        deltaRSurf.row(i) = deltaRSurfVec.block(i * sNLon, 0, sNLon, 1).transpose();
        deltaRMoho.row(i) = deltaRMohoVec.block(i * sNLon, 0, sNLon, 1).transpose();
    }
    
    //////////// plot raw data ////////////
    // std::fstream fs;
    // fs.open("/Users/kuangdai/Desktop/crust1/raw.txt", std::fstream::out);
    // fs << deltaRMoho << std::endl;
    // fs.close();
    //////////// plot raw data ////////////
    
    // Gaussian smoothing
    eigen::IColX orderRow = eigen::IColX::Constant(sNLat, mGaussianOrder);
    eigen::IColX orderCol = eigen::IColX::Constant(sNLon, mGaussianOrder);
    eigen::DColX devRow = eigen::DColX::Constant(sNLat, mGaussianDev);
    eigen::DColX devCol = eigen::DColX::Constant(sNLon, mGaussianDev);
    // smooth poles more
    // int npolar = mGaussianOrder + 1;
    // int opolar = sNLon;
    // for (int i = 0; i < npolar; i++) {
    //     double order = (double)(opolar - mGaussianOrder) / npolar * (npolar - i) + mGaussianOrder;
    //     orderRow(i) = orderRow(sNLat - 1 - i) = round(order);
    // }
    c1_tools::gaussianSmoothing(deltaRSurf, orderRow, devRow, true, orderCol, devCol, false);
    c1_tools::gaussianSmoothing(deltaRMoho, orderRow, devRow, true, orderCol, devCol, false);
    
    // cast to integer theta with unique polar values
    mDeltaRSurf = eigen::DMatXX::Zero(sNLat + 1, sNLon);
    mDeltaRMoho = eigen::DMatXX::Zero(sNLat + 1, sNLon);
    // fill north and south pole
    mDeltaRSurf.row(0).fill(deltaRSurf.row(0).sum() / sNLon);
    mDeltaRMoho.row(0).fill(deltaRMoho.row(0).sum() / sNLon);
    mDeltaRSurf.row(sNLat).fill(deltaRSurf.row(sNLat - 1).sum() / sNLon);
    mDeltaRMoho.row(sNLat).fill(deltaRMoho.row(sNLat - 1).sum() / sNLon);
    // interp at integer theta
    for (int i = 1; i < sNLat; i++) {
        mDeltaRSurf.row(i) = (deltaRSurf.row(i - 1) + deltaRSurf.row(i)) * .5;
        mDeltaRMoho.row(i) = (deltaRMoho.row(i - 1) + deltaRMoho.row(i)) * .5;
    }
    // reverse south to north
    mDeltaRSurf = mDeltaRSurf.colwise().reverse().eval();
    mDeltaRMoho = mDeltaRMoho.colwise().reverse().eval();
    
    // apply factor
    mDeltaRSurf *= mSurfFactor;
    mDeltaRMoho *= mMohoFactor;
    
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

// get undulation on an element
bool Crust1G3D::getUndulation(const eigen::DMatX3 &spz,
                              const eigen::DMat24 &nodalSZ,
                              eigen::DColX &undulation) const {
    // check center
    if (!inplaneScope<eigen::DCol2>
        (nodalSZ.rowwise().mean(),
         false, 0., 0.,
         true, mRBase, mRSurf,
         false, false)) {
        return false;
    }
    eigen::DCol2 centerSZ = nodalSZ.rowwise().mean();
    double rElemCenter = centerSZ(1);
    if (!geodesy::isCartesian()) {
        rElemCenter = sqrt(centerSZ(0) * centerSZ(0) +
                           centerSZ(1) * centerSZ(1));
    }
    
    eigen::DMatX3 crdGrid =
    coordsFromMeshToModel(spz, false, false, mEllipticity, false,
                          false, false, mModelName);
    int nCardinals = (int)spz.rows();
    undulation = eigen::DColX::Zero(nCardinals);
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        double lat = crdGrid(ipnt, 0);
        double lon = crdGrid(ipnt, 1);
        double r = crdGrid(ipnt, 2);
        undulation(ipnt) = getUndulationPoint(r, lat, lon, rElemCenter);
    }
    return true;
}

// get undulation on points
bool Crust1G3D::getUndulation(const eigen::DMatX3 &spz,
                              eigen::DColX &undulation) const {
    // compute grid coords
    eigen::DMatX3 crdGrid =
    coordsFromMeshToModel(spz, false, false, mEllipticity, false,
                          false, false, mModelName);
    int nCardinals = (int)spz.rows();
    undulation = eigen::DColX::Zero(nCardinals);
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        double lat = crdGrid(ipnt, 0);
        double lon = crdGrid(ipnt, 1);
        double r = crdGrid(ipnt, 2);
        undulation(ipnt) = getUndulationPoint(r, lat, lon, r);
    }
    return true;
}

// verbose
std::string Crust1G3D::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxSubTitle(0, mModelName + " ", '~');
    ss << boxEquals(2, 22, "class name", "Crust1G3D");
    ss << boxEquals(2, 22, "surface radius", mRSurf);
    ss << boxEquals(2, 22, "moho radius", mRMoho);
    ss << boxEquals(2, 22, "base radius", mRBase);
    ss << boxEquals(2, 22, "include sediment", mIncludeSediment);
    ss << boxEquals(2, 22, "include ice", mIncludeIce);
    ss << boxEquals(2, 22, "ellipticity correction", mEllipticity);
    ss << boxEquals(2, 22, "surface factor", mSurfFactor);
    ss << boxEquals(2, 22, "moho factor", mMohoFactor);
    ss << boxEquals(2, 22, "Gaussian order", mGaussianOrder);
    ss << boxEquals(2, 22, "Gaussian factor", mGaussianDev);
    return ss.str();
}

// get undulation on point
double Crust1G3D::getUndulationPoint(double r, double lat, double lon,
                                     double rElemCenter) const {
    if (rElemCenter < mRBase || rElemCenter > mRSurf) {
        return 0.;
    }
    
    // regularise
    if (lon > 180.) {
        lon -= 360.;
    }
    if (lon < -179.5) {
        lon += 360.;
    }
    
    // interpolation on sphere
    int llat0, llon0, llat1, llon1;
    double wlat0, wlon0, wlat1, wlon1;
    c1_tools::interpLinear(lat, mGridLat, llat0, wlat0);
    c1_tools::interpLinear(lon, mGridLon, llon0, wlon0);
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    if (llon1 == sNLon) {
        llon1 = 0;
    }
    
    double drSurf = 0.;
    drSurf += mDeltaRSurf(llat0, llon0) * wlat0 * wlon0;
    drSurf += mDeltaRSurf(llat1, llon0) * wlat1 * wlon0;
    drSurf += mDeltaRSurf(llat0, llon1) * wlat0 * wlon1;
    drSurf += mDeltaRSurf(llat1, llon1) * wlat1 * wlon1;
    
    double drMoho = 0.;
    drMoho += mDeltaRMoho(llat0, llon0) * wlat0 * wlon0;
    drMoho += mDeltaRMoho(llat1, llon0) * wlat1 * wlon0;
    drMoho += mDeltaRMoho(llat0, llon1) * wlat0 * wlon1;
    drMoho += mDeltaRMoho(llat1, llon1) * wlat1 * wlon1;
    
    // interpolation along radius
    if (rElemCenter < mRMoho) {
        return drMoho / (mRMoho - mRBase) * (r - mRBase);
    } else {
        return (drSurf - drMoho) / (mRSurf - mRMoho) * (r - mRMoho) + drMoho;
    }
}
