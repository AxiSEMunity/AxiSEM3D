//
//  Crust1O3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 17/10/2024.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

#include "Crust1O3D.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "bstring.hpp"
#include "c1_tools.hpp"
#include <cmath>

// constructor
Crust1O3D::Crust1O3D(const std::string &modelName, double waterDensity,
                     bool includeIceAsWater, bool ellipticity,
                     int gaussianOrder, double gaussianDev):
OceanLoad3D(modelName),
mWaterDensity(waterDensity),
mIncludeIceAsWater(includeIceAsWater),
mEllipticity(ellipticity),
mGaussianOrder(gaussianOrder), mGaussianDev(gaussianDev) {
    // read raw data
    int nrow = sNLat * sNLon;
    eigen::DMatXX elevation = eigen::DMatXX::Zero(nrow, sNLayer);
    if (mpi::root()) {
        std::string path =
        io::gProjectDirectory + "/src/models3D/crust1_data";
        std::fstream fs(path + "/crust1.bnds", std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Crust1O3D::Crust1O3D || "
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
    
    // water depth
    int colWaterBot = 1;
    if (mIncludeIceAsWater) {
        colWaterBot = 2;
    }
    eigen::DColX depthVec = (elevation.col(0) - elevation.col(colWaterBot)) * 1e3;
    eigen::DMatXX depth(sNLat, sNLon);
    for (int i = 0; i < sNLat; i++) {
        depth.row(i) = depthVec.block(i * sNLon, 0, sNLon, 1).transpose();
    }
    
    //////////// plot raw data ////////////
    // std::fstream fs;
    // fs.open("/Users/kuangdai/Desktop/crust1/water.txt", std::fstream::out);
    // fs << depth << std::endl;
    // fs.close();
    //////////// plot raw data ////////////
    
    // Gaussian smoothing
    eigen::IColX orderRow = eigen::IColX::Constant(sNLat, mGaussianOrder);
    eigen::IColX orderCol = eigen::IColX::Constant(sNLon, mGaussianOrder);
    eigen::DColX devRow = eigen::DColX::Constant(sNLat, mGaussianDev);
    eigen::DColX devCol = eigen::DColX::Constant(sNLon, mGaussianDev);
    c1_tools::gaussianSmoothing(depth, orderRow, devRow, true, orderCol, devCol, false);
    
    // cast to integer theta with unique polar values
    mDepth = eigen::DMatXX::Zero(sNLat + 1, sNLon);
    // fill north and south pole
    mDepth.row(0).fill(depth.row(0).sum() / sNLon);
    mDepth.row(sNLat).fill(depth.row(sNLat - 1).sum() / sNLon);
    // interp at integer theta
    for (int i = 1; i < sNLat; i++) {
        mDepth.row(i) = (depth.row(i - 1) + depth.row(i)) * .5;
    }
    // reverse south to north
    mDepth = mDepth.colwise().reverse().eval();
    
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


// get sum(rho * depth)
bool Crust1O3D::getSumRhoDepth(const eigen::DMatX3 &spz,
                               const eigen::DMat24 &nodalSZ,
                               eigen::DColX &sumRhoDepth) const {
    // compute grid coords
    const eigen::DMatX3 &crdGrid =
    coordsFromMeshToModel(spz, false, false, mEllipticity, false,
                          false, false, mModelName);
    
    int nCardinals = (int)spz.rows();
    sumRhoDepth = eigen::DColX::Zero(nCardinals);
    for (int ipnt = 0; ipnt < nCardinals; ipnt++) {
        double lat = crdGrid(ipnt, 0);
        double lon = crdGrid(ipnt, 1);
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
        
        double depth = 0.;
        depth += mDepth(llat0, llon0) * wlat0 * wlon0;
        depth += mDepth(llat1, llon0) * wlat1 * wlon0;
        depth += mDepth(llat0, llon1) * wlat0 * wlon1;
        depth += mDepth(llat1, llon1) * wlat1 * wlon1;
        sumRhoDepth(ipnt) = depth * mWaterDensity;
    }
    return true;
}

// verbose
std::string Crust1O3D::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxSubTitle(0, mModelName + " ", '~');
    ss << boxEquals(2, 22, "class name", "Crust1O3D");
    ss << boxEquals(2, 22, "water density", mWaterDensity);
    ss << boxEquals(2, 22, "include ice as water", mIncludeIceAsWater);
    ss << boxEquals(2, 22, "ellipticity correction", mEllipticity);
    ss << boxEquals(2, 22, "Gaussian order", mGaussianOrder);
    ss << boxEquals(2, 22, "Gaussian factor", mGaussianDev);
    return ss.str();
}
