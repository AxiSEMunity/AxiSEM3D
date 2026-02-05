//
//  Model3D.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of all 3D models

#include "Model3D.hpp"
#include "Quad.hpp"
#include "vicinity.hpp"

// compute spz on element
eigen::DMatX3 Model3D::computeElemSPZ(const Quad &quad, bool undulated) {
    // quad nr
    const eigen::IRowN &pointNr = quad.getPointNr();
    // GLL coordinates
    const eigen::DMat2N &pointSZ = quad.getPointSZ();
    // allocate
    eigen::DMatX3 spz(pointNr.sum(), 3);
    // structured to flattened
    int row = 0;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        int nr = pointNr(ipnt);
        // identical (s, z)
        spz.block(row, 0, nr, 1).fill(pointSZ(0, ipnt));
        spz.block(row, 2, nr, 1).fill(pointSZ(1, ipnt));
        // linearly varying phi
        spz.block(row, 1, nr, 1) =
        eigen::DColX::LinSpaced(nr, 0, 2. * numerical::dPi / nr * (nr - 1));
        row += nr;
    }
    
    // undulated geometry
    if (!undulated) {
        return spz;
    }
    // get undulation from quad
    eigen::arN_DColX und = quad.getUndulation();
    // structured to flattened
    row = 0;
    for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
        int nr = pointNr(ipnt);
        // extend to 3D
        if (und[ipnt].rows() == 1) {
            und[ipnt] = und[ipnt].replicate(nr, 1);
        }
        // change coords based on CS type
        if (geodesy::isCartesian()) {
            // add undulation to z in Cartesian
            spz.block(row, 2, nr, 1) += und[ipnt];
        } else {
            // add undulation to r in spherical
            double s = spz(row, 0);
            double z = spz(row, 2);
            double r = sqrt(s * s + z * z);
            double theta = (r < numerical::dEpsilon) ? 0. : acos(z / r);
            spz.block(row, 0, nr, 1) += und[ipnt] * sin(theta);
            spz.block(row, 2, nr, 1) += und[ipnt] * cos(theta);
        }
        row += nr;
    }
    return spz;
}

// compute spz on edge
eigen::DMatX3 Model3D::computeEdgeSPZ(const Quad &quad, int edge) {
    // edge points
    const std::vector<int> &ipnts = vicinity::constants::gEdgeIPnt[edge];
    // quad nr
    const eigen::IRowN &pointNr = quad.getPointNr();
    // GLL coordinates
    const eigen::DMat2N &pointSZ = quad.getPointSZ();
    // allocate
    eigen::DMatX3 spz(pointNr(ipnts).sum(), 3);
    // structured to flattened
    int row = 0;
    for (int ipnt: ipnts) {
        int nr = pointNr(ipnt);
        // identical (s, z)
        spz.block(row, 0, nr, 1).fill(pointSZ(0, ipnt));
        spz.block(row, 2, nr, 1).fill(pointSZ(1, ipnt));
        // linearly varying phi
        spz.block(row, 1, nr, 1) =
        eigen::DColX::LinSpaced(nr, 0, 2. * numerical::dPi / nr * (nr - 1));
        row += nr;
    }
    return spz;
}


#include "Volumetric3D.hpp"
#include "Geometric3D.hpp"
#include "OceanLoad3D.hpp"
#include "inparam.hpp"
#include "timer.hpp"

// build from inparam
void Model3D::
buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
             std::vector<std::shared_ptr<const Model3D>> &models3D,
             bool rebuilding) {
    // count
    int modelCount = inparam::gInparamModel.get<int>("list_of_3D_models:[?]");
    
    // loop over model list
    int mindexActive = -1;
    for (int mindex = 0; mindex < modelCount; mindex++) {
        // model name and key in inparam
        std::string keyInparam = "list_of_3D_models";
        const std::string &modelName = inparam::gInparamModel.
        get<std::string>(keyInparam + ":{" + bstring::toString(mindex) + "}");
        keyInparam += ":[" + bstring::toString(mindex) + "]:" + modelName;
        
        // start building
        timer::gPreloopTimer.begin("Building 3D model: " + modelName);
        
        // activated or not
        if (!inparam::gInparamModel.get<bool>(keyInparam + ":activated")) {
            timer::gPreloopTimer.message("Model is deactivated.");
            timer::gPreloopTimer.ended("Building 3D model: " + modelName);
            continue;
        }
        mindexActive++;
        
        // do not rebuild if model is mpi-independent
        if (rebuilding) {
            if (!models3D[mindexActive]->isMPI_Dependent()) {
                timer::gPreloopTimer.message("Model is MPI-independent, "
                                             "no need to rebuild.");
                timer::gPreloopTimer.ended("Building 3D model: " + modelName);
                continue;
            }
        }
        
        // try volumetric first
        std::shared_ptr<const Model3D>
        model = Volumetric3D::buildInparam(exodusMesh, localMesh,
                                           modelName, keyInparam);
        // try geometric
        if (model == nullptr) {
            model = Geometric3D::buildInparam(exodusMesh, localMesh,
                                              modelName, keyInparam);
        }
        
        // try ocean load
        if (model == nullptr) {
            model = OceanLoad3D::buildInparam(exodusMesh, localMesh,
                                              modelName, keyInparam);
        }
        
        // unknown model class
        if (model == nullptr) {
            throw std::runtime_error("Model3D::buildInparam || "
                                     "Unknown 3D model: " + modelName);
        }
        
        // replace or push back
        if (rebuilding) {
            models3D[mindexActive].reset(model.get());
        } else {
            models3D.push_back(model);
        }
        timer::gPreloopTimer.ended("Building 3D model: " + modelName);
    }
}
