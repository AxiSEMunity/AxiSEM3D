//
//  Model3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/23/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of all 3D models

#ifndef Model3D_hpp
#define Model3D_hpp

#include "geodesy.hpp"
#include "eigen_sem.hpp"
#include <memory>

class ExodusMesh;
class LocalMesh;
class Quad;

class Model3D {
public:
    // constructor
    Model3D(const std::string &modelName): mModelName(modelName) {
        // nothing
    }
    
    // destructor
    virtual ~Model3D() = default;
    
    // apply to Quad
    virtual void applyTo(std::vector<Quad> &quads) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // get model name
    const std::string &getModelName() const {
        return mModelName;
    }
    
protected:
    // mpi-dependent: rebuild after domain re-partitioning
    virtual bool isMPI_Dependent() const {
        return false;
    }
    
    // super-only: data stored only on super ranks
    virtual bool isSuperOnly() const {
        return false;
    }
    
    // model name
    const std::string mModelName;
    
    
    ////////////////////////////// static //////////////////////////////
    // compute spz on element
    static eigen::DMatX3 computeElemSPZ(const Quad &quad,
                                        bool undulated = false);
    
    // compute spz on edge
    static eigen::DMatX3 computeEdgeSPZ(const Quad &quad, int edge);
    
    
    // check inplane scope
    template <class Mat2X>
    static bool inplaneScope(const Mat2X &szMesh,
                             bool checkR, double minR, double maxR,
                             bool checkZ, double minZ, double maxZ,
                             bool useDepth, bool depthSolid) {
        // check nothing
        if ((!checkR) && (!checkZ)) {
            return true;
        }
        
        // compute coords in model CS
        Mat2X crdModel;
        if (geodesy::isCartesian()) {
            crdModel = szMesh;
        } else {
            // (s, z) -> (theta, r)
            crdModel = geodesy::sz2rtheta(szMesh, false, 0, 1, 1, 0);
        }
        
        // check R
        if (checkR) {
            if (crdModel.row(0).maxCoeff() < minR ||
                crdModel.row(0).minCoeff() > maxR) {
                return false;
            }
        }
        
        // check Z
        if (checkZ) {
            // radius -> depth
            if (useDepth) {
                double router = (depthSolid ? geodesy::getOuterSolidRadius() :
                                 geodesy::getOuterRadius());
                crdModel.row(1) = router - crdModel.row(1).array();
            }
            if (crdModel.row(1).maxCoeff() < minZ ||
                crdModel.row(1).minCoeff() > maxZ) {
                return false;
            }
        }
        return true;
    }
    
    // compute coords in model CS
    static eigen::DMatX3
    coordsFromMeshToModel(const eigen::DMatX3 &spzMesh, bool sourceCentered,
                          bool xy, bool ellipticity, bool lon360,
                          bool useDepth, bool depthSolid,
                          const std::string &modelName) {
        eigen::DMatX3 crdModel;
        if (sourceCentered) {
            if (geodesy::isCartesian()) {
                crdModel = spzMesh;
                // (s, phi) to (x, y)
                if (xy) {
                    crdModel.col(0) =
                    spzMesh.col(0).array() * spzMesh.col(1).array().cos();
                    crdModel.col(1) =
                    spzMesh.col(0).array() * spzMesh.col(1).array().sin();
                }
            } else {
                // spz -> RTZ
                crdModel = geodesy::sz2rtheta(spzMesh, true, 0, 2, 2, 0);
                crdModel.col(1) = spzMesh.col(1);
            }
        } else {
            // correct spz for Cartesian
            // NOTE: Geographic and Cartesian contradict each other.
            //       We correct spz using s as arc-length and z as radius.
            //       Without such "bending", sqrt(s*s+z*z) on the surface
            //       will exceed the outer radius.
            if (geodesy::isCartesian()) {
                // z must be significantly larger than s
                if ((spzMesh.array().col(2) <
                     spzMesh.array().col(0) * 10.).any()) {
                    throw std::runtime_error
                    ("Model3D::coordsFromMeshToModel || "
                     "Invalid geographic location in Cartesian mesh. || "
                     "Set PLANET_RADIUS in Salvus mesher data file (*.bm). || "
                     "Model name = " + modelName);
                }
                eigen::DMatX3 spzBend(spzMesh.rows(), spzMesh.cols());
                const auto &t = spzMesh.col(0).cwiseQuotient(spzMesh.col(2));
                spzBend.col(0) = spzMesh.col(2).array() * t.array().sin();
                spzBend.col(2) = spzMesh.col(2).array() * t.array().cos();
                spzBend.col(1) = spzMesh.col(1);
                // (s, phi, z) to (lat, lon, r)
                crdModel = geodesy::spz2llr(spzBend, ellipticity, lon360);
            } else {
                // (s, phi, z) to (lat, lon, r)
                crdModel = geodesy::spz2llr(spzMesh, ellipticity, lon360);
            }
        }
        
        // depth (must be in reference geometry)
        if (useDepth) {
            double R = (depthSolid ? geodesy::getOuterSolidRadius() :
                        geodesy::getOuterRadius());
            crdModel.col(2).array() = R - crdModel.col(2).array();
        }
        return crdModel;
    }
    
public:
    // build from inparam
    static void
    buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                 std::vector<std::shared_ptr<const Model3D>> &models3D,
                 bool rebuilding);
};

#endif /* Model3D_hpp */
