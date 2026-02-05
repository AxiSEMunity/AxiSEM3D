//
//  SE_Model.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  spectral element model
//  Exodus mesh ---
//                | -> mesh --------
//  local mesh ----                | -> SE model (this) -> domain
//                     3D models ---


#ifndef SE_Model_hpp
#define SE_Model_hpp

#include "eigen_sem.hpp"

// constructor
class ExodusMesh;
class ABC;
class LocalMesh;
class Model3D;

// components
#include "Quad.hpp"
#include "GLLPoint.hpp"

// release
class AttBuilder;
class Domain;

// source-receiver
class Geometric3D;
#include "RTreeND.hpp"

class SE_Model {
public:
    // step 1: constructor
    SE_Model(const ExodusMesh &exodusMesh,
             const ABC &abc, const LocalMesh &localMesh,
             const std::vector<std::shared_ptr<const Model3D>> &models3D,
             bool useLuckyNumbers);
    
    // step 2: get dt for attenuation
    double computeDt(double courant, const ABC &abc, eigen::DCol2 &sz) const;
    
    // step 3: release to domain
    void release(const ABC &abc, const LocalMesh &localMesh,
                 const std::unique_ptr<const AttBuilder> &attBuilder,
                 const TimeScheme &timeScheme, Domain &domain);
    
    // step 4: measure
    eigen::DColX
    measureCost(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                const TimeScheme &timeScheme, bool measurePoint) const;
    
    // verbose
    std::string verbose(const std::string &titleSEM) const;
    
    
    ////////////////// source-receiver //////////////////
    // get total undulation
    double getTotalUndulation(const eigen::DRow3 &spzRef) const;
    
    // undulated to reference
    eigen::DRow3 undulatedToReference(const eigen::DRow3 &spzUnd) const;
    
    // form inplane RTree
    void formInplaneRTree();
    
    // locate inplane
    int locateInplane(const eigen::DCol2 &sz, bool inFluid) const;
    
    // compute inplane factor
    eigen::DRowN computeInplaneFactor(const eigen::DCol2 &sz, int iquad) const;
    
    // get quads
    const std::vector<Quad> &getQuads() const {
        return mQuads;
    }
    
    
    ////////////////// wavefield scanning //////////////////
    // initialize scanning on points
    void initScanningOnPoints(bool vertexOnly) const;
    
private:
    template <class Vec>
    static Vec interpLagrange(double target, const Vec &bases) {
        Vec weights = Vec::Ones(bases.rows(), bases.cols());
        for (int idgr = 0; idgr < bases.size(); idgr++) {
            for (int ibase = 0; ibase < bases.size(); ibase++) {
                if (ibase != idgr) {
                    weights[idgr] *= target - bases[ibase];
                    weights[idgr] /= bases[idgr] - bases[ibase];
                }
            }
        }
        return weights;
    }
    
    ////////////////////// data //////////////////////
private:
    ///////// components /////////
    // quads
    std::vector<Quad> mQuads;
    // GLL points
    std::vector<GLLPoint> mGLLPoints;
    
    ///////// source-receiver /////////
    // distance tolerance of mesh
    const double mDistToleranceMesh;
    
    // 3D geometric models
    std::vector<std::shared_ptr<const Geometric3D>> mModelsG3D;
    
    // local mesh range
    eigen::DCol2 mRangeS = eigen::DCol2::Zero();
    eigen::DCol2 mRangeZ = eigen::DCol2::Zero();
    eigen::DCol2 mRangeR = eigen::DCol2::Zero();
    eigen::DCol2 mRangeT = eigen::DCol2::Zero();
    
    // inplane
    std::unique_ptr<RTreeND<2, 1, int>> mRTreeFluid = nullptr;
    std::unique_ptr<RTreeND<2, 1, int>> mRTreeSolid = nullptr;
};

#endif /* SE_Model_hpp */
