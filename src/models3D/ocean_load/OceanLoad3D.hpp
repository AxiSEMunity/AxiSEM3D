//
//  OceanLoad3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D ocean-load models

#ifndef OceanLoad3D_hpp
#define OceanLoad3D_hpp

#include "Model3D.hpp"

class OceanLoad3D: public Model3D {
public:
    // constructor
    OceanLoad3D(const std::string &modelName): Model3D(modelName) {
        // nothing
    }
    
    // destructor
    virtual ~OceanLoad3D() = default;
    
    // apply to Quad
    virtual void applyTo(std::vector<Quad> &quads) const;
    
protected:
    // get sum(rho * depth)
    virtual bool getSumRhoDepth(const eigen::DMatX3 &spz,
                                const eigen::DMat24 &nodalSZ,
                                eigen::DColX &sumRhoDepth) const = 0;
    
    // set sum(rho * depth) to quad
    virtual void setSumRhoDepthToQuad(const eigen::DColX &sumRhoDepth,
                                      Quad &quad) const;
    
    
    ////////////////////////////// static //////////////////////////////
public:
    // build from inparam
    static std::shared_ptr<const OceanLoad3D>
    buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                 const std::string &modelName, const std::string &keyInparam);
};

#endif /* OceanLoad3D_hpp */
