//
//  Volumetric3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/12/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  3D volumetric models
//  density, velocity, elasticity, attenuation

#ifndef Volumetric3D_hpp
#define Volumetric3D_hpp

#include "Model3D.hpp"
#include <map>

namespace eigen {
    // anisotropy
    typedef Eigen::Matrix<int, 6, 6> IMat66;
    typedef Eigen::Matrix<double, 6, 6> DMat66;
    typedef Eigen::Matrix<double, 3, 3> DMat33;
}

class Volumetric3D: public Model3D {
public:
    // reference kind
    enum class ReferenceKind {ABS, REF1D, REF3D, REF_PERTURB};
    inline static const
    std::map<ReferenceKind, std::string> sReferenceKindStr = {
        {ReferenceKind::ABS, "ABS"},
        {ReferenceKind::REF1D, "REF1D"},
        {ReferenceKind::REF3D, "REF3D"},
        {ReferenceKind::REF_PERTURB, "REF_PERTURB"}};
    
    // constructor
    Volumetric3D(const std::string &modelName): Model3D(modelName) {
        // nothing
    }
    
    // destructor
    virtual ~Volumetric3D() = default;
    
    // apply to Quad
    virtual void applyTo(std::vector<Quad> &quads) const;
    
protected:
    // using reference or undulated geometry
    virtual bool usingUndulatedGeometry() const = 0;
    
    // get property info
    virtual
    void getPropertyInfo(std::vector<std::string> &propKeys,
                         std::vector<ReferenceKind> &refKinds) const = 0;
    
    // get properties
    virtual bool getProperties(const eigen::DMatX3 &spz,
                               const eigen::DMat24 &nodalSZ,
                               eigen::IMatXX &inScopes,
                               eigen::DMatXX &propValues) const = 0;
    
    // set properties to quad
    virtual void setPropertiesToQuad(const std::vector<std::string> &propKeys,
                                     const std::vector<ReferenceKind> &refKinds,
                                     const eigen::IMatXX &inScopes,
                                     const eigen::DMatXX &propValues,
                                     Quad &quad) const;
    
    
    ////////////////////////////// static //////////////////////////////
    // Bond transformation for rotating Cijkl
    static void bondTransformation(const eigen::DMat66 &inCijkl,
                                   double alpha, double beta, double gamma,
                                   eigen::DMat66 &outCijkl);
    
public:
    // build from inparam
    static std::shared_ptr<const Volumetric3D>
    buildInparam(const ExodusMesh &exodusMesh, const LocalMesh &localMesh,
                 const std::string &modelName, const std::string &keyInparam);
};

#endif /* Volumetric3D_hpp */
