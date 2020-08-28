//
//  Source.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/17/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  source: generator of ElementSource in core

#ifndef Source_hpp
#define Source_hpp

class SE_Model;
class Domain;
class STF;
class Mechanism;
#include <memory>

#include "eigen.hpp"
namespace eigen {
    typedef Eigen::Matrix<double, 1, 3> DRow3;
    typedef Eigen::Matrix<double, 3, 3> DMat33;
}

class Source {
public:
    // constructor
    Source(const std::string &name, bool axial,
           bool sourceCentered, bool ellipticity,
           bool useDepth, bool depthSolid, bool undulatedGeometry,
           const eigen::DRow3 &crdIn):
    mName(name), mAxial(axial),
    mSourceCentered(sourceCentered), mEllipticity(ellipticity),
    mUseDepth(useDepth), mDepthSolid(depthSolid),
    mUndulatedGeometry(undulatedGeometry), mCrdIn(crdIn) {
        // nothing
    }
    
private:
    // build from inparam
    static std::shared_ptr<const Source> buildInparam(int sindex);
    
public:
    // compute spz
    static eigen::DRow3
    computeSPZ(const SE_Model &sem, const eigen::DRow3 &crdIn,
               bool sourceCentered, bool xy, bool ellipticity,
               bool useDepth, bool depthSolid, bool undulatedGeometry,
               const std::string &errInfo, bool enforceOnAxis);
    
private:
    // compute rotation matrix Q from input to (z, s, phi)
    const eigen::DMat33 &computeQzsp(const eigen::DRow3 &spz,
                                     bool ellipticity) const;
    
    // verbose
    std::string verbose(int sindex, const STF &stf,
                        const Mechanism &mechanism) const;
    
public:
    // release sources to domain
    static void release(const SE_Model &sem, Domain &domain, double dt,
                        double &minT0);
    
private:
    // data from inparam
    const std::string mName;
    const bool mAxial;
    const bool mSourceCentered;
    const bool mEllipticity;
    const bool mUseDepth;
    const bool mDepthSolid;
    const bool mUndulatedGeometry;
    const eigen::DRow3 mCrdIn;
};

#endif /* Source_hpp */
