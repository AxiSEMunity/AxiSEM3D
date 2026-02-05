//
//  Element.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  spectral element

#ifndef Element_hpp
#define Element_hpp

// components
#include "CoordTransform.hpp"
#include "GradientQuadrature.hpp"
#include "PRT.hpp"

// Real-typed GradientQuadrature
typedef GradientQuadrature<numerical::Real> GradQuad;

// point
class Point;

class Element {
public:
    // constructor
    Element(int quadTag, std::unique_ptr<const GradQuad> &gq,
            std::unique_ptr<const PRT> &prt):
    mQuadTag(quadTag), mGradQuad(gq.release()), mPRT(prt.release()) {
        mTransform = nullptr;
    }
    
    // copy constructor
    Element(const Element &other):
    mQuadTag(other.mQuadTag),
    mGradQuad(std::make_unique<GradQuad>(*(other.mGradQuad))),
    mPRT(other.mPRT == nullptr ? nullptr :
         std::make_unique<PRT>(*(other.mPRT))) {
        mTransform = nullptr;
    }
    
    // destructor
    virtual ~Element() = default;
    
    // type info
    virtual std::string typeInfo() const = 0;
    
    // measure cost
    virtual double measure(int count) const = 0;
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << typeInfo() << "$" << mNr;
        if (mGradQuad->axial()) {
            ss << "$Axial";
        }
        return ss.str();
    }
    
    
    /////////////////////////// point ///////////////////////////
    // get point
    virtual Point &getPoint(int ipnt) const = 0;
    
    // point set
    void pointSet(bool elemInFourier);
    
    // get nr
    int getNr() const {
        return mNr;
    }
    
    // get nu
    int getNu_1() const {
        return mNu_1;
    }
    
    // get quad tag
    int getQuadTag() const {
        return mQuadTag;
    }
    
    // find boundary points by tag
    std::vector<int>
    findBoundaryPointsByTag(const std::vector<int> &boundaryMeshTags) const;
    
    // find boundary points by crds
    std::vector<int>
    findBoundaryPointsByCrds(const std::vector<double> &boundaryCrdsRorZ,
                             const std::vector<double> &boundaryCrdsTorS,
                             double distTol) const;
    
    
    /////////////////////////// domain ///////////////////////////
    // set domain tag
    void setDomainTag(int domainTag) {
        mDomainTag = domainTag;
    }
    
    // get domain tag
    int getDomainTag() const {
        return mDomainTag;
    }
    
    
protected:
    /////////////////////////// crd transform ///////////////////////////
    // set to RTZ by material or PRT
    void setToRTZ_ByMaterialOrPRT();
    
    // create coordinate transform
    // internally called by setToRTZ_ByMaterialOrPRT
    // externally called by prepareWavefieldOutput
    void createCoordTransform();
    
    
protected:
    // tag in the spectral-element mesh
    // this tag is mpi-independent and unique across processors
    const int mQuadTag;
    
    // gradient operator
    const std::unique_ptr<const GradQuad> mGradQuad;
    
    // particle relabelling
    const std::unique_ptr<const PRT> mPRT;
    
    // order
    int mNr = 0;
    int mNu_1 = 0;
    
    // tag (position) in the computational domain
    // this tag is mpi-dependent
    int mDomainTag = -1;
    
    // coordinate transform
    std::unique_ptr<const CoordTransform> mTransform;
};

#endif /* Element_hpp */
