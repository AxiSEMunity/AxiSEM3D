//
//  ElementInteriorWJ.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/20/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  element on the interior boundary of wavefield injection
//  template for both solid and fluid, passing class as base

#ifndef ElementInteriorWJ_hpp
#define ElementInteriorWJ_hpp

// element
#include "SolidElement.hpp"
#include "FluidElement.hpp"
#include "mapPPvsN.hpp"

// stf
#include "SourceTimeFunctionNetCDF.hpp"

// receiver
#include "eigen_generic.hpp"
#include "eigen_point.hpp"

namespace eigen {
    typedef Eigen::Matrix<numerical::ComplexR, 1, spectral::nPEM> CRowN;
}

template <class Element, int ndim>
class ElementInteriorWJ: public Element {
    // matrix
    typedef std::vector<std::array<eigen::CMatPP_RM, ndim>> vec_arD_CMatPP_RM;
    
    // stf
    typedef SourceTimeFunctionNetCDF<numerical::Real, ndim> STF_NetCDF;
    
public:
    // constructor
    ElementInteriorWJ(const Element &elementBase,
                      const std::vector<int> &pointsOnWJB):
    Element(elementBase), mPointsOnWJB(pointsOnWJB),
    mIncidentElement(std::make_unique<Element>(elementBase)) {
        // injected displacement
        int nu_1 = this->getNu_1();
        mInjectedDispl = std::make_unique<vec_arD_CMatPP_RM>();
        mInjectedDispl->resize(nu_1);
        
        // interior stiffness
        if (sInteriorStiff.size() < nu_1) {
            sInteriorStiff.resize(nu_1);
        }
        
        // domain tag (leaving that of mIncidentElement -1)
        this->setDomainTag(elementBase.getDomainTag());
    }
    
    
    /////////////////////////// as an element ///////////////////////////
private:
    // collect displacement from points
    void collectDisplFromPoints(vec_arD_CMatPP_RM &displElem) const {
        // base
        Element::collectDisplFromPoints(displElem);
        
        // injection boundary: add injected displacement
        // NOTE: this will only affect interior elements
        for (int alpha = 0; alpha < this->getNu_1(); alpha++) {
            for (int idim = 0; idim < ndim; idim++) {
                Eigen::Map<eigen::CRowN>
                (displElem[alpha][idim].data())(mPointsOnWJB) +=
                Eigen::Map<const eigen::CRowN>
                ((*mInjectedDispl)[alpha][idim].data())(mPointsOnWJB);
            }
        }
    }
    
    // add stiffness to points
    void addStiffToPoints(vec_arD_CMatPP_RM &stiffElem) const {
        // compute injected stiffness
        mIncidentElement->displToStiff(*mInjectedDispl, sInteriorStiff);
        
        // injection boundary: subtract injected stiffness
        // NOTE: this will also affect exterior elements
        for (int alpha = 0; alpha < this->getNu_1(); alpha++) {
            for (int idim = 0; idim < ndim; idim++) {
                Eigen::Map<eigen::CRowN>
                (stiffElem[alpha][idim].data())(mPointsOnWJB) -=
                Eigen::Map<const eigen::CRowN>
                (sInteriorStiff[alpha][idim].data())(mPointsOnWJB);
            }
        }
        
        // base
        Element::addStiffToPoints(stiffElem);
    }
    
    
    /////////////////////////// as a source ///////////////////////////
public:
    std::string quadKey() const {
        std::stringstream ss;
        ss << this->mediumInfo() << "_QUAD_" << this->getQuadTag();
        return ss.str();
    }
    
    // set stf
    void setSTF(bool aligned, int bufferSize,
                const std::shared_ptr<NetCDF_Reader> &reader) {
        // stf
        mSTF = std::make_unique<STF_NetCDF>(aligned, bufferSize,
                                            reader, quadKey());
        
        // workspace
        if (sInjectedDisplFlat.rows() < mSTF->getPatternNu_1()) {
            sInjectedDisplFlat.resize(mSTF->getPatternNu_1(),
                                      spectral::nPEM * ndim);
        }
    }
    
    // apply
    void applyInjection(int timeStep, double time) const {
        // get injected displacement
        mSTF->getPatternAtTimeStep(timeStep, time, sInjectedDisplFlat);
        
        // copy to element
        int nu_1 = std::min(mSTF->getPatternNu_1(), this->getNu_1());
        mapPPvsN::N2PP(sInjectedDisplFlat, *mInjectedDispl, nu_1);
    }
    
    
    /////////////////////////// as a receiver ///////////////////////////
    // receiver info
    int receiverInfo(std::vector<std::string> &keys,
                     eigen::DMatXX &coords_spz) const {
        // number of receivers
        int nr = this->getNr();
        int nrec = spectral::nPEM * nr;
        keys.resize(nrec);
        coords_spz.resize(nrec, 3);
        
        // s, z
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            const eigen::DRow2 &sz = this->getPoint(ipnt).getCoords();
            auto seq = Eigen::seqN(ipnt * nr, nr);
            coords_spz(seq, 0).fill(sz(0));
            coords_spz(seq, 2).fill(sz(1));
        }
        
        // phi
        double dphi = 2. * numerical::dPi / nr;
        for (int iphi = 0; iphi < nr; iphi++) {
            auto seq = Eigen::seqN(iphi, Eigen::fix<spectral::nPEM>, nr);
            coords_spz(seq, 1).fill(iphi * dphi);
        }
        
        // key
        std::string prefix = quadKey() + "__";
        int irec = 0;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            for (int iphi = 0; iphi < nr; iphi++) {
                std::stringstream sp;
                sp << prefix << ipnt << "_" << iphi;
                keys[irec++] = sp.str();
            }
        }
        return nrec;
    }
    
    // get boundary points
    const std::vector<int> &pointsOnWJB() const {
        return mPointsOnWJB;
    }
    
    /////////////////////////// data ///////////////////////////
private:
    // a copy of Element to compute incident stiffness (sInteriorStiff)
    // NOTE: we need independent copies of attenuation memory variables
    //       for the incident and the scattered waves
    const std::unique_ptr<Element> mIncidentElement;
    
    // point indecies on the injection boundary
    const std::vector<int> mPointsOnWJB;
    
    // injected displacement
    // use a pointer to preserve the "const" modifier of element in domain
    std::unique_ptr<vec_arD_CMatPP_RM> mInjectedDispl = nullptr;
    
    // STF
    std::unique_ptr<STF_NetCDF> mSTF = nullptr;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // stiffness from interior to exterior
    inline static vec_arD_CMatPP_RM sInteriorStiff;
    
    // displacement read from STF
    inline static typename STF_NetCDF::CTMatXND sInjectedDisplFlat =
    typename STF_NetCDF::CTMatXND(0, spectral::nPEM * ndim);
};

#endif /* ElementInteriorWJ_hpp */
