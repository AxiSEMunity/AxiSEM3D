//
//  FluidElement.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid spectral element

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include "Element.hpp"

// point
#include <array>
class FluidPoint;

// material
#include "Acoustic.hpp"

// output
#include "channel.hpp"

class FluidElement: public Element {
public:
    // constructor
    FluidElement(int quadTag, std::unique_ptr<const GradQuad> &grad,
                 std::unique_ptr<const PRT> &prt,
                 std::unique_ptr<const Acoustic> &acoustic,
                 const std::array<std::shared_ptr<FluidPoint>,
                 spectral::nPEM> &points);
    
    // copy constructor
    FluidElement(const FluidElement &other);
    
private:
    // construct derived
    void constructDerived();
    
public:
    // type info
    std::string typeInfo() const;
    
    // medium info
    std::string mediumInfo() const {
        return "FLUID";
    }
    
    
    /////////////////////////// point ///////////////////////////
    // get point
    Point &getPoint(int ipnt) const;
    
    
    /////////////////////////// time loop ///////////////////////////
    // collect displacement from points
    virtual void
    collectDisplFromPoints(eigen::vec_ar1_CMatPP_RM &displElem) const;
    
    // displacement to stiffness
    void displToStiff(const eigen::vec_ar1_CMatPP_RM &displElem,
                      eigen::vec_ar1_CMatPP_RM &stiffElem) const;
    
    // add stiffness to points
    // allow a derived class to change stiffElem (no const)
    virtual void
    addStiffToPoints(eigen::vec_ar1_CMatPP_RM &stiffElem) const;
    
    // compute stiffness term
    void computeStiff() const;
    
    // measure cost
    double measure(int count) const;
    
    
    /////////////////////////// source ///////////////////////////
    // prepare pressure source
    void preparePressureSource() const;
    
    // add pressure source
    void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure) const;
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare wavefield output
    void prepareWavefieldOutput(const channel::fluid::ChannelOptions &chops,
                                bool enforceCoordTransform);
    
    // chi field
    void getChiField(eigen::CMatXN &chi) const;
    
    // displ field
    void getDisplField(eigen::CMatXN3 &displ) const {
        getDisplField(displ, displInRTZ());
    }
    
    // displ field
    void getDisplField(eigen::CMatXN3 &displ, bool needRTZ) const;
    
    // pressure field
    void getPressureField(eigen::CMatXN &pressure) const;
    
    // delta field
    void getDeltaField(eigen::CMatXN &delta) const;
    
    // displ crd
    inline bool displInRTZ() const {
        return (bool)mPRT;
    }
    
    
private:
    // material
    const std::unique_ptr<const Acoustic> mAcoustic;
    
    // 1D element in Fourier space
    const bool mInFourier;
    
    // points
    std::array<std::shared_ptr<FluidPoint>, spectral::nPEM> mPoints;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // expand workspace
    static void expandWorkspace(int maxNr) {
        int maxNu_1 = maxNr / 2 + 1;
        // Fourier
        sDisplSpherical_FR.resize(maxNu_1);
        sStrainSpherical_FR.resize(maxNu_1);
        sStrainUndulated_FR.resize(maxNu_1);
        sStressUndulated_FR.resize(maxNu_1);
        sStressSpherical_FR.resize(maxNu_1);
        sStiffSpherical_FR.resize(maxNu_1);
        // cardinal
        sStrainSpherical_CD.resize(maxNr, spectral::nPEM * 3);
        sStrainUndulated_CD.resize(maxNr, spectral::nPEM * 3);
        sStressUndulated_CD.resize(maxNr, spectral::nPEM * 3);
        sStressSpherical_CD.resize(maxNr, spectral::nPEM * 3);
    }
    
    // workspace
    // Fourier
    inline static eigen::vec_ar1_CMatPP_RM sDisplSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStrainSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStrainUndulated_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStressUndulated_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStressSpherical_FR;
    inline static eigen::vec_ar1_CMatPP_RM sStiffSpherical_FR;
    // cardinal
    inline static eigen::RMatXN3 sStrainSpherical_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStrainUndulated_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStressUndulated_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
    inline static eigen::RMatXN3 sStressSpherical_CD =
    eigen::RMatXN3(0, spectral::nPEM * 3);
};

#endif /* FluidElement_hpp */
