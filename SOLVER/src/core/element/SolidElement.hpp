//
//  SolidElement.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/2/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid spectral element

#ifndef SolidElement_hpp
#define SolidElement_hpp

#include "Element.hpp"

// point
#include <array>
class SolidPoint;

// material
#include "Elastic.hpp"

// output
#include "channel.hpp"

class SolidElement: public Element {
public:
    // constructor
    SolidElement(int quadTag, std::unique_ptr<const GradQuad> &grad,
                 std::unique_ptr<const PRT> &prt,
                 std::unique_ptr<const Elastic> &elastic,
                 const std::array<std::shared_ptr<SolidPoint>,
                 spectral::nPEM> &points);
    
    // copy constructor
    SolidElement(const SolidElement &other);
    
private:
    // construct derived
    void constructDerived();
    
public:
    // type info
    std::string typeInfo() const;
    
    // medium info
    std::string mediumInfo() const {
        return "SOLID";
    }
    
    
    /////////////////////////// point ///////////////////////////
    // get point
    Point &getPoint(int ipnt) const;
    
    
    /////////////////////////// time loop ///////////////////////////
    // collect displacement from points
    virtual void
    collectDisplFromPoints(eigen::vec_ar3_CMatPP_RM &displElem) const;
    
    // displacement to stiffness
    void displToStiff(const eigen::vec_ar3_CMatPP_RM &displElem,
                      eigen::vec_ar3_CMatPP_RM &stiffElem) const;
    
    // add stiffness to points
    // allow a derived class to change stiffElem (no const)
    virtual void
    addStiffToPoints(eigen::vec_ar3_CMatPP_RM &stiffElem) const;
    
    // compute stiffness term
    void computeStiff() const;
    
    // measure cost
    double measure(int count) const;
    
    
    /////////////////////////// source ///////////////////////////
    // prepare force source (force given in SPZ)
    void prepareForceSource() const;
    
    // add force source (force given in SPZ)
    void addForceSource(const eigen::CMatXN3 &force,
                        int nu_1_force) const;
    
    // prepare moment source (moment tensor given in SPZ)
    void prepareMomentSource() const;
    
    // add moment source (moment tensor given in SPZ)
    void addMomentSource(const eigen::CMatXN6 &moment,
                         int nu_1_moment) const;
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare wavefield output
    void prepareWavefieldOutput(const channel::solid::ChannelOptions &chops,
                                bool enforceCoordTransform);
    
    // displ field
    void getDisplField(eigen::CMatXN3 &displ) const {
        getDisplField(displ, displInRTZ());
    }
    
    // displ field
    void getDisplField(eigen::CMatXN3 &displ, bool needRTZ) const;
    
    // nabla field
    void getNablaField(eigen::CMatXN9 &nabla) const {
        getNablaField(nabla, nablaInRTZ());
    }
    
    // nabla field
    void getNablaField(eigen::CMatXN9 &nabla, bool needRTZ) const;
    
    // strain field
    void getStrainField(eigen::CMatXN6 &strain) const {
        getStrainField(strain, strainInRTZ());
    }
    
    // strain field
    void getStrainField(eigen::CMatXN6 &strain, bool needRTZ) const;
    
    // curl field
    void getCurlField(eigen::CMatXN3 &curl) const {
        getCurlField(curl, curlInRTZ());
    }
    
    // curl field
    void getCurlField(eigen::CMatXN3 &curl, bool needRTZ) const;
    
    // stress field
    void getStressField(eigen::CMatXN6 &stress) const {
        getStressField(stress, stressInRTZ());
    }
    
    // stress field
    void getStressField(eigen::CMatXN6 &stress, bool needRTZ) const;
    
    // displ crd
    inline bool displInRTZ() const {
        return false;
    }
    
    // nabla crd
    inline bool nablaInRTZ() const {
        return (bool)mPRT;
    }
    
    // strain crd
    inline bool strainInRTZ() const {
        return nablaInRTZ();
    }
    
    // curl crd
    inline bool curlInRTZ() const {
        return nablaInRTZ();
    }
    
    // stress crd
    inline bool stressInRTZ() const {
        return strainInRTZ() || mElastic->inRTZ();
    }
    
    
private:
    // material
    const std::unique_ptr<const Elastic> mElastic;
    
    // 1D element in Fourier space
    const bool mInFourier;
    
    // points
    std::array<std::shared_ptr<SolidPoint>, spectral::nPEM> mPoints;
    
    // stress buffer
    // stress cannot be recomputed with attenuation
    const std::unique_ptr<eigen::vec_ar6_CMatPP_RM> mStressBuffer =
    std::make_unique<eigen::vec_ar6_CMatPP_RM>();
    
    
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
        sStrainSpherical_CD.resize(maxNr, spectral::nPEM * 9);
        sStrainUndulated_CD.resize(maxNr, spectral::nPEM * 6);
        sStressUndulated_CD.resize(maxNr, spectral::nPEM * 6);
        sStressSpherical_CD.resize(maxNr, spectral::nPEM * 9);
    }
    
    // workspace
    // Fourier
    inline static eigen::vec_ar3_CMatPP_RM sDisplSpherical_FR;
    inline static eigen::vec_ar9_CMatPP_RM sStrainSpherical_FR;
    inline static eigen::vec_ar6_CMatPP_RM sStrainUndulated_FR;
    inline static eigen::vec_ar6_CMatPP_RM sStressUndulated_FR;
    inline static eigen::vec_ar9_CMatPP_RM sStressSpherical_FR;
    inline static eigen::vec_ar3_CMatPP_RM sStiffSpherical_FR;
    // cardinal
    inline static eigen::RMatXN9 sStrainSpherical_CD =
    eigen::RMatXN9(0, spectral::nPEM * 9);
    inline static eigen::RMatXN6 sStrainUndulated_CD =
    eigen::RMatXN6(0, spectral::nPEM * 6);
    inline static eigen::RMatXN6 sStressUndulated_CD =
    eigen::RMatXN6(0, spectral::nPEM * 6);
    inline static eigen::RMatXN9 sStressSpherical_CD =
    eigen::RMatXN9(0, spectral::nPEM * 9);
};

#endif /* SolidElement_hpp */
