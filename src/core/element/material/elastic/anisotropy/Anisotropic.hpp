//
//  Anisotropic.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/28/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  anisotropic

#ifndef Anisotropic_hpp
#define Anisotropic_hpp

#include "Elastic.hpp"

class Anisotropic: public Elastic {
public:
    // 1D constructor
    Anisotropic(std::unique_ptr<Attenuation> &attenuation,
                // C1
                const eigen::DMatPP_RM &C11,
                const eigen::DMatPP_RM &C12,
                const eigen::DMatPP_RM &C13,
                const eigen::DMatPP_RM &C14,
                const eigen::DMatPP_RM &C15,
                const eigen::DMatPP_RM &C16,
                // C2
                const eigen::DMatPP_RM &C22,
                const eigen::DMatPP_RM &C23,
                const eigen::DMatPP_RM &C24,
                const eigen::DMatPP_RM &C25,
                const eigen::DMatPP_RM &C26,
                // C3
                const eigen::DMatPP_RM &C33,
                const eigen::DMatPP_RM &C34,
                const eigen::DMatPP_RM &C35,
                const eigen::DMatPP_RM &C36,
                // C4
                const eigen::DMatPP_RM &C44,
                const eigen::DMatPP_RM &C45,
                const eigen::DMatPP_RM &C46,
                // C5
                const eigen::DMatPP_RM &C55,
                const eigen::DMatPP_RM &C56,
                // C6
                const eigen::DMatPP_RM &C66):
    Elastic(true, attenuation),
    mC11(C11), mC12(C12), mC13(C13), mC14(C14), mC15(C15), mC16(C16),
    mC22(C22), mC23(C23), mC24(C24), mC25(C25), mC26(C26),
    mC33(C33), mC34(C34), mC35(C35), mC36(C36),
    mC44(C44), mC45(C45), mC46(C46),
    mC55(C55), mC56(C56),
    mC66(C66) {
        // nothing
    }
    
    // 3D constructor
    Anisotropic(std::unique_ptr<Attenuation> &attenuation,
                // C1
                const eigen::DMatXN &C11,
                const eigen::DMatXN &C12,
                const eigen::DMatXN &C13,
                const eigen::DMatXN &C14,
                const eigen::DMatXN &C15,
                const eigen::DMatXN &C16,
                // C2
                const eigen::DMatXN &C22,
                const eigen::DMatXN &C23,
                const eigen::DMatXN &C24,
                const eigen::DMatXN &C25,
                const eigen::DMatXN &C26,
                // C3
                const eigen::DMatXN &C33,
                const eigen::DMatXN &C34,
                const eigen::DMatXN &C35,
                const eigen::DMatXN &C36,
                // C4
                const eigen::DMatXN &C44,
                const eigen::DMatXN &C45,
                const eigen::DMatXN &C46,
                // C5
                const eigen::DMatXN &C55,
                const eigen::DMatXN &C56,
                // C6
                const eigen::DMatXN &C66):
    Elastic(false, attenuation),
    mC11(C11), mC12(C12), mC13(C13), mC14(C14), mC15(C15), mC16(C16),
    mC22(C22), mC23(C23), mC24(C24), mC25(C25), mC26(C26),
    mC33(C33), mC34(C34), mC35(C35), mC36(C36),
    mC44(C44), mC45(C45), mC46(C46),
    mC55(C55), mC56(C56),
    mC66(C66) {
        // nothing
    }
    
    // copy constructor
    Anisotropic(const Anisotropic &other):
    Elastic(other),
    // C1
    mC11(other.mC11), mC12(other.mC12), mC13(other.mC13), mC14(other.mC14),
    mC15(other.mC15), mC16(other.mC16),
    // C2
    mC22(other.mC22), mC23(other.mC23), mC24(other.mC24), mC25(other.mC25),
    mC26(other.mC26),
    // C3
    mC33(other.mC33), mC34(other.mC34), mC35(other.mC35), mC36(other.mC36),
    // C4
    mC44(other.mC44), mC45(other.mC45), mC46(other.mC46),
    // C5
    mC55(other.mC55), mC56(other.mC56),
    // C6
    mC66(other.mC66) {
        // nothing
    }
    
    // clone for copy constructor
    std::unique_ptr<Elastic> clone() const {
        return std::make_unique<Anisotropic>(*this);
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool elemInFourier) const {
        // base
        Elastic::checkCompatibility(nr, elemInFourier);
        
        // parameters
        mC11.checkCompatibility(m1D, nr, elemInFourier, "mC11", "Anisotropic");
    }
    
    
    ///////////////////////// strain to stress /////////////////////////
    // RTZ coordinates
    bool inRTZ() const {
        return true;
    }
    
    // strain => stress in Fourier space
    void strainToStress_FR(const eigen::vec_ar6_CMatPP_RM &strain,
                           eigen::vec_ar6_CMatPP_RM &stress,
                           int nu_1) const {
        // elasticity
        for (int alpha = 0; alpha < nu_1; alpha++) {
            strainToStress<CaseFA::_1D_FR>
            (strain, stress, alpha,
             mC11, mC12, mC13, mC14, mC15, mC16,
             mC22, mC23, mC24, mC25, mC26,
             mC33, mC34, mC35, mC36,
             mC44, mC45, mC46,
             mC55, mC56,
             mC66);
        }
        
        // attenuation
        applyAttenuation(strain, stress, nu_1);
    }
    
    // strain => stress in cardinal space
    void strainToStress_CD(const eigen::RMatXN6 &strain,
                           eigen::RMatXN6 &stress,
                           int nr) const {
        // elasticity
        if (m1D) {
            strainToStress<CaseFA::_1D_CD>
            (strain, stress, nr,
             mC11, mC12, mC13, mC14, mC15, mC16,
             mC22, mC23, mC24, mC25, mC26,
             mC33, mC34, mC35, mC36,
             mC44, mC45, mC46,
             mC55, mC56,
             mC66);
        } else {
            strainToStress<CaseFA::_3D_CD>
            (strain, stress, nr,
             mC11, mC12, mC13, mC14, mC15, mC16,
             mC22, mC23, mC24, mC25, mC26,
             mC33, mC34, mC35, mC36,
             mC44, mC45, mC46,
             mC55, mC56,
             mC66);
        }
        
        // attenuation
        applyAttenuation(strain, stress, nr);
    }
    
    
    //////////////////////// properties ////////////////////////
private:
    // Cijkl
    const faN::PropertyN mC11, mC12, mC13, mC14, mC15, mC16;
    const faN::PropertyN mC22, mC23, mC24, mC25, mC26;
    const faN::PropertyN mC33, mC34, mC35, mC36;
    const faN::PropertyN mC44, mC45, mC46;
    const faN::PropertyN mC55, mC56;
    const faN::PropertyN mC66;
    
    
    /////////////////////////////////////////////////////////////////////
    /////////////////// template-based implementation ///////////////////
    /////////////////////////////////////////////////////////////////////
    
    template <CaseFA CASE, class FMat6, class PMat> static void
    strainToStress(const FMat6 &strain, FMat6 &stress, int nx,
                   // C1
                   const PMat &C11,
                   const PMat &C12,
                   const PMat &C13,
                   const PMat &C14,
                   const PMat &C15,
                   const PMat &C16,
                   // C2
                   const PMat &C22,
                   const PMat &C23,
                   const PMat &C24,
                   const PMat &C25,
                   const PMat &C26,
                   // C3
                   const PMat &C33,
                   const PMat &C34,
                   const PMat &C35,
                   const PMat &C36,
                   // C4
                   const PMat &C44,
                   const PMat &C45,
                   const PMat &C46,
                   // C5
                   const PMat &C55,
                   const PMat &C56,
                   // C6
                   const PMat &C66) {
        typedef faN::FieldArithmeticN FA;
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 0, strain,
                                 0, C11,
                                 1, C12,
                                 2, C13,
                                 3, C14,
                                 4, C15,
                                 5, C16);
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 1, strain,
                                 0, C12,
                                 1, C22,
                                 2, C23,
                                 3, C24,
                                 4, C25,
                                 5, C26);
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 2, strain,
                                 0, C13,
                                 1, C23,
                                 2, C33,
                                 3, C34,
                                 4, C35,
                                 5, C36);
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 3, strain,
                                 0, C14,
                                 1, C24,
                                 2, C34,
                                 3, C44,
                                 4, C45,
                                 5, C46);
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 4, strain,
                                 0, C15,
                                 1, C25,
                                 2, C35,
                                 3, C45,
                                 4, C55,
                                 5, C56);
        FA::F1xP1_to_F6xP6<CASE>(nx, stress, 5, strain,
                                 0, C16,
                                 1, C26,
                                 2, C36,
                                 3, C46,
                                 4, C56,
                                 5, C66);
    }
};

#endif /* Anisotropic_hpp */
