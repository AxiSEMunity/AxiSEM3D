//
//  AttBuilder.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/5/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  attenuation builder
//  generator of Attenuation in core

#ifndef AttBuilder_hpp
#define AttBuilder_hpp

#include "eigen_sem.hpp"
#include "Attenuation.hpp"
class ExodusMesh;

class AttBuilder {
public:
    // constructor
    AttBuilder(const ExodusMesh &exodusMesh, bool cg4, double dt);
    
    // verbose
    static std::string
    verbose(const std::unique_ptr<const AttBuilder> &attBuilder);
    
    // create attenuation
    std::unique_ptr<Attenuation>
    createAttenuation(const eigen::DMatXN &QKp, const eigen::DMatXN &QMu,
                      eigen::DMatXN &kp, eigen::DMatXN &mu,
                      const eigen::DRow4 &weightsCG4, bool elastic1D) const;
    
private:
    // set alpha, beta, gamma
    // NOTE: with SPECFEM legency deprecated, alpha, beta and gamma
    //       all become element-independent
    void setAlphaBetaGamma() const;
    
private:
    // from Exodus
    const double mFmin, mFmax, mFref;
    const int mNSLS;
    eigen::DColX mW, mY;
    
    // from inparam
    const bool mUseCG4;
    
    // dt
    const double mDt;
};

#endif /* AttBuilder_hpp */
