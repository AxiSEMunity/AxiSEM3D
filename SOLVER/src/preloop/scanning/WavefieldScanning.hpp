//
//  WavefieldScanning.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/28/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  wavefield scanning

#ifndef WavefieldScanning_hpp
#define WavefieldScanning_hpp

#include "numerical.hpp"
class SE_Model;
class Domain;

class WavefieldScanning {
    friend class Domain;
    
public:
    // setup
    static void setup(double dt, double period, int numTotalSteps,
                      const SE_Model &sem, Domain &domain);
    
private:
    std::string mFileName;
    numerical::Real mTolFourierH2 = 0.;
    numerical::Real mRelTolH2 = 0.;
    numerical::Real mAbsTolH2 = 0.;
    int mMaxNumPeaks = 0;
    int mScanningInterval = 0;
};

#endif /* WavefieldScanning_hpp */
