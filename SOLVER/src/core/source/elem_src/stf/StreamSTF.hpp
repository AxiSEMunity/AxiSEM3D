//
//  StreamSTF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  stream source-time function

#ifndef StreamSTF_hpp
#define StreamSTF_hpp

#include "STF.hpp"
#include <vector>

class StreamSTF: public STF {
public:
    // constructor
    StreamSTF(const std::string &fileName, PaddingMode padding,
              numerical::Real left, numerical::Real right);
    
    // get start time
    double getStartTime() const {
        return mTimes.front();
    }
    
    // get value
    numerical::Real getValue(double time);
    
    // verbose
    std::string verbose() const;
    
private:
    // file
    const std::string mFileName;
    
    // data
    std::vector<double> mTimes;
    std::vector<numerical::Real> mData;
    
    // padding
    const bool mPadding;
    numerical::Real mLeftPadding = 0.;
    numerical::Real mRightPadding = 0.;
};

#endif /* StreamSTF_hpp */
