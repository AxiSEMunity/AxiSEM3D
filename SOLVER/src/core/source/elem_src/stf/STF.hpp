//
//  STF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/13/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of source-time function

#ifndef STF_hpp
#define STF_hpp

#include "numerical.hpp"
#include <memory>

class STF {
public:
    // padding mode
    enum class PaddingMode {None, FirstLast, Specified};
    
    // destructor
    virtual ~STF() = default;
    
    // get start time
    virtual double getStartTime() const = 0;
    
    // get value
    virtual numerical::Real getValue(double time) = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // verbose key width
    virtual int verboseKeyWidth() const {
        return 18;
    }
    
    ////////////////////////////// static //////////////////////////////
public:
    // build from inparam
    static std::unique_ptr<STF>
    buildInparam(int sindex, const std::string &sourceName, double dt);
};

#endif /* STF_hpp */
