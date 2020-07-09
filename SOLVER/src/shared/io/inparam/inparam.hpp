//
//  inparam.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global input parameters

#ifndef inparam_hpp
#define inparam_hpp

#include "InparamYAML.hpp"

namespace inparam {
    // global input parameters
    extern InparamYAML gInparamModel;
    extern InparamYAML gInparamNr;
    extern InparamYAML gInparamSource;
    extern InparamYAML gInparamOutput;
    extern InparamYAML gInparamAdvanced;
    
    // setup
    void setup();
    
    // verbose
    std::string verbose();
}

#endif /* inparam_hpp */
