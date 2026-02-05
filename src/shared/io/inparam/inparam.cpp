//
//  inparam.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  global input parameters

#include "inparam.hpp"
#include "io.hpp"

namespace inparam {
    // global input parameters
    InparamYAML gInparamModel("Model");
    InparamYAML gInparamNr("Nr");
    InparamYAML gInparamSource("Source");
    InparamYAML gInparamOutput("Output");
    InparamYAML gInparamAdvanced("Advanced");
    
    // setup
    void setup() {
        // parse
        gInparamModel.parse(io::gInputDirectory + "/inparam.model.yaml");
        gInparamNr.parse(io::gInputDirectory + "/inparam.nr.yaml");
        gInparamSource.parse(io::gInputDirectory + "/inparam.source.yaml");
        gInparamOutput.parse(io::gInputDirectory + "/inparam.output.yaml");
        gInparamAdvanced.parse(io::gInputDirectory + "/inparam.advanced.yaml");
    }
    
    // verbose
    std::string verbose() {
        std::stringstream ss;
        ss << bstring::boxTitle("Parameters");
        ss << gInparamModel.verbose();
        ss << gInparamNr.verbose();
        ss << gInparamSource.verbose();
        ss << gInparamOutput.verbose();
        ss << gInparamAdvanced.verbose();
        ss << bstring::boxBaseline() << "\n\n";
        return ss.str();
    }
}
