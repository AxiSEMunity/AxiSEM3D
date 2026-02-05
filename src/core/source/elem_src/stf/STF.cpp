//
//  STF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of source-time function

#include "STF.hpp"
#include "inparam.hpp"
#include "GaussianSTF.hpp"
#include "StreamSTF.hpp"
#include "NetCDF_STF.hpp"

// build from inparam
std::unique_ptr<STF> STF::
buildInparam(int sindex, const std::string &sourceName, double dt) {
    // short alias
    const InparamYAML &gm = inparam::gInparamSource;
    const std::string &root = ("list_of_sources:[" + bstring::toString(sindex)
                               + "]:" + sourceName + ":source_time_function");
    
    // class name
    const std::string &className = gm.get<std::string>(root + ":class_name");
    
    // init class
    if (className == "GaussianSTF") {
        double hdur = gm.getWithBounds(root + ":half_duration", 0.);
        // use dt * 5. as the minimum (SPECFEM3D_GLOBE)
        hdur = std::max(hdur, dt * 5.);
        double decay = gm.getWithBounds(root + ":decay_factor", 0.);
        double shift = gm.get<double>(root + ":time_shift");
        int order = gm.getWithLimits<int>(root + ":use_derivative_integral", {
            {"ERF", -1},
            {"GAUSSIAN", 0},
            {"FIRST_DERIVATIVE", 1},
            {"RICKER", 2}});
        return std::make_unique<GaussianSTF>(hdur, decay, shift, order);
    } else if (className == "StreamSTF" || className == "NetCDF_STF") {
        // padding
        const std::string &padding = gm.get<std::string>(root + ":padding");
        StreamSTF::PaddingMode pm;
        double left = 0., right = 0.;
        if (padding == "NONE") {
            pm = StreamSTF::PaddingMode::None;
        } else if (padding == "FIRST_LAST") {
            pm = StreamSTF::PaddingMode::FirstLast;
        } else {
            pm = StreamSTF::PaddingMode::Specified;
            left = gm.get<double>(root + ":padding:[0]");
            right = gm.get<double>(root + ":padding:[1]");
        }
        
        // file type
        if (className == "StreamSTF") {
            const std::string &fname =
            gm.get<std::string>(root + ":ascii_data_file");
            return std::make_unique<StreamSTF>(fname, pm, left, right);
        } else {
            const std::string &fname =
            gm.get<std::string>(root + ":nc_data_file");
            const std::string &varTime =
            gm.get<std::string>(root + ":nc_var_times");
            const std::string &varData =
            gm.get<std::string>(root + ":nc_var_data");
            int chunkSize = gm.getWithOptions<int>
            (root + ":chunk_size", {{"NONE", std::numeric_limits<int>::max()}});
            return std::make_unique<NetCDF_STF>(fname, varTime, varData,
                                                chunkSize, pm, left, right);
        }
    } else {
        throw std::runtime_error("STF::buildInparam || "
                                 "Unknown class of source-time function."
                                 " || Source name: " + sourceName +
                                 " || STF class name: " + className);
    }
}
