//
//  NrFieldAnalytical.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  analytical Nr(s,z)

#include "NrFieldAnalytical.hpp"
#include "geodesy.hpp"
#include "vector_tools.hpp"
#include "bstring.hpp"
#include "inparam.hpp"

// TODO: for reproducibility, edit sCodeID to indentify your code
std::string NrFieldAnalytical::sCodeID = "depth-dependent (AxiSEM3D default)";

// constructor
NrFieldAnalytical::NrFieldAnalytical() {
    // check code ID
    const std::string &codeID =
    inparam::gInparamNr.get<std::string>("analytical:code_ID");
    if (sCodeID != codeID) {
        throw
        std::runtime_error("NrFieldAnalytical::NrFieldAnalytical || "
                           "Inconsistent code ID's: || "
                           "in inparam.nr.yaml:       " + codeID + " ||"
                           "in NrFieldAnalytical.cpp: " + sCodeID);
    }
    
    // TODO: add your code here to initialize required data
    // below are data initialization for
    // sCodeID = "depth-dependent (AxiSEM3D default)"
    
    // read data
    mControlDepths = inparam::gInparamNr.getVector<double>
    ("analytical:depth_dependent_AxiSEM3D_default:control_depths");
    mControlNrs = inparam::gInparamNr.getVector<double>
    ("analytical:depth_dependent_AxiSEM3D_default:Nr_at_control_depths");
    // verify data
    if (mControlDepths.size() != mControlNrs.size()) {
        throw std::runtime_error("NrFieldAnalytical::NrFieldAnalytical || "
                                 "Depths and Nr's do not pair up. || "
                                 "Code ID: " + sCodeID);
    }
    if (mControlDepths.size() < 2) {
        throw std::runtime_error("NrFieldAnalytical::NrFieldAnalytical || "
                                 "Insufficient control depths, at least 2. || "
                                 "Code ID: " + sCodeID);
    }
    if (!vector_tools::isSortedUnique(mControlDepths)) {
        throw std::runtime_error("NrFieldAnalytical::NrFieldAnalytical || "
                                 "Depths must be ascendingly sorted. || "
                                 "Code ID: " + sCodeID);
    }
}

// get nr by (s, z)
eigen::IColX NrFieldAnalytical::
getNrAtPoints(const eigen::DMatX2_RM &sz) const {
    // input: coordinates (s, z, r, theta, depth) and mUserParameters
    [[maybe_unused]]
    const eigen::DColX &s = sz.col(0);
    const eigen::DColX &z = sz.col(1);
    const eigen::DMatX2_RM &rtheta = geodesy::sz2rtheta(sz, true);
    const eigen::DColX &r = rtheta.col(0);
    [[maybe_unused]]
    const eigen::DColX &theta = rtheta.col(1);
    [[maybe_unused]]
    const eigen::DColX &depth = geodesy::isCartesian() ?
    geodesy::getOuterRadius() - z.array() :
    geodesy::getOuterRadius() - r.array();
    
    // output: azimuthal dimension Nr
    eigen::IColX nr = eigen::IColX::Zero(sz.rows());
    
    
    ///////////////////////////////////////////////////////////////////
    ///////////////////// Only edit the box below /////////////////////
    ///////////////////////////////////////////////////////////////////
    
    // TODO: compute Nr from s, z, r, theta, depth and mUserParameters
    // NOTE: Add your own code in this box.
    //       If you never edited this part of code, the following example
    //       implements a depth-dependent Nr field with a depth-Nr profile
    //       initialized in the constructor.
    //       sCodeID = "depth-dependent (AxiSEM3D default)"
    
    // linear interpolation in the depth-Nr profile
    nr = depth.unaryExpr([this](double _dep) {
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        try {
            vector_tools::linearInterpSorted(this->mControlDepths, _dep,
                                             index0, index1,
                                             factor0, factor1);
        } catch (...) {
            throw std::runtime_error("NrFieldAnalytical::getNrAtPoints || "
                                     "Target depth is out of range. || "
                                     "Code ID: " + sCodeID);
        }
        return (int)round(this->mControlNrs[index0] * factor0 +
                          this->mControlNrs[index1] * factor1);
    });
    
    ///////////////////////////////////////////////////////////////////
    ///////////////////// Only edit the box above /////////////////////
    ///////////////////////////////////////////////////////////////////
    
    // return
    return nr;
}

// verbose
std::string NrFieldAnalytical::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxTitle("Nr(s,z)");
    ss << boxEquals(0, 20, "type", "ANALYTICAL");
    ss << boxEquals(0, 20, "code ID", sCodeID);
    {
        // TODO: verbose you model parameters here
        // below are data verbose for
        // sCodeID = "depth-dependent (AxiSEM3D default)"
        ss << boxEquals(0, 20, "control depths", mControlDepths, "=", true);
        ss << boxEquals(0, 20, "Nr at control depths", mControlNrs, "=", true);
    }
    ss << boxBaseline() << "\n\n";
    return ss.str();
}
