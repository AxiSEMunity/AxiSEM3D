//
//  NrFieldConstant.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  constant Nr(s,z)

#include "NrFieldConstant.hpp"
#include "bstring.hpp"

// get nr by (s, z)
eigen::IColX NrFieldConstant::getNrAtPoints(const eigen::DMatX2_RM &sz) const {
    return eigen::IColX::Constant(sz.rows(), mNr);
}

// verbose
std::string NrFieldConstant::verbose() const {
    std::stringstream ss;
    ss << bstring::boxTitle("Nr(s,z)");
    ss << bstring::boxEquals(0, 5, "type", "CONSTANT");
    ss << bstring::boxEquals(0, 5, "value", mNr);
    ss << bstring::boxBaseline() << "\n\n";
    return ss.str();
}
