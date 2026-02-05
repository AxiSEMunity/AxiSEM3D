//
//  SolidFluidBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#include "SolidFluidBoundary.hpp"

// domain
#include "Messaging.hpp"
#include "SolidPoint.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"

// count info
std::map<std::string, int> SolidFluidBoundary::
countInfo(const Messaging &msg) const {
    std::map<std::string, int> countMap;
    for (const auto &sfc: mSFCs) {
        if (!msg.pointInSmallerRank(sfc->getSolidPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*sfc), 1);
        }
    }
    return countMap;
}
