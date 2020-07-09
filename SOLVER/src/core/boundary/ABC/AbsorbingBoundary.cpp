//
//  AbsorbingBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  absorbing boundary condition

#include "AbsorbingBoundary.hpp"

// domain
#include "Messaging.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"
#include "bstring.hpp"
#include "vector_tools.hpp"

// count info
std::map<std::string, int> AbsorbingBoundary::
countInfo(const Messaging &msg) const {
    std::map<std::string, int> countMap;
    // Clayton in solid
    for (const auto &clayton: mClaytonSolids) {
        if (!msg.pointInSmallerRank(clayton->getPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*clayton), 1);
        }
    }
    // Clayton in fluid
    for (const auto &clayton: mClaytonFluids) {
        if (!msg.pointInSmallerRank(clayton->getPoint())) {
            vector_tools::aggregate(countMap, bstring::typeName(*clayton), 1);
        }
    }
    // sponge in solid
    for (const auto &sponge: mSpongeSolids) {
        if (!msg.pointInSmallerRank(sponge->getPoint())) {
            vector_tools::aggregate(countMap, "SpongeSolid", 1);
        }
    }
    // sponge in fluid
    for (const auto &sponge: mSpongeFluids) {
        if (!msg.pointInSmallerRank(sponge->getPoint())) {
            vector_tools::aggregate(countMap, "SpongeFluid", 1);
        }
    }
    return countMap;
}
