//
//  MessageRank.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/12/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1-to-1 mpi communication

#include "MessageRank.hpp"
#include "SolidPoint.hpp"
#include "FluidPoint.hpp"

// constructor
MessageRank::
MessageRank(int rankOther, const std::vector<MeshPoint> &meshPoints):
mRankOther(rankOther), mMeshPoints(meshPoints) {
    // allocate buffer
    allocateBuffer();
}

// allocate buffer
void MessageRank::allocateBuffer() {
    // allocate buffers
    int sizeComm = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        if (const auto &sp = std::get<1>(*iter)) {
            sizeComm += sp->sizeComm();
        }
        if (const auto &fp = std::get<2>(*iter)) {
            sizeComm += fp->sizeComm();
        }
    }
    mBufferSend = eigen::CColX::Zero(sizeComm);
    mBufferRecv = eigen::CColX::Zero(sizeComm);
}

// gather from points
void MessageRank::gatherFromPoints() {
    int row = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        if (const auto &sp = std::get<1>(*iter)) {
            sp->feedComm(mBufferSend, row);
        }
        if (const auto &fp = std::get<2>(*iter)) {
            fp->feedComm(mBufferSend, row);
        }
    }
}

// scatter to points
void MessageRank::scatterToPoints() const {
    int row = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        if (const auto &sp = std::get<1>(*iter)) {
            sp->extractComm(mBufferRecv, row);
        }
        if (const auto &fp = std::get<2>(*iter)) {
            fp->extractComm(mBufferRecv, row);
        }
    }
}

// check if a point exists on this domain boundary
bool MessageRank::
contains(const std::shared_ptr<const Point> &target) const {
    // here we use a lambda expression
    // [&target] means that our lambda captures 'target' by ref
    return std::find_if
    (mMeshPoints.begin(), mMeshPoints.end(),
     [&target](const MeshPoint &mpnt) {
        return std::get<0>(mpnt) == target->getMeshTag();
    }) != mMeshPoints.end();
}
