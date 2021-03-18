//
//  MessageRank.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/12/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1-to-1 mpi communication

#include "MessageRank.hpp"
#include "WindowSum.hpp"

// constructor
MessageRank::
MessageRank(int rankOther, const std::vector<MeshWindowSum> &meshWindowSums):
mRankOther(rankOther), mMeshWindowSums(meshWindowSums) {
    // allocate buffer
    allocateBuffer();
}

// allocate buffer
void MessageRank::allocateBuffer() {
    // allocate buffers
    int sizeComm = 0;
    for (auto iter = mMeshWindowSums.begin(); iter != mMeshWindowSums.end(); iter++) {
        for (const auto &ws: std::get<1>(*iter)) {
            sizeComm += ws->sizeComm();
        }
    }
    mBufferSend = eigen::RColX::Zero(sizeComm);
    mBufferRecv = eigen::RColX::Zero(sizeComm);
}

// gather from points
void MessageRank::gatherFromPointWindows() {
    int row = 0;
    for (auto iter = mMeshWindowSums.begin(); iter != mMeshWindowSums.end(); iter++) {
        for (auto iter = mMeshWindowSums.begin(); iter != mMeshWindowSums.end(); iter++) {
            for (const auto &ws: std::get<1>(*iter)) {
                ws->feedComm(mBufferSend, row);
            }
        }
    }
}

// scatter to points
void MessageRank::scatterToPointWindows() const {
    int row = 0;
    for (auto iter = mMeshWindowSums.begin(); iter != mMeshWindowSums.end(); iter++) {
        for (auto iter = mMeshWindowSums.begin(); iter != mMeshWindowSums.end(); iter++) {
            for (const auto &ws: std::get<1>(*iter)) {
                ws->extractComm(mBufferRecv, row);
            }
        }
    }
}

// check if a point exists on this domain boundary
bool MessageRank::
contains(const std::shared_ptr<const PointWindow> &target) const {
    // here we use a lambda expression
    // [&target] means that our lambda captures 'target' by ref
    return std::find_if
    (mMeshWindowSums.begin(), mMeshWindowSums.end(),
     [&target](const MeshWindowSum &mws) {
        return std::get<0>(mws) == target->getMeshTag();
    }) != mMeshWindowSums.end();
}
