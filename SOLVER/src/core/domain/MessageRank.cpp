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
#include <iostream>
// constructor
MessageRank::
MessageRank(int rankOther, const std::vector<MeshPoint> &meshPoints):
mRankOther(rankOther), mMeshPoints(meshPoints) {
    // allocate buffer
    setUpWindowSums();
}

void MessageRank::setUpWindowSums() {
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        for (const auto &ws: std::get<1>(*iter)) {
            ws->allocateComm();
        }
    }
    
    // allocate buffers
    int sizeComm = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        for (const auto &ws: std::get<1>(*iter)) {
            sizeComm += ws->sizeComm();
        }
    }
    mBufferSend = eigen::CColX::Zero(sizeComm);
    mBufferRecv = eigen::CColX::Zero(sizeComm);
}

void MessageRank::finalizeComms() {
    //nothing
}

// gather from points
void MessageRank::gatherFromPointWindows() {
    int row = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        for (const auto &ws: std::get<1>(*iter)) {
            ws->feedComm(mBufferSend, row);
        }
    }
}

// scatter to points
void MessageRank::scatterToPointWindows() const {
    int row = 0;
    for (auto iter = mMeshPoints.begin(); iter != mMeshPoints.end(); iter++) {
        for (const auto &ws: std::get<1>(*iter)) {
            ws->extractComm(mBufferRecv, row);
        }
    }
}

// check if a point exists on this domain boundary
bool MessageRank::
contains(const std::shared_ptr<const PointWindow> &target) const {
    // here we use a lambda expression
    // [&target] means that our lambda captures 'target' by ref
    return std::find_if
    (mMeshPoints.begin(), mMeshPoints.end(),
     [&target](const MeshPoint &mws) {
        return std::get<0>(mws) == target->getMeshTag();
    }) != mMeshPoints.end();
}
