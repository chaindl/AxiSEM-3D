//
//  MessageRank.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/12/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  1-to-1 mpi communication

#ifndef MessageRank_hpp
#define MessageRank_hpp

// point
#include <memory>
class PointWindow;
class WindowSum;

// buffer
#include "eigen_generic.hpp"

// mpi
#include "mpi.hpp"

class MessageRank {
public:
    // a mesh point on mpi boundary
    // NOTE: it may contain a solid point or a fluid point or both
    //       std::map cannot be used because it changes the order of insertion
    typedef std::tuple<int,
    std::vector<std::shared_ptr<WindowSum>>> MeshPoint;
    
    // constructor
    MessageRank(int rankOther, const std::vector<MeshPoint> &MeshPoints);
    
private:
    // allocate buffer
    void setUpWindowSums();
    
public:
    void finalizeComms();

    // gather from points
    void gatherFromPointWindows();
    
    // scatter to points
    void scatterToPointWindows() const;
    
    // send buffer to the other rank
    void sendBuffer(MPI_Request &request) const {
        mpi::isend(mRankOther, mBufferSend, request);
    }
    
    // recv buffer from the other rank
    void recvBuffer(MPI_Request &request) {
        mpi::irecv(mRankOther, mBufferRecv, request);
    }
    
    // check if a point exists on this domain boundary
    bool contains(const std::shared_ptr<const PointWindow> &target) const;
    
    // get the other rank
    inline int getRankOther() const {
        return mRankOther;
    }
    
private:
    // the other rank to communicate with
    const int mRankOther;
    
    // my points involved in this 1-to-1 communication
    std::vector<MeshPoint> mMeshPoints;
    
    // buffers
    eigen::RColX mBufferSend = eigen::RColX(0);
    eigen::RColX mBufferRecv = eigen::RColX(0);
};

#endif /* MessageRank_hpp */
