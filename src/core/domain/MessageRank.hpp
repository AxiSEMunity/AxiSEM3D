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
class Point;
class SolidPoint;
class FluidPoint;

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
    std::shared_ptr<SolidPoint>,
    std::shared_ptr<FluidPoint>> MeshPoint;
    
    // constructor
    MessageRank(int rankOther, const std::vector<MeshPoint> &meshPoints);
    
private:
    // allocate buffer
    void allocateBuffer();
    
public:
    // gather from points
    void gatherFromPoints();
    
    // scatter to points
    void scatterToPoints() const;
    
    // send buffer to the other rank
    void sendBuffer(MPI_Request &request) const {
        mpi::isend(mRankOther, mBufferSend, request);
    }
    
    // recv buffer from the other rank
    void recvBuffer(MPI_Request &request) {
        mpi::irecv(mRankOther, mBufferRecv, request);
    }
    
    // check if a point exists on this domain boundary
    bool contains(const std::shared_ptr<const Point> &target) const;
    
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
    eigen::CColX mBufferSend = eigen::CColX(0);
    eigen::CColX mBufferRecv = eigen::CColX(0);
};

#endif /* MessageRank_hpp */
