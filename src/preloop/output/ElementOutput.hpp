//
//  ElementOutput.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output: generator of ElementOpGroup in core

#ifndef ElementOutput_hpp
#define ElementOutput_hpp

#include "channel.hpp"
#include <memory>

class SE_Model;
class Domain;

class ElementOutput {
public:
    // constructor
    ElementOutput(const std::string &groupName,
                  double minR, double maxR, double minZ, double maxZ,
                  int edgeDim, double edgeCoord, const std::vector<int> &ipols,
                  const std::vector<double> &phis,
                  const std::vector<double> &lats,
                  const std::vector<double> &lons, int naSpace,
                  channel::WavefieldCS wcs, bool fluid,
                  const std::vector<std::string> &userChannels,
                  double samplingPeriod, double tmin, double tmax,
                  int bufferSize, bool flush):
    mGroupName(groupName),
    mMinR(minR), mMaxR(maxR), mMinZ(minZ), mMaxZ(maxZ),
    mEdgeDim(edgeDim), mEdgeCoord(edgeCoord), mIPols(ipols),
    mPhis(phis), mLats(lats), mLons(lons), mNaSpace(naSpace),
    mWCS(wcs), mFluid(fluid), mUserChannels(userChannels),
    mSamplingPeriod(samplingPeriod), mTmin(tmin), mTmax(tmax),
    mBufferSize(bufferSize), mFlush(flush) {
        // nothing
    }
    
private:
    // build from inparam
    static std::shared_ptr<const ElementOutput>
    buildInparam(int gindex, double dt, double tmin_simu, double tmax_simu);
    
    // verbose
    std::string verbose(double dt, int numRecordSteps, int npnts,
                        int nphis, int numElements) const;
    
public:
    // release ElementOp to domain
    static void release(const SE_Model &sem, Domain &domain, double dt,
                        double tmin_simu, double tmax_simu,
                        int nTotalSteps, double distTol);
    
private:
    /////////////// raw data ///////////////
    // name
    const std::string mGroupName;
    
    // elements
    double mMinR, mMaxR;
    double mMinZ, mMaxZ;
    
    // inplane
    const int mEdgeDim;
    const double mEdgeCoord;
    const std::vector<int> mIPols;
    
    // azimuth
    const std::vector<double> mPhis;
    const std::vector<double> mLats;
    const std::vector<double> mLons;
    const int mNaSpace;
    
    // fields
    const channel::WavefieldCS mWCS;
    const bool mFluid;
    const std::vector<std::string> mUserChannels;
    
    // temporal
    const double mSamplingPeriod;
    const double mTmin, mTmax;
    
    // file
    const int mBufferSize;
    const bool mFlush;
};

#endif /* ElementOutput_hpp */
