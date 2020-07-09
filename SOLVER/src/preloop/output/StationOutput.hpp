//
//  StationOutput.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/25/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  stations: generator of StationGroup in core

#ifndef StationOutput_hpp
#define StationOutput_hpp

#include "channel.hpp"
#include <memory>

class SE_Model;
class Domain;

class StationOutput {
    enum class Format {AsciiStation, AsciiChannel, NetCDF};
    
public:
    // constructor
    StationOutput(const std::string &groupName, const std::string &fileName,
                  bool sourceCentered, bool ellipticity,
                  bool useDepth, bool depthSolid, bool undulatedGeometry,
                  channel::WavefieldCS wcs, bool fluid,
                  const std::vector<std::string> &userChannels,
                  Format format, double samplingPeriod, int bufferTimeSteps):
    mGroupName(groupName), mFileName(fileName),
    mSourceCentered(sourceCentered), mEllipticity(ellipticity),
    mUseDepth(useDepth), mDepthSolid(depthSolid),
    mUndulatedGeometry(undulatedGeometry),
    mWCS(wcs), mFluid(fluid), mUserChannels(userChannels),
    mFormat(format), mSamplingPeriod(samplingPeriod),
    mBufferTimeSteps(bufferTimeSteps) {
        // nothing
    }
    
private:
    // build from inparam
    static std::shared_ptr<const StationOutput>
    buildInparam(int gindex, double dt);
    
    // verbose
    std::string verbose(double dt, int numRecordSteps, int numStations) const;
    
public:
    // release stations to domain
    static void release(const SE_Model &sem, Domain &domain, double dt,
                        int nTotalSteps);
    
private:
    /////////////// raw data ///////////////
    // name
    const std::string mGroupName;
    
    // locations
    const std::string mFileName;
    const bool mSourceCentered;
    const bool mEllipticity;
    const bool mUseDepth;
    const bool mDepthSolid;
    const bool mUndulatedGeometry;
    
    // fields
    const channel::WavefieldCS mWCS;
    const bool mFluid;
    const std::vector<std::string> mUserChannels;
    
    // format
    const Format mFormat;
    const double mSamplingPeriod;
    const int mBufferTimeSteps;
};

#endif /* StationOutput_hpp */
