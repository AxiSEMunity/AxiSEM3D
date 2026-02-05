//
//  StationGroup.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/7/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station group

#ifndef StationGroup_hpp
#define StationGroup_hpp

#include "StationIO.hpp"
#include "StationSolid.hpp"
#include "StationFluid.hpp"

///////////////////// type inference /////////////////////
template <class StationT>
struct StationInference {
};

template <>
struct StationInference<StationSolid> {
    typedef channel::solid::ChannelOptions ChannelOptions;
    inline static const auto &gChannelMap = channel::solid::gChannelMap;
};

template <>
struct StationInference<StationFluid> {
    typedef channel::fluid::ChannelOptions ChannelOptions;
    inline static const auto &gChannelMap = channel::fluid::gChannelMap;
};


///////////////////// station group /////////////////////
template <class StationT>
class StationGroup {
public:
    // constructor
    StationGroup(const std::string &groupName, int numRecordSteps,
                 int sampleIntv, double tmin, double tmax,
                 int dumpIntv, channel::WavefieldCS wcs,
                 const std::vector<std::string> &userChannels,
                 std::unique_ptr<StationIO> &io):
    mGroupName(groupName), mNumRecordSteps(numRecordSteps),
    mSampleIntv(sampleIntv), mTmin(tmin), mTmax(tmax),
    mDumpIntv(std::min(dumpIntv, numRecordSteps)),
    mChannelOptions(wcs, userChannels), mIO(io.release()) {
        // nothing
    }
    
    // add a station
    void addStation(std::unique_ptr<StationT> &station) {
        mStations.push_back(std::move(station));
    }
    
    // initialize after adding all stations
    void initialize() {
        // channels
        int nch = (int)mChannelOptions.mStdChannels.size();
        std::vector<std::string> channels;
        channels.reserve(nch);
        for (int ich = 0; ich < nch; ich++) {
            const std::string &chstr =
            std::get<0>(StationInference<StationT>::gChannelMap.
                        at(mChannelOptions.mStdChannels[ich]));
            channels.push_back(chstr);
        }
        
        // station keys
        int nst = (int)mStations.size();
        std::vector<std::string> stKeys;
        stKeys.reserve(nst);
        for (int ist = 0; ist < nst; ist++) {
            stKeys.push_back(mStations[ist]->getKey());
        }
        
        // IO
        mIO->initialize(mGroupName, mNumRecordSteps, channels, stKeys);
        
        // buffers
        if (nst > 0) {
            // time buffer is not needed without stations
            mBufferTime.resize(mDumpIntv);
        }
        mBufferFields.resize(nst, nch, mDumpIntv);
        mBufferLine = 0;
        
        // stations
        for (const std::unique_ptr<StationT> &st: mStations) {
            st->setInGroup(mDumpIntv, mChannelOptions);
        }
    }
    
    // finalize after time loop
    void finalize() const {
        mIO->finalize();
    }
    
    // record at a time step
    void record(int timeStep, double time) {
        // interval and time window
        if (timeStep % mSampleIntv != 0 || time < mTmin || time > mTmax) {
            return;
        }
        
        // time
        if (mStations.size() > 0) {
            mBufferTime(mBufferLine) = time;
        }
        
        // stations
        for (const std::unique_ptr<StationT> &st: mStations) {
            st->record(mBufferLine, mChannelOptions);
        }
        
        // increment buffer line
        mBufferLine++;
        
        // dump to file
        if (mBufferLine == mDumpIntv) {
            dumpToFile();
        }
    }
    
    // dump to file
    void dumpToFile() {
        // redundant call at the end of the time loop
        if (mBufferLine == 0) {
            return;
        }
        
        // stations
        for (int ist = 0; ist < mStations.size(); ist++) {
            mStations[ist]->processReport(mBufferLine, mChannelOptions,
                                          ist, mBufferFields);
        }
        
        // IO
        // check zero station inside
        mIO->dumpToFile(mBufferTime, mBufferFields, mBufferLine);
        
        // rewind buffer line
        mBufferLine = 0;
    }
    
private:
    //////////////// const ////////////////
    // name
    const std::string mGroupName;
    
    // steps
    const int mNumRecordSteps;
    const int mSampleIntv;
    const double mTmin, mTmax;
    const int mDumpIntv;
    
    // channels
    const typename StationInference<StationT>::ChannelOptions mChannelOptions;
    
    // IO
    const std::unique_ptr<StationIO> mIO;
    
    //////////////// non-const ////////////////
    // buffer
    eigen::DColX mBufferTime = eigen::DColX(0);
    eigen::RTensor3 mBufferFields = eigen::RTensor3(0, 0, 0);
    int mBufferLine = 0;
    
    // stations
    std::vector<std::unique_ptr<StationT>> mStations;
};

#endif /* StationGroup_hpp */
