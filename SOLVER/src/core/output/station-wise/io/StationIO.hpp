//
//  StationIO.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/23/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  IO for station output

#ifndef StationIO_hpp
#define StationIO_hpp

#include "eigen_station.hpp"

class StationIO {
public:
    // destructor
    virtual ~StationIO() = default;
    
    // initialize
    virtual void initialize(const std::string &groupName,
                            int numRecordSteps,
                            const std::vector<std::string> &channels,
                            const std::vector<std::string> &stKeys);
    
    // finalize
    virtual void finalize() = 0;
    
    // dump to file
    virtual void
    dumpToFile(const eigen::DColX &bufferTime,
               const eigen::RTensor3 &bufferFields,
               int bufferLine) = 0;
    
    // set flush
    void setFlush(bool flush) {
        mFlush = flush;
    }
    
protected:
    // create rank_station.info
    void createRankStation(const std::string &groupName,
                           const std::vector<std::string> &stKeys) const;
    
    //////////////// mpi global properties ////////////////
protected:
    // total number of stations
    int mNumStationsGlobal = 0;
    
    // mpi rank with maximum number of stations
    int mRankWithMaxNumStations = -1;
    
    // global index of first station
    int mGlobalIndexFirstStation = 0;
    
    // flush
    bool mFlush = false;
};

#endif /* StationIO_hpp */
