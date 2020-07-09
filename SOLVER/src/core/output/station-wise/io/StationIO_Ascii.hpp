//
//  StationIO_Ascii.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/23/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Ascii IO for station output
//  only for a small number of stations

#ifndef StationIO_Ascii_hpp
#define StationIO_Ascii_hpp

#include "StationIO.hpp"
#include <memory>

class StationIO_Ascii: public StationIO {
public:
    // constructor
    StationIO_Ascii(bool stationCentric):
    mStationCentric(stationCentric) {
        // nothing
    }
    
    // initialize
    void initialize(const std::string &groupName,
                    int numRecordSteps,
                    const std::vector<std::string> &channels,
                    const std::vector<std::string> &stKeys);
    
    // finalize
    void finalize();
    
    // dump to file
    void dumpToFile(const eigen::DColX &bufferTime,
                    const eigen::RTensor3 &bufferFields,
                    int bufferLine);
    
private:
    // create file stream
    void createFileStream(const std::string &varKey, int precision);
    
    // station centric
    const bool mStationCentric;
    
    // file streams
    std::vector<std::unique_ptr<std::ofstream>> mFileStreams;
};

#endif /* StationIO_Ascii_hpp */
