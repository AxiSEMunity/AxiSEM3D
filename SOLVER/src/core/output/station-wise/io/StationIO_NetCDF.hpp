//
//  StationIO_NetCDF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  serial NetCDF IO for station output

#ifndef StationIO_NetCDF_hpp
#define StationIO_NetCDF_hpp

#include "StationIO.hpp"
#include "NetCDF_Writer.hpp"

class StationIO_NetCDF: public StationIO {
public:
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
    //////////////////// const ////////////////////
    // file
    std::unique_ptr<NetCDF_Writer> mNcFile = nullptr;
    
    // variable id
    int mVarID_Time = -1;
    int mVarID_Data = -1;
    
    // current line in time dimension
    int mFileLineTime = 0;
};

#endif /* StationIO_NetCDF_hpp */
