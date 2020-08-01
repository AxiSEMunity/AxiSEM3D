//
//  ElementIO_ParNetCDF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  parallel NetCDF IO for element output

#ifndef ElementIO_ParNetCDF_hpp
#define ElementIO_ParNetCDF_hpp

#include "ElementIO.hpp"
#include "NetCDF_Writer.hpp"

class ElementIO_ParNetCDF: public ElementIO {
public:
    // initialize
    void initialize(const std::string &groupName,
                    int numRecordSteps,
                    const std::vector<std::string> &channels,
                    int npnts, const std::vector<int> &naGrid,
                    const eigen::IMatX4_RM &elemNaInfo,
                    const eigen::DMatXX_RM &elemCoords);
    
    // finalize
    void finalize();
    
    // dump to file
    void dumpToFile(const eigen::DColX &bufferTime,
                    const std::vector<eigen::RTensor5> &bufferFields,
                    int bufferLine);
    
private:
    //////////////////// const ////////////////////
    // file
    std::unique_ptr<NetCDF_Writer> mNcFile = nullptr;
    
    // variable id
    int mVarID_Time = -1;
    std::vector<int> mVarID_Data;
    
    // current line in time dimension
    int mFileLineTime = 0;
    
    // first element index on na-grid
    std::vector<int> mFirstElemNaGrid;
};

#endif /* ElementIO_ParNetCDF_hpp */
