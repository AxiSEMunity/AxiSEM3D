//
//  StationIO_ParNetCDF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  parallel NetCDF IO for station output

#include "StationIO_ParNetCDF.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include "vector_tools.hpp"

// initialize
void StationIO_ParNetCDF::initialize(const std::string &groupName,
                                     int numRecordSteps,
                                     const std::vector<std::string> &channels,
                                     const std::vector<std::string> &stKeys) {
    // finalize
    finalize();
    
    // base
    StationIO::initialize(groupName, numRecordSteps, channels, stKeys);
    if (mNumStationsGlobal == 0) {
        return;
    }
    
    // filename
    const std::string &gdir = io::gOutputDirectory + "/stations/" + groupName;
    const std::string &fname = gdir + "/axisem3d_synthetics.nc";
    
    // gather stations
    std::vector<std::vector<std::string>> stKeysRanks;
    mpi::gather(stKeys, stKeysRanks, mRankWithMaxNumStations);
    
    /////////////////////// only on max rank ///////////////////////
    if (mpi::rank() == mRankWithMaxNumStations) {
        // flatten keys
        std::vector<std::string> stKeysAll;
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            stKeysAll.insert(stKeysAll.end(),
                             stKeysRanks[irank].begin(),
                             stKeysRanks[irank].end());
        }
        
        // create and open
        mNcFile = std::make_unique<NetCDF_Writer>();
        mNcFile->open(fname, true);
        
        ///////////////////// define variables /////////////////////
        mNcFile->defModeOn();
        
        // time
        mVarID_Time = mNcFile->defineVariable("time_points", {
            {"dim_time", numRecordSteps}}, numerical::dErr);
        
        // data
        mVarID_Data = mNcFile->defineVariable("data_wave", {
            {"dim_time", numRecordSteps},
            {"dim_channel", channels.size()},
            {"dim_station", stKeysAll.size()}
        }, (numerical::Real)numerical::dErr);
        
        // channels
        mNcFile->defineVariable("channel_order", {
            {"dim_channel", channels.size()},
            {"dim_channel_str_length", vector_tools::maxLength(channels)}
        }, (char)0);
        
        // stations
        mNcFile->defineVariable("station_order", {
            {"dim_station", stKeysAll.size()},
            {"dim_station_str_length", vector_tools::maxLength(stKeysAll)}
        }, (char)0);
        
        // end defining variables
        mNcFile->defModeOff();
        
        ///////////////////// write info /////////////////////
        // write channels
        for (int ich = 0; ich < channels.size(); ich++) {
            mNcFile->writeVariable("channel_order", channels[ich],
                                   {ich, 0}, {1, (int)channels[ich].size()});
        }
        
        // write station keys
        for (int ist = 0; ist < stKeysAll.size(); ist++) {
            mNcFile->writeVariable("station_order", stKeysAll[ist],
                                   {ist, 0}, {1, (int)stKeysAll[ist].size()});
        }
        
        // close serial
        mNcFile->close();
    }
    
    // bcast variable IDs
    // variable ID should not change upon re-opening
    mpi::bcast(mVarID_Time);
    mpi::bcast(mVarID_Data);
    
    // open parallel
    mpi::barrier();
    mNcFile->openParallel(fname);
    
    // reset time line
    mFileLineTime = 0;
}

// finalize
void StationIO_ParNetCDF::finalize() {
    // close files
    if (mNcFile) {
        mNcFile->close();
        mNcFile = nullptr;
    }
    mFileLineTime = 0;
}

// dump to file
void StationIO_ParNetCDF::dumpToFile(const eigen::DColX &bufferTime,
                                     const eigen::RTensor3 &bufferFields,
                                     int bufferLine) {
    // no station
    int nst = (int)bufferFields.dimensions()[2];
    if (nst == 0) {
        return;
    }
    
    // no line
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    mNcFile->writeVariable(mVarID_Time, "time_points", bufferTime,
                           {mFileLineTime}, {bufferLine});
    
    // write data
    int nch = (int)bufferFields.dimensions()[1];
    mNcFile->writeVariable(mVarID_Data, "data", bufferFields,
                           {mFileLineTime, 0, mGlobalIndexFirstStation},
                           {bufferLine, nch, nst});
    
    // flush?
    // mNcFile->flush();
    
    // update record postion in file
    mFileLineTime += bufferLine;
}
