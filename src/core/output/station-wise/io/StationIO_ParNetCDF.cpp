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
    
    // return only if no station in group
    if (mNumStationsGlobal == 0) {
        return;
    }
    
    // filename
    const std::string &gdir = io::gOutputDirectory + "/stations/" + groupName;
    const std::string &fname = gdir + "/axisem3d_synthetics.nc";
    
    // gather stations
    std::vector<std::vector<std::string>> stKeysRanks;
    mpi::gather(stKeys, stKeysRanks, mRankWithMaxNumStations);
    
    // file pointer
    mNcFile = std::make_unique<NetCDF_Writer>();
    
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
        mNcFile->open(fname, true);
        
        ///////////////////// define variables /////////////////////
        mNcFile->defModeOn();
        
        // data time
        mVarID_Time = mNcFile->defineVariable("data_time", {
            {"dim_time", numRecordSteps}}, numerical::dErr);
        
        // data wave
        mVarID_Data = mNcFile->defineVariable("data_wave", {
            {"dim_station", stKeysAll.size()},
            {"dim_channel", channels.size()},
            {"dim_time", numRecordSteps}
        }, (numerical::Real)numerical::dErr);
        
        // list channel
        mNcFile->defineVariable("list_channel", {
            {"dim_channel", channels.size()},
            {"dim_channel_str_length", vector_tools::maxLength(channels)}
        }, (char)0);
        
        // list station
        mNcFile->defineVariable("list_station", {
            {"dim_station", stKeysAll.size()},
            {"dim_station_str_length", vector_tools::maxLength(stKeysAll)}
        }, (char)0);
        
        // end defining variables
        mNcFile->defModeOff();
        
        ///////////////////// write info /////////////////////
        // list channel
        for (int ich = 0; ich < channels.size(); ich++) {
            mNcFile->writeVariable("list_channel", channels[ich],
                                   {ich, 0}, {1, (int)channels[ich].size()});
        }
        
        // list station
        for (int ist = 0; ist < stKeysAll.size(); ist++) {
            mNcFile->writeVariable("list_station", stKeysAll[ist],
                                   {ist, 0}, {1, (int)stKeysAll[ist].size()});
        }
        
        // close serial
        mNcFile->close();
    }
    
    // bcast variable IDs
    // variable ID should not change upon re-opening
    mpi::bcast(mVarID_Time, mRankWithMaxNumStations);
    mpi::bcast(mVarID_Data, mRankWithMaxNumStations);
    
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
        mNcFile.reset(nullptr);
    }
    mFileLineTime = 0;
}

// dump to file
void StationIO_ParNetCDF::dumpToFile(const eigen::DColX &bufferTime,
                                     const eigen::RTensor3 &bufferFields,
                                     int bufferLine) {
    // write time
    if (mpi::rank() == mRankWithMaxNumStations) {
        mNcFile->writeVariable(mVarID_Time, "data_time", bufferTime,
                               {mFileLineTime}, {bufferLine});
    }
    
    // write data
    int nst = (int)bufferFields.dimensions()[0];
    if (nst > 0) {
        int nch = (int)bufferFields.dimensions()[1];
        
        // write to file
        if (bufferLine == bufferFields.dimension(2)) {
            // full buffer
            mNcFile->writeVariable(mVarID_Data, "data_wave", bufferFields,
                                   {mGlobalIndexFirstStation, 0, mFileLineTime},
                                   {nst, nch, bufferLine});
        } else {
            // must truncate by copy because time is the fastest varying
            // dimension; occuring only at the end of the simulation
            eigen::IArray3 loc = {0, 0, 0};
            eigen::IArray3 len = {nst, nch, bufferLine};
            Eigen::internal::set_is_malloc_allowed(true);
            eigen::RTensor3 timeTruncated = bufferFields.slice(loc, len);
            Eigen::internal::set_is_malloc_allowed(false);
            mNcFile->writeVariable(mVarID_Data, "data_wave", timeTruncated,
                                   {mGlobalIndexFirstStation, 0, mFileLineTime},
                                   {nst, nch, bufferLine});
        }
    }
    
    // flush
    if (mFlush) {
        mNcFile->flush();
    }
    
    // update record postion in file
    mFileLineTime += bufferLine;
}
