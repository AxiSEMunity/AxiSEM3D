//
//  ElementIO_NetCDF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  serial NetCDF IO for element output

#include "ElementIO_NetCDF.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include "vector_tools.hpp"

// initialize
void ElementIO_NetCDF::initialize(const std::string &groupName,
                                  int numRecordSteps,
                                  const std::vector<std::string> &channels,
                                  int npnts, const std::vector<int> &naGrid,
                                  const eigen::IMatX4_RM &elemNaInfo,
                                  const eigen::DMatXX_RM &elemCoords) {
    // finalize
    finalize();
    
    // base
    ElementIO::initialize(groupName, numRecordSteps, channels,
                          npnts, naGrid, elemNaInfo, elemCoords);
    
    
    ///////////////////////////////////////////////////////////////////////
    // the fifth column: element index in data after concat
    // gather elements
    std::vector<eigen::IMatX4_RM> elemNaInfoRanks;
    mpi::gatherEigen(elemNaInfo, elemNaInfoRanks, MPI_ALL);
    
    // form dict of naGrid for fast search
    std::map<int, int> naGridIndexDict;
    for (int inag = 0; inag < naGrid.size(); inag++) {
        naGridIndexDict[naGrid[inag]] = inag;
    }
    
    // first element index on na-grid
    std::vector<std::vector<int>> firstElemNaGridRanks;
    firstElemNaGridRanks.assign
    (mpi::nproc(), std::vector<int>(naGrid.size(), 0));
    // loop over ranks
    for (int irank = 1; irank < mpi::nproc(); irank++) {
        int nelocal = (int)elemNaInfoRanks[irank - 1].rows();
        // count elements on this rank
        for (int ielem = 0; ielem < nelocal; ielem++) {
            // four columns: tag, actual na, grid na, element index in data
            int nag = elemNaInfoRanks[irank - 1](ielem, 2);
            int nagIndex = naGridIndexDict[nag];
            firstElemNaGridRanks[irank][nagIndex]++;
        }
        // add counts on previous ranks
        for (int inag = 0; inag <naGrid.size(); inag++) {
            firstElemNaGridRanks[irank][inag] +=
            firstElemNaGridRanks[irank - 1][inag];
        }
    }
    
    // the fifth column
    int nelem = (int)elemNaInfo.rows();
    eigen::IMatX5_RM elemNaInfo5(nelem, 5);
    elemNaInfo5.leftCols(4) = elemNaInfo;
    for (int ielem = 0; ielem < nelem; ielem++) {
        int nag = elemNaInfo5(ielem, 2);
        int nagIndex = naGridIndexDict[nag];
        elemNaInfo5(ielem, 4) = elemNaInfo5(ielem, 3) +
        firstElemNaGridRanks[mpi::rank()][nagIndex];
    }
    ///////////////////////////////////////////////////////////////////////
    
    // nothing locally without elements
    if (nelem == 0) {
        return;
    }
    
    // filename
    const std::string &gdir = io::gOutputDirectory + "/elements/" + groupName;
    const std::string &fname =
    gdir + "/axisem3d_synthetics.nc.rank" + mpi::strRank();
    
    // create and open
    mNcFile = std::make_unique<NetCDF_Writer>();
    mNcFile->open(fname, true);
    
    // reset time line
    mFileLineTime = 0;
    
    ///////////////////// define variables /////////////////////
    mNcFile->defModeOn();
    
    //////// time ////////
    mVarID_Time = mNcFile->defineVariable("data_time", {
        {"dim_time", numRecordSteps}}, numerical::dErr);
    
    //////// wave ////////
    // classify elements on na-grid
    std::vector<std::vector<int>> elemsNaGrid;
    elemsNaGrid.assign(naGrid.size(), std::vector<int>());
    for (int ielem = 0; ielem < nelem; ielem++) {
        // four columns: tag, actual na, grid na, index
        int nag = elemNaInfo(ielem, 2);
        int elemTag = elemNaInfo(ielem, 0);
        int nagIndex = naGridIndexDict[nag];
        elemsNaGrid[nagIndex].push_back(elemTag);
    }
    
    // create variable for each grid na
    for (int inag = 0; inag < naGrid.size(); inag++) {
        const std::string &strNag = bstring::toString(naGrid[inag]);
        int varID = mNcFile->defineVariable("data_wave__NaG=" + strNag, {
            {"dim_element__NaG=" + strNag, elemsNaGrid[inag].size()},
            {"dim_na__NaG=" + strNag, naGrid[inag]},
            {"dim_GLL", npnts},
            {"dim_channel", channels.size()},
            {"dim_time", numRecordSteps}
        }, (numerical::Real)numerical::dErr);
        mVarID_Data.push_back(varID);
    }
    
    //////// info ////////
    // channel
    mNcFile->defineVariable("list_channel", {
        {"dim_channel", channels.size()},
        {"dim_channel_str_length", vector_tools::maxLength(channels)}
    }, (char)0);
    
    // element
    for (int inag = 0; inag < naGrid.size(); inag++) {
        const std::string &strNag = bstring::toString(naGrid[inag]);
        mNcFile->defineVariable("list_element__NaG=" + strNag, {
            {"dim_element__NaG=" + strNag, elemsNaGrid[inag].size()}
        }, (int)-1);
    }
    
    // na-grid
    mNcFile->defineVariable("list_na_grid", {
        {"dim_na_grid", naGrid.size()}
    }, (int)-1);
    
    // element-na info
    mNcFile->defineVariable("list_element_na", {
        {"dim_element", nelem}, {"dim_5", 5}
    }, (int)-1);
    
    // element coords
    mNcFile->defineVariable("list_element_coords", {
        {"dim_element", nelem}, {"dim_GLL", npnts}, {"dim_2", 2}
    }, numerical::dErr);
    
    // end defining variables
    mNcFile->defModeOff();
    
    ///////////////////// write info /////////////////////
    // channel
    for (int ich = 0; ich < channels.size(); ich++) {
        mNcFile->writeVariable("list_channel", channels[ich],
                               {ich, 0}, {1, (int)channels[ich].size()});
    }
    
    // element
    for (int inag = 0; inag < naGrid.size(); inag++) {
        const std::string &strNag = bstring::toString(naGrid[inag]);
        mNcFile->writeWholeVariable("list_element__NaG=" + strNag,
                                    elemsNaGrid[inag]);
    }
    
    // na-grid
    mNcFile->writeWholeVariable("list_na_grid", naGrid);
    
    // element-na info
    mNcFile->writeWholeVariable("list_element_na", elemNaInfo5);
    
    // element coords
    mNcFile->writeWholeVariable("list_element_coords", elemCoords);
}

// finalize
void ElementIO_NetCDF::finalize() {
    // close files
    if (mNcFile) {
        mNcFile->close();
        mNcFile.reset(nullptr);
    }
    mFileLineTime = 0;
}

// dump to file
void ElementIO_NetCDF::
dumpToFile(const eigen::DColX &bufferTime,
           const std::vector<eigen::RTensor5> &bufferFields,
           int bufferLine) {
    // no element
    if (!mNcFile) {
        return;
    }
    
    // write time
    mNcFile->writeVariable(mVarID_Time, "data_time", bufferTime,
                           {mFileLineTime}, {bufferLine});
    
    // write data
    for (int inag = 0; inag < bufferFields.size(); inag++) {
        int nelem = (int)bufferFields[inag].dimensions()[0];
        if (nelem == 0) {
            // no element on this rank uses this grid na
            continue;
        }
        int nag = (int)bufferFields[inag].dimensions()[1];
        int npnts = (int)bufferFields[inag].dimensions()[2];
        int nch = (int)bufferFields[inag].dimensions()[3];
        const std::string &strNag = bstring::toString(nag);
        
        // write to file
        if (bufferLine == bufferFields[inag].dimension(4)) {
            // full buffer
            mNcFile->writeVariable(mVarID_Data[inag],
                                   "data_wave__NaG=" + strNag,
                                   bufferFields[inag],
                                   {0, 0, 0, 0, mFileLineTime},
                                   {nelem, nag, npnts, nch, bufferLine});
        } else {
            // must truncate by copy because time is the fastest varying
            // dimension; occuring only at the end of the simulation
            eigen::IArray5 loc = {0, 0, 0, 0, 0};
            eigen::IArray5 len = {nelem, nag, npnts, nch, bufferLine};
            Eigen::internal::set_is_malloc_allowed(true);
            eigen::RTensor5 timeTruncated = bufferFields[inag].slice(loc, len);
            Eigen::internal::set_is_malloc_allowed(false);
            mNcFile->writeVariable(mVarID_Data[inag],
                                   "data_wave__NaG=" + strNag,
                                   timeTruncated,
                                   {0, 0, 0, 0, mFileLineTime},
                                   {nelem, nag, npnts, nch, bufferLine});
        }
    }
    
    // flush
    if (mFlush) {
        mNcFile->flush();
    }
    
    // update record postion in file
    mFileLineTime += bufferLine;
}
