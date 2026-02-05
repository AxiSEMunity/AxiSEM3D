//
//  ElementIO_ParNetCDF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  parallel NetCDF IO for element output

#include "ElementIO_ParNetCDF.hpp"
#include "io.hpp"
#include "mpi.hpp"
#include "vector_tools.hpp"

// initialize
void ElementIO_ParNetCDF::initialize(const std::string &groupName,
                                     int numRecordSteps,
                                     const std::vector<std::string> &channels,
                                     int npnts, const std::vector<int> &naGrid,
                                     const eigen::IMatX4_RM &elemNaInfoL,
                                     const eigen::DMatXX_RM &elemCoordsL) {
    // finalize
    finalize();
    
    // base
    ElementIO::initialize(groupName, numRecordSteps, channels,
                          npnts, naGrid, elemNaInfoL, elemCoordsL);
    
    // return only if no element in group
    if (mNumElementsGlobal == 0) {
        return;
    }
    
    // filename
    const std::string &gdir = io::gOutputDirectory + "/elements/" + groupName;
    const std::string &fname = gdir + "/axisem3d_synthetics.nc";
    
    // gather elements
    std::vector<eigen::IMatX4_RM> elemNaInfoRanks;
    std::vector<eigen::DMatXX_RM> elemCoordsRanks;
    mpi::gatherEigen(elemNaInfoL, elemNaInfoRanks, mRankWithMaxNumElements);
    mpi::gatherEigen(elemCoordsL, elemCoordsRanks, mRankWithMaxNumElements);
    
    // first element index on na-grid
    std::vector<std::vector<int>> firstElemNaGridRanks;
    
    // file pointer
    mNcFile = std::make_unique<NetCDF_Writer>();
    
    /////////////////////// only on max rank ///////////////////////
    if (mpi::rank() == mRankWithMaxNumElements) {
        // flatten
        int ncrd = npnts * 2;
        eigen::IMatX5_RM elemNaInfoAll(mNumElementsGlobal, 5);
        eigen::DMatXX_RM elemCoordsAll(mNumElementsGlobal, ncrd);
        int row = 0;
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            int nelocal = (int)elemNaInfoRanks[irank].rows();
            elemNaInfoAll.block(row, 0, nelocal, 4) = elemNaInfoRanks[irank];
            elemCoordsAll.block(row, 0, nelocal, ncrd) = elemCoordsRanks[irank];
            row += nelocal;
        }
        
        // form dict of naGrid for fast search
        std::map<int, int> naGridIndexDict;
        for (int inag = 0; inag < naGrid.size(); inag++) {
            naGridIndexDict[naGrid[inag]] = inag;
        }
        
        // first element index on na-grid
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
        
        // the fifth column: element index in data after concat
        row = 0;
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            int nelocal = (int)elemNaInfoRanks[irank].rows();
            for (int ielem = 0; ielem < nelocal; ielem++) {
                int nag = elemNaInfoRanks[irank](ielem, 2);
                int nagIndex = naGridIndexDict[nag];
                elemNaInfoAll(row, 4) = elemNaInfoAll(row, 3) +
                firstElemNaGridRanks[irank][nagIndex];
                row++;
            }
        }
        
        
        //////////////////////////////////////////////////////////////
        //////////////////////////// file ////////////////////////////
        //////////////////////////////////////////////////////////////
        
        // to reuse code in ElementIO_NetCDF.cpp
        int nelem = mNumElementsGlobal;
        
        // create and open
        mNcFile->open(fname, true);
        
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
            int nag = elemNaInfoAll(ielem, 2);
            int elemTag = elemNaInfoAll(ielem, 0);
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
        mNcFile->writeWholeVariable("list_element_na", elemNaInfoAll);
        
        // element coords
        mNcFile->writeWholeVariable("list_element_coords", elemCoordsAll);
    }
    
    // bcast variable IDs
    // variable ID should not change upon re-opening
    mpi::bcast(mVarID_Time, mRankWithMaxNumElements);
    mpi::bcast(mVarID_Data, mRankWithMaxNumElements);
    
    // open parallel
    mpi::barrier();
    mNcFile->openParallel(fname);
    
    // reset time line
    mFileLineTime = 0;
    
    // scatter first element index on na-grid
    mpi::scatter(firstElemNaGridRanks, mFirstElemNaGrid,
                 mRankWithMaxNumElements);
}

// finalize
void ElementIO_ParNetCDF::finalize() {
    // close files
    if (mNcFile) {
        mNcFile->close();
        mNcFile.reset(nullptr);
    }
    mFileLineTime = 0;
}

// dump to file
void ElementIO_ParNetCDF::
dumpToFile(const eigen::DColX &bufferTime,
           const std::vector<eigen::RTensor5> &bufferFields,
           int bufferLine) {
    // write time
    if (mpi::rank() == mRankWithMaxNumElements) {
        mNcFile->writeVariable(mVarID_Time, "data_time", bufferTime,
                               {mFileLineTime}, {bufferLine});
    }
    
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
            mNcFile->writeVariable
            (mVarID_Data[inag], "data_wave__NaG=" + strNag, bufferFields[inag],
             {mFirstElemNaGrid[inag], 0, 0, 0, mFileLineTime},
             {nelem, nag, npnts, nch, bufferLine});
        } else {
            // must truncate by copy because time is the fastest varying
            // dimension; occuring only at the end of the simulation
            eigen::IArray5 loc = {0, 0, 0, 0, 0};
            eigen::IArray5 len = {nelem, nag, npnts, nch, bufferLine};
            Eigen::internal::set_is_malloc_allowed(true);
            eigen::RTensor5 timeTruncated = bufferFields[inag].slice(loc, len);
            Eigen::internal::set_is_malloc_allowed(false);
            mNcFile->writeVariable
            (mVarID_Data[inag], "data_wave__NaG=" + strNag, timeTruncated,
             {mFirstElemNaGrid[inag], 0, 0, 0, mFileLineTime},
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
