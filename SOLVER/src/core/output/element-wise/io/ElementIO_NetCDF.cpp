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
                                  const std::vector<int> &ipnts,
                                  const std::vector<int> &naGrid,
                                  const eigen::IMatX4_RM &elemNaInfo,
                                  const eigen::DMatXX_RM &elemCoords) {
    // finalize
    finalize();
    
    // base
    ElementIO::initialize(groupName, numRecordSteps, channels,
                          ipnts, naGrid, elemNaInfo, elemCoords);
    
    // nothing locally without elements
    int nelem = (int)elemNaInfo.rows();
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
    // form dict of naGrid for fast search
    std::map<int, int> naGridIndexDict;
    for (int inag = 0; inag < naGrid.size(); inag++) {
        naGridIndexDict[naGrid[inag]] = inag;
    }
    
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
            {"dim_GLL", ipnts.size()},
            {"dim_channel", channels.size()},
            {"dim_time", numRecordSteps}
        }, (numerical::Real)0.); // must fill with zeros
        mVarID_Data.push_back(varID);
    }
    
    //////// info ////////
    // channel
    mNcFile->defineVariable("list_channel", {
        {"dim_channel", channels.size()},
        {"dim_channel_str_length", vector_tools::maxLength(channels)}
    }, (char)0);
    
    // GLL
    mNcFile->defineVariable("list_GLL", {
        {"dim_GLL", ipnts.size()}
    }, (int)-1);
    
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
    }, (int)0);
    
    // element-na info
    mNcFile->defineVariable("list_element_na", {
        {"dim_element", nelem}, {"dim_4", 4}
    }, (int)0);
    
    // element coords
    mNcFile->defineVariable("list_element_coords", {
        {"dim_element", nelem}, {"dim_GLL", ipnts.size()}, {"dim_2", 2}
    }, (double)0.);
    
    // end defining variables
    mNcFile->defModeOff();
    
    ///////////////////// write info /////////////////////
    // channel
    for (int ich = 0; ich < channels.size(); ich++) {
        mNcFile->writeVariable("list_channel", channels[ich],
                               {ich, 0}, {1, (int)channels[ich].size()});
    }
    
    // GLL
    mNcFile->writeWholeVariable("list_GLL", ipnts);
    
    // element
    for (int inag = 0; inag < naGrid.size(); inag++) {
        const std::string &strNag = bstring::toString(naGrid[inag]);
        mNcFile->writeWholeVariable("list_element__NaG=" + strNag,
                                    elemsNaGrid[inag]);
    }
    
    // na-grid
    mNcFile->writeWholeVariable("list_na_grid", naGrid);
    
    // element-na info
    mNcFile->writeWholeVariable("list_element_na", elemNaInfo);
    
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
    
    // write datas
    for (int inag = 0; inag < bufferFields.size(); inag++) {
        int nelem = (int)bufferFields[inag].dimensions()[0];
        int nag = (int)bufferFields[inag].dimensions()[1];
        int npnts = (int)bufferFields[inag].dimensions()[2];
        int nch = (int)bufferFields[inag].dimensions()[3];
        const std::string &strNag = bstring::toString(nag);
        mNcFile->writeVariable(mVarID_Data[inag], "data_wave__NaG=" + strNag,
                               bufferFields[inag],
                               {0, 0, 0, 0, mFileLineTime},
                               {nelem, nag, npnts, nch, bufferLine});
    }
    
    // flush?
    // mNcFile->flush();
    
    // update record postion in file
    mFileLineTime += bufferLine;
}
