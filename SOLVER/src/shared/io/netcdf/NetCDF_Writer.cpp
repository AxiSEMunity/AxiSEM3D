//
//  NetCDF_Writer.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  NetCDF writer

#include "NetCDF_Writer.hpp"

////////////////// file system //////////////////
// open
void NetCDF_Writer::open(const std::string &fname, bool overwrite) {
    close();
    if (overwrite) {
        if (nc_create(fname.c_str(), NC_NETCDF4, &mFileID) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::open || "
                                     "Error creating NetCDF file: || " + fname);
        }
        defModeOff();
    } else {
        if (nc_open(fname.c_str(),  NC_WRITE | NC_NETCDF4, &mFileID)
            != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::open || "
                                     "Error opening NetCDF file: || " + fname);
        }
    }
    mPWD = mFileID;
    mFileName = fname;
}

// open parallel
void NetCDF_Writer::openParallel(const std::string &fname) {
#ifdef _USE_PARALLEL_NETCDF
    close();
    if (nc_open_par(fname.c_str(),  NC_MPIIO | NC_WRITE | NC_NETCDF4,
                    MPI_COMM_WORLD, MPI_INFO_NULL, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Writer::openParallel || "
                                 "Error opening NetCDF file: || " + fname);
    }
    mPWD = mFileID;
    mFileName = fname;
#else
    throw std::runtime_error("NetCDF_Writer::openParallel || "
                             "CMakeLists.txt has disabled Parallel NetCDF.");
#endif
}

// close
void NetCDF_Writer::close() {
    if (isOpen()) {
        netcdf::error(nc_close(mFileID), "nc_close", mFileName);
        mPWD = mFileID = -1;
        mFileName = "";
    }
}


////////////////// group //////////////////
// create group
void NetCDF_Writer::createGroup(const std::string &gname) const {
    int grpid = -1;
    netcdf::error(nc_def_grp(mPWD, gname.c_str(), &grpid),
                  "nc_def_grp", mFileName);
}

// go to group
void NetCDF_Writer::goToGroup(const std::string &gname) {
    int grpid = -1;
    netcdf::error(nc_inq_grp_ncid(mPWD, gname.c_str(), &grpid),
                  "nc_inq_grp_ncid", mFileName);
    mPWD = grpid;
}


////////////////// specialization //////////////////
// add string attribute
template <>
void NetCDF_Writer::addAttribute<std::string>(const std::string &attname,
                                              const std::string &attvalue,
                                              const std::string &vname) const {
    if (vname == "") {
        // file attribute
        netcdf::error(nc_put_att_text(mFileID, NC_GLOBAL, attname.c_str(),
                                      attvalue.length(), attvalue.c_str()),
                      "nc_put_att_text", mFileName);
        
        
    } else {
        // variable attribute
        int varid = netcdf::varID(mPWD, vname, mFileName);
        netcdf::error(nc_put_att_text(mPWD, varid, attname.c_str(),
                                      attvalue.length(), attvalue.c_str()),
                      "nc_put_att_text", mFileName);
    }
}
