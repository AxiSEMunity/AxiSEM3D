//
//  NetCDF_Reader.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  NetCDF reader

#include "NetCDF_Reader.hpp"

////////////////// file system //////////////////
// open
void NetCDF_Reader::open(const std::string &fname) {
    close();
    if (nc_open(fname.c_str(), NC_NETCDF4, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::open || "
                                 "Error opening NetCDF file: || " + fname);
    }
    // this must come after nc_open
    mFileName = fname;
}

// open parallel
void NetCDF_Reader::openParallel(const std::string &fname) {
#ifdef _USE_PARALLEL_NETCDF
    close();
    if (nc_open_par(fname.c_str(), NC_MPIIO | NC_NETCDF4, MPI_COMM_WORLD,
                    MPI_INFO_NULL, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::openParallel || "
                                 "Error opening NetCDF file: || " + fname);
    }
    mFileName = fname;
#else
    throw std::runtime_error("NetCDF_Reader::openParallel || "
                             "CMakeLists.txt has disabled Parallel NetCDF.");
#endif
}

// close
void NetCDF_Reader::close() {
    if (isOpen()) {
        netcdf::error(nc_close(mFileID), "nc_close", mFileName);
        mFileID = -1;
        mFileName = "";
    }
}


////////////////// read + allocation //////////////////
// read string
void NetCDF_Reader::readString(const std::string &vname,
                               std::vector<std::string> &dest) const {
    // access variable
    int varid = getVariableID(vname, std::string(""), true);;
    
    // get ndims
    std::vector<numerical::Int> dims;
    getVariableDimensions(varid, dims);
    if (dims.size() != 2) {
        throw std::runtime_error("NetCDF_Reader::readString || "
                                 "Number of dimensions must be 2. || "
                                 "Variable name: " + vname + " || "
                                 "NetCDF file: " + mFileName);
    }
    
    // long buffer
    size_t numString = dims[0];
    size_t lenString = dims[1];
    std::vector<char> cstr;
    cstr.resize(numString * lenString);
    
    // read and split
    netcdf::error(nc_get_var_text(mFileID, varid, cstr.data()),
                  "nc_get_var_text", mFileName);
    dest.resize(numString);
    for (int istr = 0; istr < numString; istr++) {
        dest[istr] = std::string(&cstr[istr * lenString]);
    }
}
