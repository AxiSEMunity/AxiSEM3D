//
//  NetCDF_Writer.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  NetCDF writer

#ifndef NetCDF_Writer_hpp
#define NetCDF_Writer_hpp

#include "numerical.hpp"
#include "netcdf.hpp"

class NetCDF_Writer {
public:
    // constructor
    NetCDF_Writer() {
        // nothing
    }
    
    // constructor
    NetCDF_Writer(const std::string &fname, bool overwrite) {
        open(fname, overwrite);
    }
    
    // destructor
    ~NetCDF_Writer() {
        close();
    }
    
    
    ////////////////// file system //////////////////
    // open
    void open(const std::string &fname, bool overwrite);
    
    // open parallel
    void openParallel(const std::string &fname);
    
    // close
    void close();
    
    // is open
    bool isOpen() const {
        return mFileName != "";
    }
    
    
    ////////////////// define variable //////////////////
    // turn on def mode
    void defModeOn() const {
        netcdf::error(nc_redef(mFileID), "nc_redef", mFileName);
    }
    
    // turn off def mode
    void defModeOff() const {
        netcdf::error(nc_enddef(mFileID), "nc_enddef", mFileName);
    }
    
    // define variable
    template <class scalar>
    int defineVariable(const std::string &vname,
                       const std::vector<
                       std::pair<std::string, numerical::Int>> &dims,
                       const scalar &fill, bool collective = true) const {
        // define dimensions
        std::vector<int> dimids(dims.size(), -1);
        for (int idim = 0; idim < dims.size(); idim++) {
            // find dimension ID
            if (nc_inq_dimid(mFileID, dims[idim].first.c_str(),
                             &(dimids[idim])) != NC_NOERR) {
                // if not found, create it
                netcdf::error(nc_def_dim(mFileID, dims[idim].first.c_str(),
                                         (size_t)dims[idim].second,
                                         &dimids[idim]),
                              "nc_def_dim", mFileName);
            }
        }
        
        // define the variable
        int varid = -1;
        netcdf::error(nc_def_var(mPWD, vname.c_str(),
                                 netcdf::typeNC_T<scalar>(),
                                 (int)dims.size(), dimids.data(), &varid),
                      "nc_def_var", mFileName);
        
        // fill with constant
        netcdf::error(nc_def_var_fill(mPWD, varid, NC_FILL, &fill),
                      "nc_def_var_fill", mFileName);
        
        // this call to nc_var_par_access is not valid
        // // collective
        // #ifdef _USE_PARALLEL_NETCDF
        // if (collective) {
        //     netcdf::error(nc_var_par_access(mPWD, varid, NC_COLLECTIVE),
        //                   "nc_var_par_access", mFileName);
        // } else {
        //     netcdf::error(nc_var_par_access(mPWD, varid, NC_INDEPENDENT),
        //                   "nc_var_par_access", mFileName);
        // }
        // #endif
        
        // return id
        return varid;
    }
    
    
    ////////////////// write variable //////////////////
    // write data to variable chunk
    // this is UNSAFE without checking
    // * ID existence
    // * datatype consistency
    // * container structure (e.g., row/col-major)
    template <class Container>
    void writeVariable(int varid, const std::string &vname,
                       const Container &val,
                       const std::vector<numerical::Int> &starts,
                       const std::vector<numerical::Int> &counts) const {
        // sizeof(numerical::Int) = sizeof(size_t), cast directly
        if (nc_put_vara(mPWD, varid, (size_t *)starts.data(),
                        (size_t *)counts.data(), val.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::writeVariable || "
                                     "Error writing variable to file. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    // write data to variable chunk
    // this is SAFE with everything checked
    template <class Container>
    void writeVariable(const std::string &vname, const Container &val,
                       const std::vector<numerical::Int> &starts,
                       const std::vector<numerical::Int> &counts) const {
        int varid = netcdf::varID_CheckType(mPWD, vname, mFileName,
                                            netcdf::typeNC_V(val));
        writeVariable(varid, vname, val, starts, counts);
    }
    
    // write data to the whole variable
    // this is SAFE with everything checked except size match
    template <class Container>
    void writeWholeVariable(const std::string &vname,
                            const Container &val) const {
        int varid = netcdf::varID_CheckType(mPWD, vname, mFileName,
                                            netcdf::typeNC_V(val));
        if (nc_put_var(mPWD, varid, val.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::writeWholeVariable || "
                                     "Error writing variable to file. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    // flush
    void flush() const {
        netcdf::error(nc_sync(mFileID), "nc_sync", mFileName);
    }
    
    
    ////////////////// group //////////////////
    // create group
    void createGroup(const std::string &gname) const;
    
    // go to group
    void goToGroup(const std::string &gname);
    
    // go to root
    void goToFileRoot() {
        mPWD = mFileID;
    }
    
    
    ////////////////// attribute //////////////////
    // add attribute
    template <class scalar>
    void addAttribute(const std::string &attname, const scalar &attvalue,
                      const std::string &vname = "") const {
        if (vname == "") {
            // file attribute
            netcdf::error(nc_put_att(mFileID, NC_GLOBAL, attname.c_str(),
                                     netcdf::typeNC_T<scalar>(), 1, &attvalue),
                          "nc_put_att", mFileName);
        } else {
            // variable attribute
            int varid = netcdf::varID(mPWD, vname, mFileName);
            netcdf::error(nc_put_att(mPWD, varid, attname.c_str(),
                                     netcdf::typeNC_T<scalar>(), 1, &attvalue),
                          "nc_put_att", mFileName);
        }
    }
    
    
    ////////////////// data //////////////////
private:
    // file ID
    int mFileID = -1;
    
    // file name
    std::string mFileName = "";
    
    // current target ID
    int mPWD = -1;
};


////////////////// specialization //////////////////
// add string attribute
template <>
void NetCDF_Writer::addAttribute<std::string>(const std::string &attname,
                                              const std::string &attvalue,
                                              const std::string &vname) const;

#endif /* NetCDF_Writer_hpp */
