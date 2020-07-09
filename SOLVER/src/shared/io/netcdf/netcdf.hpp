//
//  netcdf.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/21/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  netcdf header

#ifndef netcdf_hpp
#define netcdf_hpp

extern "C" {
#include <netcdf.h>
}

#ifdef _USE_PARALLEL_NETCDF
extern "C" {
#include <netcdf_par.h>
}
#endif

#include "eigen.hpp"
#include "eigen_tensor.hpp"
#include "bstring.hpp"

namespace netcdf {
    ///////////////////////// datatype /////////////////////////
    // datatype covertion, e.g., int -> NC_INT
    template <typename T>
    inline nc_type typeNC_T() = delete;
    
    // specializations
    template <>
    inline nc_type typeNC_T<char>() {
        return NC_CHAR;
    }
    
    template <>
    inline nc_type typeNC_T<int>() {
        return NC_INT;
    }
    
    template <>
    inline nc_type typeNC_T<long>() {
        return NC_INT64;
    }
    
    template <>
    inline nc_type typeNC_T<float>() {
        return NC_FLOAT;
    }
    
    template <>
    inline nc_type typeNC_T<double>() {
        return NC_DOUBLE;
    }
    
    
    ///////////////////////// allowed containers /////////////////////////
    // vector
    template <typename T>
    inline nc_type typeNC_V(const std::vector<T> &val) {
        return typeNC_T<T>();
    }
    
    // string
    inline nc_type typeNC_V(const std::string &val) {
        return typeNC_T<char>();
    }
    
    // Eigen::Matrix, RowMajor or column vector
    template <typename Derived>
    inline typename std::enable_if<Eigen::DenseBase<Derived>::IsRowMajor ||
    Eigen::DenseBase<Derived>::ColsAtCompileTime == 1, nc_type>::type
    typeNC_V(const Eigen::DenseBase<Derived> &val) {
        return typeNC_T<typename Eigen::DenseBase<Derived>::Scalar>();
    }
    
    // Eigen::Tensor, RowMajor
    template <class Scalar, int R>
    inline nc_type typeNC_V(const Eigen::Tensor<Scalar, R,
                            Eigen::RowMajor> &val) {
        return typeNC_T<Scalar>();
    }
    
    // Eigen::Tensor, ColMajor of rank one
    template <class Scalar>
    inline nc_type typeNC_V(const Eigen::Tensor<Scalar, 1,
                            Eigen::ColMajor> &val) {
        return typeNC_T<Scalar>();
    }
    
    
    ///////////////////////// error /////////////////////////
    // error handler
    inline void error(int retval, const std::string &funcName,
                      const std::string &fname) {
        if (retval != NC_NOERR) {
            throw std::
            runtime_error("netcdf::error || "
                          "Error in NetCDF function: " + funcName + " || "
                          "Error code = " + bstring::toString(retval) + " || "
                          "NetCDF file: " + fname);
        }
    }
    
    
    ///////////////////////// type and id /////////////////////////
    // get id
    inline int varID(int fid, const std::string &vname,
                     const std::string &fname) {
        int varid = -1;
        if (nc_inq_varid(fid, vname.c_str(), &varid) != NC_NOERR) {
            throw std::runtime_error("netcdf::getVariableID || "
                                     "Error finding variable. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + fname);
        }
        return varid;
    }
    
    // get id and check type
    inline int varID_CheckType(int fid, const std::string &vname,
                               const std::string &fname,
                               const nc_type &requiredTypeNC) {
        // get id
        int varid = varID(fid, vname, fname);
        
        // get nc datatype
        nc_type typeInFile;
        error(nc_inq_vartype(fid, varid, &typeInFile),
              "nc_inq_vartype", fname);
        
        // compare datatype
        if (typeInFile != requiredTypeNC) {
            throw std::runtime_error("netcdf::checkGetVariableID || "
                                     "Inconsistent C++ and netcdf datatype. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + fname);
        }
        return varid;
    }
}

#endif /* netcdf_hpp */
