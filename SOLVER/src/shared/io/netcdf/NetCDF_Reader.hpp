//
//  NetCDF_Reader.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/13/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  NetCDF reader

#ifndef NetCDF_Reader_hpp
#define NetCDF_Reader_hpp

#include "numerical.hpp"
#include "netcdf.hpp"

class NetCDF_Reader {
public:
    // constructor
    NetCDF_Reader() {
        // nothing
    }
    
    // constructor
    NetCDF_Reader(const std::string &fname) {
        open(fname);
    }
    
    // destructor
    ~NetCDF_Reader() {
        close();
    }
    
    
    ////////////////// file system //////////////////
    // open
    void open(const std::string &fname);
    
    // open parallel
    void openParallel(const std::string &fname);
    
    // close
    void close();
    
    // is open
    bool isOpen() const {
        return mFileName != "";
    }
    
    
    ///////////// workflow: info -> allocation -> read /////////////
    // get variable ID
    // ID and Container checked, datatype optionally checked
    template <class Container>
    int getVariableID(const std::string &vname, const Container &val,
                      bool checkDateType = true) const {
        if (checkDateType) {
            return netcdf::varID_CheckType(mFileID, vname, mFileName,
                                           netcdf::typeNC_V(val));
        } else {
            return netcdf::varID(mFileID, vname, mFileName);
        }
    }
    
    // get dimensions
    void getVariableDimensions(int varid,
                               std::vector<numerical::Int> &dims) const {
        // get number of dimensions
        int ndims = -1;
        netcdf::error(nc_inq_varndims(mFileID, varid, &ndims),
                      "nc_inq_varndims", mFileName);
        dims.resize(ndims);
        
        // get id of dimensions
        std::vector<int> id_dims(ndims, -1);
        netcdf::error(nc_inq_vardimid(mFileID, varid, id_dims.data()),
                      "nc_inq_vardimid", mFileName);
        
        // get dimensions
        for (int idim = 0; idim < ndims; idim++) {
            netcdf::error(nc_inq_dimlen(mFileID, id_dims[idim],
                                        (size_t *)(&dims[idim])),
                          "nc_inq_dimlen", mFileName);
        }
    }
    
    // read data from variable chunk
    // ID, datatype and container checked by getVariableID
    template <class Container>
    void readVariable(int varid, const std::string &vname, Container &val,
                      const std::vector<numerical::Int> &starts,
                      const std::vector<numerical::Int> &counts) const {
        // sizeof(numerical::Int) = sizeof(size_t), cast directly
        if (nc_get_vara(mFileID, varid, (size_t *)starts.data(),
                        (size_t *)counts.data(), val.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Reader::readVariable || "
                                     "Error reading variable from file. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    
    ///////////// read with type check /////////////
    // Eigen::Matrix
    template <class Derived>
    void readMatrix(const std::string &vname,
                    Eigen::DenseBase<Derived> &mat) const {
        // access variable with datatype check
        int varid = getVariableID(vname, mat, true);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() == 1) {
            // vector to matrix
            dims.push_back(1);
        }
        if (dims.size() != 2) {
            throw std::runtime_error("NetCDF_Reader::readMatrix || "
                                     "NC variable is not a matrix. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // must check compile-time size before allocation
        int rowsFix = Eigen::DenseBase<Derived>::RowsAtCompileTime;
        int colsFix = Eigen::DenseBase<Derived>::ColsAtCompileTime;
        if ((rowsFix == Eigen::Dynamic || rowsFix == dims[0]) &&
            (colsFix == Eigen::Dynamic || colsFix == dims[1])) {
            // allocate and read
            mat.derived().resize(dims[0], dims[1]);
            netcdf::error(nc_get_var(mFileID, varid, mat.derived().data()),
                          "nc_get_var", mFileName);
        } else {
            throw std::runtime_error("NetCDF_Reader::readMatrix || "
                                     "Inconsistent variable dimensions. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    // std::vector
    template <class T>
    void readVector(const std::string &vname, std::vector<T> &vec) const {
        // access variable with datatype check
        int varid = getVariableID(vname, vec, true);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() != 1) {
            throw std::runtime_error("NetCDF_Reader::readVector || "
                                     "NC variable is not a vector. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // allocate and read
        vec.resize(dims[0]);
        netcdf::error(nc_get_var(mFileID, varid, vec.data()),
                      "nc_get_var", mFileName);
    }
    
    // Eigen::Tensor
    template <class T, int R>
    void readTensor(const std::string &vname,
                    Eigen::Tensor<T, R, Eigen::RowMajor> &tensor) const {
        // access variable with datatype check
        int varid = getVariableID(vname, tensor, true);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() != R) {
            throw std::runtime_error("NetCDF_Reader::readTensor || "
                                     "Inconsistent tensor ranks. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // allocate and read
        std::array<numerical::Int, R> dimsArray;
        std::copy(dims.begin(), dims.end(), dimsArray.begin());
        tensor.resize(dimsArray);
        netcdf::error(nc_get_var(mFileID, varid, tensor.data()),
                      "nc_get_var", mFileName);
    }
    
    // read string
    void readString(const std::string &vname,
                    std::vector<std::string> &dest) const;
    
    
    ///////////// read with enforced type of double /////////////
    // Eigen::Matrix
    template <class Derived,
    typename T = typename Eigen::DenseBase<Derived>::Scalar>
    typename std::enable_if<std::is_floating_point<T>::value &&
    sizeof(T) == sizeof(double), void>::type
    readMatrixDouble(const std::string &vname,
                     Eigen::DenseBase<Derived> &mat) const {
        // access variable without datatype check
        int varid = getVariableID(vname, mat, false);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() == 1) {
            // vector to matrix
            dims.push_back(1);
        }
        if (dims.size() != 2) {
            throw std::runtime_error("NetCDF_Reader::readMatrixDouble || "
                                     "NC variable is not a matrix. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // must check compile-time size before allocation
        int rowsFix = Eigen::DenseBase<Derived>::RowsAtCompileTime;
        int colsFix = Eigen::DenseBase<Derived>::ColsAtCompileTime;
        if ((rowsFix == Eigen::Dynamic || rowsFix == dims[0]) &&
            (colsFix == Eigen::Dynamic || colsFix == dims[1])) {
            // allocate and read
            mat.derived().resize(dims[0], dims[1]);
            netcdf::error(nc_get_var_double(mFileID, varid, mat.derived().data()),
                          "nc_get_var_double", mFileName);
        } else {
            throw std::runtime_error("NetCDF_Reader::readMatrixDouble || "
                                     "Inconsistent variable dimensions. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
    }
    
    // std::vector
    void readVectorDouble(const std::string &vname,
                          std::vector<double> &vec) const {
        // access variable without datatype check
        int varid = getVariableID(vname, vec, false);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() != 1) {
            throw std::runtime_error("NetCDF_Reader::readVectorDouble || "
                                     "NC variable is not a vector. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // allocate and read
        vec.resize(dims[0]);
        netcdf::error(nc_get_var_double(mFileID, varid, vec.data()),
                      "nc_get_var_double", mFileName);
    }
    
    // Eigen::Tensor
    template <int R>
    void readTensorDouble(const std::string &vname,
                          Eigen::Tensor<double, R,
                          Eigen::RowMajor> &tensor) const {
        // access variable without datatype check
        int varid = getVariableID(vname, tensor, false);
        
        // get dimensions
        std::vector<numerical::Int> dims;
        getVariableDimensions(varid, dims);
        if (dims.size() != R) {
            throw std::runtime_error("NetCDF_Reader::readTensorDouble || "
                                     "Inconsistent tensor ranks. || "
                                     "Variable name: " + vname + " || "
                                     "NetCDF file: " + mFileName);
        }
        
        // allocate and read
        std::array<numerical::Int, R> dimsArray;
        std::copy(dims.begin(), dims.end(), dimsArray.begin());
        tensor.resize(dimsArray);
        netcdf::error(nc_get_var_double(mFileID, varid, tensor.data()),
                      "nc_get_var_double", mFileName);
    }
    
    
    ////////////////// data //////////////////
private:
    // file ID
    int mFileID = -1;
    
    // file name
    std::string mFileName = "";
};

#endif /* NetCDF_Reader_hpp */
