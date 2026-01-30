//
//  StructuredGrid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  data on a structured grid
//  class parameters:
//  - RANK: number of dimensions
//  - DTYPE: datatype of data

#ifndef StructuredGrid_hpp
#define StructuredGrid_hpp

#include "NetCDF_Reader.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "eigen_tools.hpp"
#include "vector_tools.hpp"
#include <set>

template <int D, typename T>
class StructuredGrid {
public:
    // constructor
    StructuredGrid() = delete;
    
    // constructor
    StructuredGrid(const std::string &ncFile,
                   const std::array<std::string, D> &coordVarNames,
                   const std::vector<std::pair<std::string, double>> &dataInfo,
                   const std::array<int, D> &shuffleData) {
        // verify shuffle of data
        std::array<int, D> shuffleSorted = shuffleData;
        std::sort(shuffleSorted.begin(), shuffleSorted.end());
        for (int idim = 0; idim < D; idim++) {
            if (shuffleSorted[idim] != idim) {
                throw std::runtime_error("StructuredGrid::StructuredGrid || "
                                         "Invalid data ranks to shuffle.");
            }
        }
        
        // find unique data variable names
        std::set<std::string> uniqueDataVarNames;
        for (auto it = dataInfo.begin(); it != dataInfo.end(); ++it) {
            uniqueDataVarNames.insert(it->first);
        }
        
        // data rank and factor
        for (auto it = dataInfo.begin(); it != dataInfo.end(); ++it) {
            int loc = (int)std::distance(uniqueDataVarNames.begin(),
                                         uniqueDataVarNames.find(it->first));
            mDataRankFactor.push_back({loc, it->second});
        }
        
        // read data
        timer::gPreloopTimer.begin("Reading grid data");
        timer::gPreloopTimer.message("data file: " + io::popInputDir(ncFile));
        std::vector<std::vector<double>> gridCoords;
        if (mpi::root()) {
            // open
            NetCDF_Reader reader(io::popInputDir(ncFile));
            
            // read coords
            for (const std::string &cvName: coordVarNames) {
                std::vector<double> crd;
                reader.readVectorDouble(cvName, crd);
                gridCoords.push_back(crd);
            }
            
            // coord dimensions
            Eigen::array<Eigen::DenseIndex, D> dimsCrd;
            for (int idim = 0; idim < D; idim++) {
                dimsCrd[idim] = gridCoords[idim].size();
            }
            
            // allocate data
            Eigen::array<Eigen::DenseIndex, 1 + D> dimsData;
            dimsData[0] = uniqueDataVarNames.size(); // # data rank comes first
            std::copy(dimsCrd.begin(), dimsCrd.end(), dimsData.begin() + 1);
            mGridData = Eigen::Tensor<T, 1 + D, Eigen::RowMajor>(dimsData);
            
            // shape
            Eigen::array<Eigen::DenseIndex, 1 + D> start, count;
            start.fill(0);
            count = mGridData.dimensions();
            count[0] = 1;
            
            // read data
            for (const std::string &vname: uniqueDataVarNames) {
                // read to double
                Eigen::Tensor<double, D, Eigen::RowMajor> dGridData;
                reader.readTensorDouble(vname, dGridData);
                // round
                if (std::is_integral<T>::value) {
                    dGridData = dGridData.round();
                }
                // copy data
                mGridData.slice(start, count).reshape(dimsCrd)
                = dGridData.shuffle(shuffleData).template cast<T>();
                // next
                start[0]++;
            }
            
            // check
            for (int idim = 0; idim < D; idim++) {
                // size
                if (gridCoords[idim].size() != mGridData.dimension(idim + 1)) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Incompatible grid sizes in coordinates and data. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncFile);
                }
                
                // size
                if (gridCoords[idim].size() < 2) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Too few grid points; at least two are required. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncFile);
                }
                
                // sorted
                if (!vector_tools::isSortedUnique(gridCoords[idim])) {
                    throw std::runtime_error
                    ("StructuredGrid::StructuredGrid || "
                     "Coordinates are not ascendingly sorted. || "
                     "Dimension: " + coordVarNames[idim] + " || "
                     "NetCDF file: " + ncFile);
                }
            }
        }
        timer::gPreloopTimer.ended("Reading grid data");
        
        // broadcast
        timer::gPreloopTimer.begin("Broadcasting grid data");
        mpi::bcast(gridCoords);
        // cast to array
        for (int idim = 0; idim < D; idim++) {
            mGridCoords[idim] = gridCoords[idim];
        }
        
        // below we broadcast the tensor manually
        // add mpi::bcastEigenTensor if needed in the future
        Eigen::array<Eigen::DenseIndex, 1 + D> dims;
        dims[0] = uniqueDataVarNames.size();
        for (int idim = 0; idim < D; idim++) {
            dims[idim + 1] = mGridCoords[idim].size();
        }
        // allocate
        if (mpi::rank() != 0) {
            mGridData.resize(dims);
        }
        // broadcast tensor data
        mpi::bcast(mGridData.data(), (int)mGridData.size());
        
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfoTensor(mGridData, "unique grid data"));
        timer::gPreloopTimer.ended("Broadcasting grid data");
    }
    
    // in scope
    template <typename ArrayCS>
    bool inScope(const ArrayCS &crdTarget) const {
        // linear interp 0-1 points
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        for (int idim = 0; idim < D; idim++) {
            try {
                vector_tools::
                linearInterpSorted(mGridCoords[idim], crdTarget[idim],
                                   index0, index1, factor0, factor1);
            } catch (...) {
                // location out of range
                return false;
            }
        }
        return true;
    }
    
    // compute
    typedef Eigen::Matrix<T, 1, Eigen::Dynamic> TRowV;
    template <typename ArrayCS>
    TRowV compute(const ArrayCS &crdTarget,
                  const TRowV &valOutOfRange, int &dimOutOfRange) const {
        // linear interp 0-1 points
        static std::array<std::array<int, D>, 2> index01;
        static std::array<std::array<double, D>, 2> factor01;
        for (int idim = 0; idim < D; idim++) {
            try {
                vector_tools::
                linearInterpSorted(mGridCoords[idim], crdTarget[idim],
                                   index01[0][idim], index01[1][idim],
                                   factor01[0][idim], factor01[1][idim]);
            } catch (...) {
                // location out of range
                dimOutOfRange = idim;
                return valOutOfRange;
            }
        }
        
        // shapes
        Eigen::array<Eigen::DenseIndex, 1> shapeV = {numUniqueData()};
        Eigen::array<Eigen::DenseIndex, 1 + D> start, count;
        start.fill(0);
        count.fill(1);
        count[0] = numUniqueData();
        
        // loop over 2^D points
        typedef Eigen::Matrix<T, 1, Eigen::Dynamic> DRowV;
        DRowV uData = DRowV::Zero(numUniqueData());
        for (int deci = 0; deci < pow(2, D); deci++) {
            int left = deci;
            T factorCombine = (T)1.;
            for (int idim = 0; idim < D; idim++) {
                start[idim + 1] = index01[left % 2][idim];
                factorCombine *= factor01[left % 2][idim];
                left /= 2;
            }
            const Eigen::Tensor<T, 1, Eigen::RowMajor> &row =
            factorCombine * mGridData.slice(start, count).reshape(shapeV);
            uData += Eigen::Map<const DRowV>(row.data(), numUniqueData());
        }
        
        // map to non-unique and apply factor
        DRowV data = DRowV::Zero(mDataRankFactor.size());
        for (int ivar = 0; ivar < numData(); ivar++) {
            int index = mDataRankFactor[ivar].first;
            T factor = mDataRankFactor[ivar].second;
            data(ivar) = uData(index) * factor;
        }
        
        // round
        if (std::is_integral<T>::value) {
            data = data.array().round();
        }
        
        // return
        dimOutOfRange = -1;
        return data.template cast<T>();
    }
    
    // compute
    template <typename ArrayCS>
    TRowV compute(const ArrayCS &crdTarget,
                  const TRowV &valOutOfRange) const {
        static int dimOutOfRange = 0;
        return compute(crdTarget, valOutOfRange, dimOutOfRange);
    }
    
    // compute single
    template <typename ArrayCS>
    T compute(const ArrayCS &crdTarget, T valOutOfRange,
              int varIndex = 0) const {
        TRowV oneVal = TRowV::Zero(numData());
        oneVal(varIndex) = valOutOfRange;
        return compute(crdTarget, oneVal)(varIndex);
    }
    
    // get grid coords
    const std::array<std::vector<double>, D> &getGridCoords() const {
        return mGridCoords;
    }
    
    // get grid coords
    std::array<std::vector<double>, D> &getGridCoords() {
        return mGridCoords;
    }
    
    // get grid data
    const Eigen::Tensor<T, 1 + D, Eigen::RowMajor> &getGridData() const {
        return mGridData;
    }
    
    // get data range
    eigen::DMatXX getDataRange() const {
        // min/max of unique data
        eigen::DMatXX uniqueMinMax(numUniqueData(), 2);
        Eigen::array<Eigen::DenseIndex, 1 + D> start, count;
        start.fill(0);
        count = mGridData.dimensions();
        count[0] = 1;
        typedef Eigen::Tensor<double, 0, Eigen::RowMajor> Tensor0;
        for (; start[0] < numUniqueData(); start[0]++) {
            const auto &sliced = mGridData.slice(start, count);
            uniqueMinMax(start[0], 0) = ((Tensor0)sliced.minimum())(0);
            uniqueMinMax(start[0], 1) = ((Tensor0)sliced.maximum())(0);
        }
        
        // min/max of data
        eigen::DMatXX dataMinMax(numData(), 2);
        for (int ivar = 0; ivar < numData(); ivar++) {
            int index = mDataRankFactor[ivar].first;
            double factor = mDataRankFactor[ivar].second;
            dataMinMax(ivar, 0) = uniqueMinMax(index, 0) * factor;
            dataMinMax(ivar, 1) = uniqueMinMax(index, 1) * factor;
        }
        return dataMinMax;
    }
    
    // get number of data
    inline int numData() const {
        return (int)mDataRankFactor.size();
    }
    
    // get number of unique data
    inline int numUniqueData() const {
        return (int)mGridData.dimension(0);
    }
    
private:
    // grid coords
    std::array<std::vector<double>, D> mGridCoords;
    // grid data
    Eigen::Tensor<T, 1 + D, Eigen::RowMajor> mGridData;
    // data rank and factor
    // NOTE: to avoid store duplicated data using the same nc variable
    //       but different factors
    std::vector<std::pair<int, double>> mDataRankFactor;
};

#endif /* StructuredGrid_hpp */
