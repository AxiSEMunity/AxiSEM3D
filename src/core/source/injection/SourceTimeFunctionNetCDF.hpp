//
//  SourceTimeFunctionNetCDF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source-time function from NetCDF

#ifndef SourceTimeFunctionNetCDF_hpp
#define SourceTimeFunctionNetCDF_hpp

#include "eigen.hpp"
#include "spectral.hpp"
#include "vector_tools.hpp"

// reader
#include <memory>
class NetCDF_Reader;

template <typename T, int ndim, int nCols = ndim * spectral::nPEM>
class SourceTimeFunctionNetCDF {
public:
    // pattern matrix
    typedef Eigen::Matrix<T, 1, Eigen::Dynamic> TRowX;
    typedef Eigen::Matrix<std::complex<T>, Eigen::Dynamic, nCols> CTMatXND;
    
    // constructor
    SourceTimeFunctionNetCDF(bool alignedToTimeStep, int bufferSize,
                             const std::shared_ptr<NetCDF_Reader> &reader,
                             const std::string &variableName):
    mAlignedToTimeStep(alignedToTimeStep),
    mReader(reader), mVariableName(variableName) {
        // construct derived
        constructDerived(bufferSize);
    }
    
    // get order
    int getPatternNu_1() const {
        return mNu_1;
    }
    
    // get value at a time step
    void getPatternAtTimeStep(int timeStep, double time, CTMatXND &result) {
        if (mAlignedToTimeStep) {
            // check buffer
            while (timeStep > mTimeStepOfPatternLast) {
                loadNextBufferChunk();
            }
            // retrieve by index
            result.topRows(mNu_1) =
            mPatterns[timeStep - mTimeStepOfPattern0];
        } else {
            // check buffer
            while (time > mTimePoints[mTimeStepOfPatternLast -
                                      mTimeStepOfPattern0]) {
                loadNextBufferChunk();
            }
            // retrieve by linear interpolation
            int index0 = -1, index1 = -1;
            double factor0 = 0., factor1 = 0.;
            try {
                vector_tools::linearInterpSorted(mTimePoints, time,
                                                 index0, index1,
                                                 factor0, factor1, 0,
                                                 mTimeStepOfPatternLast -
                                                 mTimeStepOfPattern0 + 1);
            } catch (...) {
                throw std::runtime_error
                ("SourceTimeFunctionNetCDF::getPatternAtTimeStep || "
                 "Target time point is out of range.");
            }
            
            result.topRows(mNu_1) = (mPatterns[index0] * (T)factor0 +
                                     mPatterns[index1] * (T)factor1);
        }
    }
    
    // these private functions avoid including NetCDF_Reader.hpp
private:
    // construct derived
    void constructDerived(int bufferSize);
    
    // load next buffer chunk
    void loadNextBufferChunk();
    
    
private:
    ////////////// alignment //////////////
    // aligned to time steps
    // true: directly retrieve stf using time step as index
    // false: interpolation in time is needed
    const bool mAlignedToTimeStep;
    
    ////////////// file //////////////
    // reader
    const std::shared_ptr<NetCDF_Reader> mReader;
    
    // variable name
    const std::string mVariableName;
    
    // variable id
    int mTimeID = -1;
    int mPatternReID = -1;
    int mPatternImID = -1;
    
    // total time steps in file
    int mTotalTimeStepsInFile = 0;
    
    ////////////// pattern series //////////////
    // buffer of time points
    std::vector<double> mTimePoints;
    
    // buffer of patterns
    std::vector<CTMatXND> mPatterns;
    
    // time step of mTimePoints[0] and mTimePoints[-1]
    int mTimeStepOfPattern0 = -1;
    int mTimeStepOfPatternLast = -1;
    
    // nu_1
    int mNu_1 = 0;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // static buffer for reading
    inline static TRowX sReadRe = TRowX(0);
    inline static TRowX sReadIm = TRowX(0);
};

#endif /* SourceTimeFunctionNetCDF_hpp */
