//
//  SourceTimeFunctionNetCDF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/21/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source-time function from NetCDF

#include "SourceTimeFunctionNetCDF.hpp"
#include "NetCDF_Reader.hpp"

// construct derived
template <typename T, int ndim, int nCols> void
SourceTimeFunctionNetCDF<T, ndim, nCols>::constructDerived(int bufferSize) {
    // variable id
    mTimeID = mReader->getVariableID("time_points", mTimePoints);
    mPatternReID = mReader->getVariableID(mVariableName + "_RE", sReadRe);
    mPatternImID = mReader->getVariableID(mVariableName + "_IM", sReadIm);
    
    // dimensions
    std::vector<numerical::Int> dims;
    mReader->getVariableDimensions(mPatternReID, dims);
    if (dims.size() != 2) {
        throw std::runtime_error("SourceTimeFunctionNetCDF::"
                                 "SourceTimeFunctionNetCDF || "
                                 "STF data must have 2 dimensions.");
    }
    mTotalTimeStepsInFile = (int)dims[0];
    mNu_1 = (int)(dims[1] / nCols);
    if (dims[1] % nCols != 0 || mNu_1 == 0 || mTotalTimeStepsInFile < 2) {
        throw std::runtime_error("SourceTimeFunctionNetCDF::"
                                 "SourceTimeFunctionNetCDF || "
                                 "Invalid dimensions of STF data.");
    }
    
    // buffers
    if (!mAlignedToTimeStep && bufferSize < 2) {
        throw std::runtime_error("SourceTimeFunctionNetCDF::"
                                 "SourceTimeFunctionNetCDF || "
                                 "Buffer size must be greater than 2 "
                                 "for unaligned STF.");
    }
    // must fill time with min double
    mTimePoints.assign(bufferSize, std::numeric_limits<double>::lowest());
    mPatterns.resize(bufferSize);
    for (CTMatXND &pattern: mPatterns) {
        pattern.resize(mNu_1, nCols);
    }
    
    // static buffer for reading
    int totalSize = bufferSize * (int)dims[1];
    if (sReadRe.cols() < totalSize) {
        sReadRe.resize(totalSize);
        sReadIm.resize(totalSize);
    }
}

// load next buffer chunk
template <typename T, int ndim, int nCols> void
SourceTimeFunctionNetCDF<T, ndim, nCols>::loadNextBufferChunk() {
    // check remaining
    int haveBeenRead = mTimeStepOfPatternLast + 1;
    int remainInFile = mTotalTimeStepsInFile - haveBeenRead;
    if (remainInFile == 0) {
        throw std::runtime_error("SourceTimeFunctionNetCDF::"
                                 "loadNextBufferChunk || EOF reached.");
    }
    
    // update start
    if (mAlignedToTimeStep) {
        mTimeStepOfPattern0 = mTimeStepOfPatternLast + 1;
    } else {
        // special at 0
        if (mTimeStepOfPatternLast == -1) {
            mTimeStepOfPattern0 = 0;
        } else {
            mTimeStepOfPattern0 = mTimeStepOfPatternLast;
            remainInFile += 1;
        }
    }
    // count
    int readCount = std::min((int)mTimePoints.size(), remainInFile);
    mTimeStepOfPatternLast = mTimeStepOfPattern0 + readCount - 1;
    
    // read time
    mReader->readVariable(mTimeID, "time_points", mTimePoints,
                          {mTimeStepOfPattern0}, {readCount});
    
    // read data
    mReader->readVariable(mPatternReID, mVariableName + "_RE", sReadRe,
                          {mTimeStepOfPattern0, 0},
                          {readCount, mNu_1 * nCols});
    mReader->readVariable(mPatternImID, mVariableName + "_IM", sReadIm,
                          {mTimeStepOfPattern0, 0},
                          {readCount, mNu_1 * nCols});
    // copy data
    int pos = 0;
    for (int itime = 0; itime < readCount; itime++) {
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            mPatterns[itime].row(alpha).real() =
            sReadRe.block(0, pos, 1, nCols);
            mPatterns[itime].row(alpha).imag() =
            sReadIm.block(0, pos, 1, nCols);
            pos += nCols;
        }
    }
}
