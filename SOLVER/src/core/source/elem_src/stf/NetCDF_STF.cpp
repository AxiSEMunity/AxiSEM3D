//
//  NetCDF_STF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  netcdf source-time function

#include "NetCDF_STF.hpp"
#include "io.hpp"
#include "vector_tools.hpp"
#include "eigen_generic.hpp"
#include "bstring.hpp"

// constructor
NetCDF_STF::
NetCDF_STF(const std::string &fileName, const std::string &varTime,
           const std::string &varData, int chunkSize,
           PaddingMode padding, double left, double right):
mFileName(fileName), mVarTime(varTime), mVarData(varData),
mPadding(padding != PaddingMode::None),
mLeftPadding(left), mRightPadding(right) {
    // file
    if (sReaders.find(mFileName) == sReaders.end()) {
        sReaders.insert({mFileName, NetCDF_Reader()});
        sReaders.at(mFileName).open(io::popInputDir(mFileName));
    }
    
    // variable id
    mVarID_Time = sReaders.at(mFileName).getVariableID(mVarTime, mTimes);
    mVarID_Data = sReaders.at(mFileName).getVariableID(mVarData, mData);
    
    // read all data for verification
    std::vector<double> timesAll;
    std::vector<numerical::Real> dataAll;
    sReaders.at(mFileName).readVector(mVarTime, timesAll);
    sReaders.at(mFileName).readVector(mVarData, dataAll);
    mTotalTimeStepsInFile = (int)timesAll.size();
    chunkSize = std::min(chunkSize, mTotalTimeStepsInFile);
    
    // check time size
    if (chunkSize < 2) {
        throw std::
        runtime_error("NetCDF_STF::NetCDF_STF || "
                      "Insufficient time points, at least 2."
                      " || NetCDF file: " + io::popInputDir(mFileName) +
                      " || NetCDF variable for time: " + mVarTime +
                      " || NetCDF variable for data: " + mVarData);
    }
    // the time points must be strictly increasing
    if (!vector_tools::isSortedUnique(timesAll)) {
        throw std::
        runtime_error("NetCDF_STF::NetCDF_STF || "
                      "Time points are not ascendingly sorted."
                      " || NetCDF file: " + io::popInputDir(mFileName) +
                      " || NetCDF variable for time: " + mVarTime +
                      " || NetCDF variable for data: " + mVarData);
    }
    // size
    if (timesAll.size() != dataAll.size()) {
        throw std::
        runtime_error("NetCDF_STF::NetCDF_STF || "
                      "Times and values have different lengths."
                      " || NetCDF file: " + io::popInputDir(mFileName) +
                      " || NetCDF variable for time: " + mVarTime +
                      " || NetCDF variable for data: " + mVarData);
    }
    
    // range
    mTmin = timesAll.front();
    mTmax = timesAll.back();
    mVTmin = dataAll.front();
    mVTmax = dataAll.back();
    const eigen::RColX &vals =
    Eigen::Map<const eigen::RColX>(dataAll.data(), dataAll.size());
    mVmin = vals.minCoeff();
    mVmax = vals.maxCoeff();
    
    // buffer (time must be filled with lowest)
    mTimes.assign(chunkSize, std::numeric_limits<double>::lowest());
    mData.assign(chunkSize, 0.);
    
    // padding
    if (padding == PaddingMode::FirstLast) {
        mLeftPadding = dataAll.front();
        mRightPadding = dataAll.back();
    }
}

// get value
numerical::Real NetCDF_STF::getValue(double time) {
    // check padding first
    if (mPadding) {
        if (time < mTmin) {
            return mLeftPadding;
        }
        if (time > mTmax) {
            return mRightPadding;
        }
    }
    
    // check buffer
    while (time > mTimes[mTimeStepOfTimesLast -
                         mTimeStepOfTimes0]) {
        if (!loadNextBufferChunk()) {
            throw std::runtime_error
            ("NetCDF_STF::getValue || "
             "Time is out of range with padding disabled."
             " || Requested time = " + bstring::toString(time) +
             " || STF time range = " + bstring::range(mTmin, mTmax) +
             " || NetCDF file: " + io::popInputDir(mFileName) +
             " || NetCDF variable for time: " + mVarTime +
             " || NetCDF variable for data: " + mVarData);
        }
    }
    try {
        // interpolate in time
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        vector_tools::linearInterpSorted(mTimes, time, index0, index1,
                                         factor0, factor1, 0,
                                         mTimeStepOfTimesLast -
                                         mTimeStepOfTimes0 + 1);
        return (numerical::Real)(factor0 * mData[index0] +
                                 factor1 * mData[index1]);
    } catch (...) {
        throw std::runtime_error
        ("NetCDF_STF::getValue || "
         "Time is out of range with padding disabled."
         " || Requested time = " + bstring::toString(time) +
         " || STF time range = " + bstring::range(mTmin, mTmax) +
         " || NetCDF file: " + io::popInputDir(mFileName) +
         " || NetCDF variable for time: " + mVarTime +
         " || NetCDF variable for data: " + mVarData);
    }
}

// verbose
std::string NetCDF_STF::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxSubTitle(2, "Source-time function");
    ss << boxEquals(4, 18, "class name", "NetCDF_STF");
    ss << boxEquals(4, 18, "NetCDF file", mFileName);
    ss << boxEquals(4, 18, "NC-var for time", mVarTime);
    ss << boxEquals(4, 18, "NC-var for data", mVarData);
    ss << boxEquals(4, 18, "# time points", mTotalTimeStepsInFile);
    ss << boxEquals(4, 18, "chunk size", mTimes.size());
    // range
    ss << boxEquals(4, 18, "1st (time, value)",
                    range(mTmin, mVTmin, '(', ')'));
    ss << boxEquals(4, 18, "last (time, value)",
                    range(mTmax, mVTmax, '(', ')'));
    ss << boxEquals(4, 18, "[min, max] values", range(mVmin, mVmax));
    // padding
    if (mPadding) {
        ss << boxEquals(4, 18, "padding", range(mLeftPadding, mRightPadding));
    } else {
        ss << boxEquals(4, 18, "padding", "disabled");
    }
    return ss.str();
}

// load next buffer chunk
bool NetCDF_STF::loadNextBufferChunk() {
    // check remaining
    int haveBeenRead = mTimeStepOfTimesLast + 1;
    int remainInFile = mTotalTimeStepsInFile - haveBeenRead;
    if (remainInFile == 0) {
        // EOF
        return false;
    }
    
    // update start
    if (mTimeStepOfTimesLast == -1) {
        mTimeStepOfTimes0 = 0;
    } else {
        mTimeStepOfTimes0 = mTimeStepOfTimesLast;
        remainInFile += 1;
    }
    // count
    int readCount = std::min((int)mTimes.size(), remainInFile);
    mTimeStepOfTimesLast = mTimeStepOfTimes0 + readCount - 1;
    
    // read time
    sReaders.at(mFileName).readVariable(mVarID_Time, mVarTime, mTimes,
                                        {mTimeStepOfTimes0}, {readCount});
    sReaders.at(mFileName).readVariable(mVarID_Data, mVarData, mData,
                                        {mTimeStepOfTimes0}, {readCount});
    return true;
}
