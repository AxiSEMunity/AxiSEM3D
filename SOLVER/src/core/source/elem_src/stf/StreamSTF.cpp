//
//  StreamSTF.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  stream source-time function

#include "StreamSTF.hpp"
#include "io.hpp"
#include "vector_tools.hpp"
#include "eigen_generic.hpp"
#include "bstring.hpp"

// constructor
StreamSTF::StreamSTF(const std::string &fileName, PaddingMode padding,
                     numerical::Real left, numerical::Real right):
mFileName(fileName), mPadding(padding != PaddingMode::None),
mLeftPadding(left), mRightPadding(right) {
    // open and read
    const std::vector<std::string> &lines =
    bstring::readLines(io::popInputDir(mFileName), "StreamSTF::StreamSTF");
    // cast to vector
    for (const std::string &line: lines) {
        const std::vector<std::string> &words = bstring::split(line, " \t");
        // empty or comment line
        if (words[0].front() == '#' || words[0] == "") {
            continue;
        }
        // size check
        if (words.size() != 2) {
            throw std::
            runtime_error("StreamSTF::StreamSTF || "
                          "Input file for STF must contain two columns."
                          " || Ascii file: " + io::popInputDir(mFileName));
        }
        // valid line
        mTimes.push_back(bstring::cast<double>
                         (words[0], "StreamSTF::StreamSTF"));
        mData.push_back(bstring::cast<numerical::Real>
                        (words[1], "StreamSTF::StreamSTF"));
    }
    
    // check time size
    if (mTimes.size() < 2) {
        throw std::
        runtime_error("StreamSTF::StreamSTF || "
                      "Insufficient time points, at least 2."
                      " || Ascii file: " + io::popInputDir(mFileName));
    }
    // the time points must be strictly increasing
    if (!vector_tools::isSortedUnique(mTimes)) {
        throw std::
        runtime_error("StreamSTF::StreamSTF || "
                      "Time points are not ascendingly sorted."
                      " || Ascii file: " + io::popInputDir(mFileName));
    }
    
    // padding
    if (padding == PaddingMode::FirstLast) {
        mLeftPadding = mData.front();
        mRightPadding = mData.back();
    }
}

// get value
numerical::Real StreamSTF::getValue(double time) {
    try {
        // interpolate in time
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        vector_tools::linearInterpSorted(mTimes, time, index0, index1,
                                         factor0, factor1);
        return (numerical::Real)(factor0 * mData[index0] +
                                 factor1 * mData[index1]);
    } catch (...) {
        if (!mPadding) {
            throw std::
            runtime_error("StreamSTF::getValue || "
                          "Time is out of range with padding disabled."
                          " || Requested time = " + bstring::toString(time) +
                          " || STF time range = " +
                          bstring::range(mTimes.front(), mTimes.back()) +
                          " || Ascii file: " + io::popInputDir(mFileName));
        }
        // use padding if time is out of range
        if (time < mTimes.front()) {
            return mLeftPadding;
        } else if (time > mTimes.back()) {
            return mRightPadding;
        } else {
            throw std::
            runtime_error("StreamSTF::getValue || "
                          "Program should not reach here."
                          " || Ascii file: " + io::popInputDir(mFileName));
        }
    }
}

// verbose
std::string StreamSTF::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxSubTitle(2, "Source-time function");
    ss << boxEquals(4, 18, "class name", "StreamSTF");
    ss << boxEquals(4, 18, "Ascii file", mFileName);
    // range
    ss << boxEquals(4, 18, "1st (time, value)",
                    range(mTimes.front(), (double)mData.front(), '(', ')'));
    ss << boxEquals(4, 18, "last (time, value)",
                    range(mTimes.back(), (double)mData.back(), '(', ')'));
    const eigen::RColX &d =
    Eigen::Map<const eigen::RColX>(mData.data(), mData.size());
    ss << boxEquals(4, 18, "[min, max] values",
                    range(d.minCoeff(), d.maxCoeff()));
    // padding
    if (mPadding) {
        ss << boxEquals(4, 18, "padding", range(mLeftPadding, mRightPadding));
    } else {
        ss << boxEquals(4, 18, "padding", "disabled");
    }
    return ss.str();
}
