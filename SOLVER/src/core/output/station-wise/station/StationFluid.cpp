//
//  StationFluid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station in fluid

#include "StationFluid.hpp"
#include "FluidElement.hpp"
#include "geodesy.hpp"

/////////////////////////// setup ///////////////////////////
// set element
void StationFluid::
setElement(const std::shared_ptr<FluidElement> &element,
           const eigen::DRowN &weights) {
    // element
    mElement = element;
    
    // base
    Station::setElement(weights, mElement->getNu_1());
}

// set in group
void StationFluid::
setInGroup(int dumpIntv, const channel::fluid::ChannelOptions &chops) {
    // member buffers
    if (chops.mNeedBufferX) {
        mBufferX.resize(dumpIntv, 1);
    }
    if (chops.mNeedBufferU) {
        mBufferU.resize(dumpIntv, 3);
    }
    if (chops.mNeedBufferP) {
        mBufferP.resize(dumpIntv, 1);
    }
    if (chops.mNeedBufferD) {
        mBufferD.resize(dumpIntv, 1);
    }
    
    // element
    mElement->prepareWavefieldOutput(chops, false);
    
    // workspace
    expandWorkspaceRecord(mElement->getNu_1(), chops);
    expandWorkspaceProcess(dumpIntv, false);
}


/////////////////////////// record ///////////////////////////
// record
void StationFluid::
record(int bufferLine, const channel::fluid::ChannelOptions &chops) {
    int nu_1 = mElement->getNu_1();
    // chi
    if (chops.mNeedBufferX) {
        mElement->getChiField(sXXN1);
        interpolate<1>(sXXN1, sXX1, sX1, nu_1);
        mBufferX.row(bufferLine) = sX1;
    }
    // displacement
    if (chops.mNeedBufferU) {
        mElement->getDisplField(sUXN3);
        interpolate<3>(sUXN3, sUX3, sU3, nu_1);
        mBufferU.row(bufferLine) = sU3;
    }
    // pressure
    if (chops.mNeedBufferP) {
        mElement->getPressureField(sPXN1);
        interpolate<1>(sPXN1, sPX1, sP1, nu_1);
        mBufferP.row(bufferLine) = sP1;
    }
    // delta
    if (chops.mNeedBufferD) {
        mElement->getDeltaField(sDXN1);
        interpolate<1>(sDXN1, sDX1, sD1, nu_1);
        mBufferD.row(bufferLine) = sD1;
    }
}


/////////////////////////// process ///////////////////////////
// process and report to group
void StationFluid::
processReport(int bufferLine,
              const channel::fluid::ChannelOptions &chops,
              int stationIndex, eigen::RTensor3 &bufferFields) {
    // rotate
    rotate(bufferLine, chops);
    
    // channels
    for (int ich = 0; ich < chops.mStdChannels.size(); ich++) {
        // find field and index of channel
        const int &cha = chops.mStdChannels[ich];
        const auto &tup = channel::fluid::gChannelMap.at(cha);
        channel::fluid::FieldType ftype = std::get<2>(tup);
        int fieldIndex = std::get<3>(tup);
        // compute and feed
        if (ftype == channel::fluid::FieldType::Chi) {
            computeFeedChannel<1>(mBufferX, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::fluid::FieldType::Displ) {
            computeFeedChannel<3>(mBufferU, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::fluid::FieldType::Pressure) {
            computeFeedChannel<1>(mBufferP, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::fluid::FieldType::Delta) {
            computeFeedChannel<1>(mBufferD, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else {
            throw std::runtime_error("StationFluid::processReport || "
                                     "Unknown field type.");
        }
    }
}

// process 1: rotate
void StationFluid::
rotate(int bufferLine, const channel::fluid::ChannelOptions &chops) {
    bool cartesian = geodesy::isCartesian();
    if (chops.mNeedBufferU) {
        rotateField<3>(mBufferU, bufferLine, mElement->displInRTZ(),
                       chops.mWCS, cartesian);
    }
}
