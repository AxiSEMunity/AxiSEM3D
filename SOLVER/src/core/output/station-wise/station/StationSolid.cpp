//
//  StationSolid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station in solid

#include "StationSolid.hpp"
#include "SolidElement.hpp"
#include "geodesy.hpp"

/////////////////////////// setup ///////////////////////////
// set element
void StationSolid::
setElement(const std::shared_ptr<SolidElement> &element,
           const eigen::DRowN &weights) {
    // element
    mElement = element;
    
    // base
    Station::setElement(weights, mElement->getNu_1());
}

// set in group
void StationSolid::
setInGroup(int dumpIntv, const channel::solid::ChannelOptions &chops) {
    // member buffers
    if (chops.mNeedBufferU) {
        mBufferU.resize(dumpIntv, 3);
    }
    if (chops.mNeedBufferG) {
        mBufferG.resize(dumpIntv, 9);
    }
    if (chops.mNeedBufferE) {
        mBufferE.resize(dumpIntv, 6);
    }
    if (chops.mNeedBufferR) {
        mBufferR.resize(dumpIntv, 3);
    }
    if (chops.mNeedBufferS) {
        mBufferS.resize(dumpIntv, 6);
    }
    
    // element
    mElement->prepareWavefieldOutput(chops, false);
    
    // workspace
    expandWorkspaceRecord(mElement->getNu_1(), chops);
    expandWorkspaceProcess(dumpIntv, chops.mNeedBufferE || chops.mNeedBufferS);
}


/////////////////////////// record ///////////////////////////
// record
void StationSolid::
record(int bufferLine, const channel::solid::ChannelOptions &chops) {
    int nu_1 = mElement->getNu_1();
    // displ
    if (chops.mNeedBufferU) {
        mElement->getDisplField(sUXN3);
        interpolate<3>(sUXN3, sUX3, sU3, nu_1);
        mBufferU.row(bufferLine) = sU3;
    }
    // nabla
    if (chops.mNeedBufferG) {
        mElement->getNablaField(sGXN9);
        interpolate<9>(sGXN9, sGX9, sG9, nu_1);
        mBufferG.row(bufferLine) = sG9;
    }
    // strain
    if (chops.mNeedBufferE) {
        mElement->getStrainField(sEXN6);
        interpolate<6>(sEXN6, sEX6, sE6, nu_1);
        mBufferE.row(bufferLine) = sE6;
    }
    // curl
    if (chops.mNeedBufferR) {
        mElement->getCurlField(sRXN3);
        interpolate<3>(sRXN3, sRX3, sR3, nu_1);
        mBufferR.row(bufferLine) = sR3;
    }
    // stress
    if (chops.mNeedBufferS) {
        mElement->getStressField(sSXN6);
        interpolate<6>(sSXN6, sSX6, sS6, nu_1);
        mBufferS.row(bufferLine) = sS6;
    }
}


/////////////////////////// process ///////////////////////////
// process and report to group
void StationSolid::
processReport(int bufferLine,
              const channel::solid::ChannelOptions &chops,
              int stationIndex, eigen::RTensor3 &bufferFields) {
    // rotate
    rotate(bufferLine, chops);
    
    // loop over channels
    for (int ich = 0; ich < chops.mStdChannels.size(); ich++) {
        // find field and index of channel
        const int &cha = chops.mStdChannels[ich];
        const auto &tup = channel::solid::gChannelMap.at(cha);
        channel::solid::FieldType ftype = std::get<2>(tup);
        int fieldIndex = std::get<3>(tup);
        // compute and feed
        if (ftype == channel::solid::FieldType::Displ) {
            computeFeedChannel<3>(mBufferU, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Nabla) {
            computeFeedChannel<9>(mBufferG, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Strain) {
            computeFeedChannel<6>(mBufferE, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Curl) {
            computeFeedChannel<3>(mBufferR, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else if (ftype == channel::solid::FieldType::Stress) {
            computeFeedChannel<6>(mBufferS, fieldIndex,
                                  bufferLine, ich, stationIndex,
                                  bufferFields);
        } else {
            throw std::runtime_error("StationSolid::processReport || "
                                     "Unknown field type.");
        }
    }
}

// process 1: rotate
void StationSolid::
rotate(int bufferLine, const channel::solid::ChannelOptions &chops) {
    bool cartesian = geodesy::isCartesian();
    if (chops.mNeedBufferU) {
        rotateField<3>(mBufferU, bufferLine, mElement->displInRTZ(),
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferG) {
        rotateField<9>(mBufferG, bufferLine, mElement->nablaInRTZ(),
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferE) {
        // halve off-diagonal components before rotation
        mBufferE.rightCols(3) *= (numerical::Real).5;
        rotateField<6>(mBufferE, bufferLine, mElement->strainInRTZ(),
                       chops.mWCS, cartesian);
        mBufferE.rightCols(3) *= (numerical::Real)2.;
    }
    if (chops.mNeedBufferR) {
        rotateField<3>(mBufferR, bufferLine, mElement->curlInRTZ(),
                       chops.mWCS, cartesian);
    }
    if (chops.mNeedBufferS) {
        rotateField<6>(mBufferS, bufferLine, mElement->stressInRTZ(),
                       chops.mWCS, cartesian);
    }
}
