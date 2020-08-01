//
//  ElementOpFluid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output in fluid

#include "ElementOpFluid.hpp"

/////////////////////////// setup ///////////////////////////
// set in group
void ElementOpFluid::
setInGroup(int dumpIntv, const channel::fluid::ChannelOptions &chops,
           int nphis) {
    // sizes
    int na = getNa(nphis);
    int npnts = (int)mIPnts.size();
    
    // member buffers
    if (chops.mNeedBufferX) {
        mBufferX.resize(na, npnts, 1, dumpIntv);
    }
    if (chops.mNeedBufferU) {
        mBufferU.resize(na, npnts, 3, dumpIntv);
    }
    if (chops.mNeedBufferP) {
        mBufferP.resize(na, npnts, 1, dumpIntv);
    }
    if (chops.mNeedBufferD) {
        mBufferD.resize(na, npnts, 1, dumpIntv);
    }
    
    // element
    mElement->prepareWavefieldOutput(chops, true);
    
    // workspace
    expandWorkspaceRecord(mElement->getNu_1(), na, chops);
}


/////////////////////////// record ///////////////////////////
// record
void ElementOpFluid::
record(int bufferLine, const channel::fluid::ChannelOptions &chops,
       const eigen::CMatXX &expIAlphaPhi) {
    int nu_1 = mElement->getNu_1();
    bool needRTZ = (chops.mWCS == channel::WavefieldCS::RTZ);
    
    if (chops.mNeedBufferX) {
        mElement->getChiField(sCXXN1);
        recordToElem<1>(sCXXN1, nu_1, sRXXN1, expIAlphaPhi,
                        mBufferX, bufferLine);
    }
    if (chops.mNeedBufferU) {
        mElement->getDisplField(sCUXN3, needRTZ);
        recordToElem<3>(sCUXN3, nu_1, sRUXN3, expIAlphaPhi,
                        mBufferU, bufferLine);
    }
    if (chops.mNeedBufferP) {
        mElement->getPressureField(sCPXN1);
        recordToElem<1>(sCPXN1, nu_1, sRPXN1, expIAlphaPhi,
                        mBufferP, bufferLine);
    }
    if (chops.mNeedBufferD) {
        mElement->getDeltaField(sCDXN1);
        recordToElem<1>(sCDXN1, nu_1, sRDXN1, expIAlphaPhi,
                        mBufferD, bufferLine);
    }
}


/////////////////////////// process ///////////////////////////
// process and report to group
void ElementOpFluid::
processReport(int bufferLine, const channel::fluid::ChannelOptions &chops,
              int elemIndexNaGrid, int naGridIndex,
              std::vector<eigen::RTensor5> &ioBuffers) {
    // loop over channels
    for (int ich = 0; ich < chops.mStdChannels.size(); ich++) {
        // find field and index of channel
        const int &cha = chops.mStdChannels[ich];
        const auto &tup = channel::fluid::gChannelMap.at(cha);
        channel::fluid::FieldType ftype = std::get<2>(tup);
        int fieldIndex = std::get<3>(tup);
        
        // compute and feed
        if (ftype == channel::fluid::FieldType::Chi) {
            dumpToIO<1>(mBufferX, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::fluid::FieldType::Displ) {
            dumpToIO<3>(mBufferU, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::fluid::FieldType::Pressure) {
            dumpToIO<1>(mBufferP, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::fluid::FieldType::Delta) {
            dumpToIO<1>(mBufferD, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else {
            throw std::runtime_error("ElementOpFluid::processReport || "
                                     "Unknown field type.");
        }
    }
}
