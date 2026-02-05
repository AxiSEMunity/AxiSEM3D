//
//  ElementOpSolid.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output in solid

#include "ElementOpSolid.hpp"

/////////////////////////// setup ///////////////////////////
// set in group
void ElementOpSolid::
setInGroup(int dumpIntv, const channel::solid::ChannelOptions &chops,
           int nphis) {
    // sizes
    int na = getNa(nphis);
    int npnts = (int)mIPnts.size();
    
    // member buffers
    if (chops.mNeedBufferU) {
        mBufferU.resize(na, npnts, 3, dumpIntv);
    }
    if (chops.mNeedBufferG) {
        mBufferG.resize(na, npnts, 9, dumpIntv);
    }
    if (chops.mNeedBufferE) {
        mBufferE.resize(na, npnts, 6, dumpIntv);
    }
    if (chops.mNeedBufferR) {
        mBufferR.resize(na, npnts, 3, dumpIntv);
    }
    if (chops.mNeedBufferS) {
        mBufferS.resize(na, npnts, 6, dumpIntv);
    }
    
    // element
    mElement->prepareWavefieldOutput(chops, true);
    
    // workspace
    expandWorkspaceRecord(mElement->getNu_1(), na, chops);
}


/////////////////////////// record ///////////////////////////
// record
void ElementOpSolid::
record(int bufferLine, const channel::solid::ChannelOptions &chops,
       const eigen::CMatXX &expIAlphaPhi) {
    int nu_1 = mElement->getNu_1();
    bool needRTZ = (chops.mWCS == channel::WavefieldCS::RTZ);
    
    if (chops.mNeedBufferU) {
        mElement->getDisplField(sCUXN3, needRTZ);
        recordToElem<3>(sCUXN3, nu_1, sRUXN3, expIAlphaPhi,
                        mBufferU, bufferLine);
    }
    if (chops.mNeedBufferG) {
        mElement->getNablaField(sCGXN9, needRTZ);
        recordToElem<9>(sCGXN9, nu_1, sRGXN9, expIAlphaPhi,
                        mBufferG, bufferLine);
    }
    if (chops.mNeedBufferE) {
        mElement->getStrainField(sCEXN6, needRTZ);
        recordToElem<6>(sCEXN6, nu_1, sREXN6, expIAlphaPhi,
                        mBufferE, bufferLine);
        
    }
    if (chops.mNeedBufferR) {
        mElement->getCurlField(sCRXN3, needRTZ);
        recordToElem<3>(sCRXN3, nu_1, sRRXN3, expIAlphaPhi,
                        mBufferR, bufferLine);
    }
    if (chops.mNeedBufferS) {
        mElement->getStressField(sCSXN6, needRTZ);
        recordToElem<6>(sCSXN6, nu_1, sRSXN6, expIAlphaPhi,
                        mBufferS, bufferLine);
    }
}


/////////////////////////// process ///////////////////////////
// process and report to group
void ElementOpSolid::
processReport(int bufferLine, const channel::solid::ChannelOptions &chops,
              int elemIndexNaGrid, int naGridIndex,
              std::vector<eigen::RTensor5> &ioBuffers) {
    // loop over channels
    for (int ich = 0; ich < chops.mStdChannels.size(); ich++) {
        // find field and index of channel
        const int &cha = chops.mStdChannels[ich];
        const auto &tup = channel::solid::gChannelMap.at(cha);
        channel::solid::FieldType ftype = std::get<2>(tup);
        int fieldIndex = std::get<3>(tup);
        
        // compute and feed
        if (ftype == channel::solid::FieldType::Displ) {
            dumpToIO<3>(mBufferU, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::solid::FieldType::Nabla) {
            dumpToIO<9>(mBufferG, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::solid::FieldType::Strain) {
            dumpToIO<6>(mBufferE, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::solid::FieldType::Curl) {
            dumpToIO<3>(mBufferR, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else if (ftype == channel::solid::FieldType::Stress) {
            dumpToIO<6>(mBufferS, fieldIndex,
                        bufferLine, ich, elemIndexNaGrid,
                        naGridIndex, ioBuffers);
        } else {
            throw std::runtime_error("ElementOpSolid::processReport || "
                                     "Unknown field type.");
        }
    }
}
