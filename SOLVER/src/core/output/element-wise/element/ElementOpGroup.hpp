//
//  ElementOpGroup.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output group

#ifndef ElementOpGroup_hpp
#define ElementOpGroup_hpp

#include "ElementIO.hpp"
#include "ElementOpSolid.hpp"
#include "ElementOpFluid.hpp"
#include "vector_tools.hpp"
#include "mpi.hpp"

///////////////////// type inference /////////////////////
template <class ElementOpT>
struct ElementOpInference {
};

template <>
struct ElementOpInference<ElementOpSolid> {
    typedef channel::solid::ChannelOptions ChannelOptions;
    inline static const auto &gChannelMap = channel::solid::gChannelMap;
};

template <>
struct ElementOpInference<ElementOpFluid> {
    typedef channel::fluid::ChannelOptions ChannelOptions;
    inline static const auto &gChannelMap = channel::fluid::gChannelMap;
};


///////////////////// element group /////////////////////
template <class ElementOpT>
class ElementOpGroup {
public:
    // constructor
    ElementOpGroup(const std::string &groupName,
                   int numRecordSteps,
                   int sampleIntv, double tmin, double tmax,
                   int dumpIntv, channel::WavefieldCS wcs,
                   const std::vector<std::string> &userChannels,
                   int npnts, const std::vector<double> &phis, int naSpace,
                   std::unique_ptr<ElementIO> &io):
    mGroupName(groupName), mNumRecordSteps(numRecordSteps),
    mSampleIntv(sampleIntv), mTmin(tmin), mTmax(tmax),
    mDumpIntv(std::min(dumpIntv, numRecordSteps)),
    mChannelOptions(wcs, userChannels), mNPnts(npnts), mPhis(phis),
    mNaSpace(naSpace), mIO(io.release()) {
        // nothing
    }
    
    // add an element
    void addElementOp(std::unique_ptr<ElementOpT> &elementOp) {
        mElementOps.push_back(std::move(elementOp));
    }
    
    // initialize after adding all elements
    void initialize() {
        // channels
        int nch = (int)mChannelOptions.mStdChannels.size();
        std::vector<std::string> channels;
        channels.reserve(nch);
        for (int ich = 0; ich < nch; ich++) {
            const std::string &chstr =
            std::get<0>(ElementOpInference<ElementOpT>::gChannelMap.
                        at(mChannelOptions.mStdChannels[ich]));
            channels.push_back(chstr);
        }
        
        // Fourier
        int nphis = (int)mPhis.size();
        if (nphis > 0) {
            // find max nu_1
            int maxNu_1 = 0;
            for (const std::unique_ptr<ElementOpT> &eop: mElementOps) {
                maxNu_1 = std::max(maxNu_1, eop->getNu_1());
            }
            // allocate and compute
            mExpIAlphaPhi = eigen::CMatXX::Zero(maxNu_1, nphis);
            for (int iphi = 0; iphi < nphis; iphi++) {
                eigen::CColX temp(maxNu_1);
                eigen_tools::computeTwoExpIAlphaPhi(maxNu_1, mPhis[iphi], temp);
                mExpIAlphaPhi.col(iphi) = temp;
            }
        }
        
        //////////// na grid ////////////
        // all na on elements
        std::vector<int> allNa, sortedNa, allTag;
        int nelem = (int)mElementOps.size();
        allNa.reserve(nelem);
        allTag.reserve(nelem);
        for (const std::unique_ptr<ElementOpT> &eop: mElementOps) {
            allNa.push_back(eop->getNa(nphis));
            allTag.push_back(eop->getElementTag());
        }
        sortedNa = allNa;
        vector_tools::sortUnique(sortedNa);
        
        // merge the sorted na
        mNaGrid.clear();
        if (sortedNa.size() > 0) {
            mNaGrid.push_back(sortedNa.back());
            for (int isort = (int)sortedNa.size() - 2; isort >= 0; isort--) {
                if (sortedNa[isort] <= mNaGrid.back() - mNaSpace) {
                    mNaGrid.push_back(sortedNa[isort]);
                }
            }
        }
        // reverse order (descending to ascending)
        vector_tools::sortUnique(mNaGrid);
        
        // na-grid must be the same on ranks
        std::vector<std::vector<int>> naGridGlobal;
        mpi::gather(mNaGrid, naGridGlobal, MPI_ALL);
        mNaGrid.clear();
        for (int irank = 0; irank < mpi::nproc(); irank++) {
            mNaGrid.insert(mNaGrid.end(), naGridGlobal[irank].begin(),
                           naGridGlobal[irank].end());
        }
        vector_tools::sortUnique(mNaGrid);
        mNaGrid.shrink_to_fit();
        
        // form dict of na-grid for fast search
        for (int inag = 0; inag < mNaGrid.size(); inag++) {
            mNaGridIndexDict[mNaGrid[inag]] = inag;
        }
        
        //////////// element-na info and coords ////////////
        std::vector<int> numElemNaGrid = std::vector<int>(mNaGrid.size(), 0);
        mElemNaInfo = eigen::IMatX4_RM(nelem, 4);
        eigen::DMatXX_RM elemCoords(nelem, mNPnts * 2);
        for (int ielem = 0; ielem < nelem; ielem++) {
            mElemNaInfo(ielem, 0) = allTag[ielem];
            mElemNaInfo(ielem, 1) = allNa[ielem];
            // find na gird for this na
            auto it = std::lower_bound(mNaGrid.begin(), mNaGrid.end(),
                                       allNa[ielem]);
            mElemNaInfo(ielem, 2) = *it;
            int naGridIndex = (int)(it - mNaGrid.begin());
            mElemNaInfo(ielem, 3) = numElemNaGrid[naGridIndex];
            // increment element number of this grid-na
            numElemNaGrid[naGridIndex]++;
            // coords
            elemCoords.row(ielem) = mElementOps[ielem]->getCoords();
        }
        
        // buffers
        if (nelem > 0) {
            // time buffer is not needed without stations
            mBufferTime.resize(mDumpIntv);
        }
        for (int inag = 0; inag < mNaGrid.size(); inag++) {
            eigen::RTensor5 bufferField(numElemNaGrid[inag],
                                        mNaGrid[inag], mNPnts,
                                        nch, mDumpIntv);
            mBufferFields.push_back(bufferField);
        }
        mBufferLine = 0;
        
        // IO
        mIO->initialize(mGroupName, mNumRecordSteps, channels,
                        mNPnts, mNaGrid, mElemNaInfo, elemCoords);
        
        // elementOps
        for (const std::unique_ptr<ElementOpT> &eop: mElementOps) {
            eop->setInGroup(mDumpIntv, mChannelOptions, (int)mPhis.size());
        }
    }
    
    // finalize after time loop
    void finalize() const {
        mIO->finalize();
    }
    
    // record at a time step
    void record(int timeStep, double time) {
        // interval and time window
        if (timeStep % mSampleIntv != 0 || time < mTmin || time > mTmax) {
            return;
        }
        
        // time
        if (mElementOps.size() > 0) {
            mBufferTime(mBufferLine) = time;
        }
        
        // stations
        for (const std::unique_ptr<ElementOpT> &eop: mElementOps) {
            eop->record(mBufferLine, mChannelOptions, mExpIAlphaPhi);
        }
        
        // increment buffer line
        mBufferLine++;
        
        // dump to file
        if (mBufferLine == mDumpIntv) {
            dumpToFile();
        }
    }
    
    // dump to file
    void dumpToFile() {
        // redundant call at the end of the time loop
        if (mBufferLine == 0) {
            return;
        }
        
        // stations
        for (int ie = 0; ie < mElementOps.size(); ie++) {
            int naGridIndex = mNaGridIndexDict.at(mElemNaInfo(ie, 2));
            int elemIndexNaGrid = mElemNaInfo(ie, 3);
            mElementOps[ie]->processReport(mBufferLine, mChannelOptions,
                                           elemIndexNaGrid, naGridIndex,
                                           mBufferFields);
        }
        
        // IO
        // check zero station inside
        mIO->dumpToFile(mBufferTime, mBufferFields, mBufferLine);
        
        // rewind buffer line
        mBufferLine = 0;
    }
    
private:
    //////////////// const ////////////////
    // name
    const std::string mGroupName;
    
    // steps
    const int mNumRecordSteps;
    const int mSampleIntv;
    const double mTmin, mTmax;
    const int mDumpIntv;
    
    // channels
    const typename
    ElementOpInference<ElementOpT>::ChannelOptions mChannelOptions;
    
    // inplane
    const int mNPnts;
    
    // azimuth
    const std::vector<double> mPhis;
    
    // na space
    const int mNaSpace;
    
    // IO
    const std::unique_ptr<ElementIO> mIO;
    
    //////////////// non-const ////////////////
    // elements
    std::vector<std::unique_ptr<ElementOpT>> mElementOps;
    
    // Fourier
    eigen::CMatXX mExpIAlphaPhi = eigen::CMatXX(0, 0);
    
    // na
    std::vector<int> mNaGrid;
    std::map<int, int> mNaGridIndexDict;
    eigen::IMatX4_RM mElemNaInfo;
    
    // buffer
    eigen::DColX mBufferTime = eigen::DColX(0);
    std::vector<eigen::RTensor5> mBufferFields;
    int mBufferLine = 0;
};

#endif /* ElementOpGroup_hpp */
