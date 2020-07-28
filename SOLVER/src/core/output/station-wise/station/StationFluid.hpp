//
//  StationFluid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station in fluid

#ifndef StationFluid_hpp
#define StationFluid_hpp

#include "Station.hpp"
class FluidElement;

class StationFluid: public Station {
public:
    // constructor
    StationFluid(const std::string &key, double phi,
                 double theta, double backAzimuth):
    Station(key, phi, theta, backAzimuth) {
        // nothing
    }
    
    
    /////////////////////////// setup ///////////////////////////
    // set element
    void setElement(const std::shared_ptr<FluidElement> &element,
                    const eigen::DRowN &weights);
    
    // set in group
    void setInGroup(int dumpIntv, const channel::fluid::ChannelOptions &chops);
    
    
    /////////////////////////// record ///////////////////////////
public:
    // record
    void record(int bufferLine, const channel::fluid::ChannelOptions &chops);
    
private:
    // element
    std::shared_ptr<FluidElement> mElement = nullptr;
    
    // buffer
    eigen::RMatX1_RM mBufferX = eigen::RMatX1_RM(0, 1);
    eigen::RMatX3_RM mBufferU = eigen::RMatX3_RM(0, 3);
    eigen::RMatX1_RM mBufferP = eigen::RMatX1_RM(0, 1);
    eigen::RMatX1_RM mBufferD = eigen::RMatX1_RM(0, 1);
    
    
    /////////////////////////// process ///////////////////////////
public:
    // process and report to group
    void processReport(int bufferLine,
                       const channel::fluid::ChannelOptions &chops,
                       int stationIndex, eigen::RTensor3 &bufferFields);
    
private:
    // process 1: rotate
    void rotate(int bufferLine, const channel::fluid::ChannelOptions &chops);
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // expand workspace for record
    static void
    expandWorkspaceRecord(int nu_1,
                          const channel::fluid::ChannelOptions &chops) {
        // for element
        if (chops.mNeedBufferX && sXXN1.rows() < nu_1) {
            sXXN1.resize(nu_1, spectral::nPEM * 1);
            sXX1.resize(nu_1, 1);
        }
        
        if (chops.mNeedBufferU && sUXN3.rows() < nu_1) {
            sUXN3.resize(nu_1, spectral::nPEM * 3);
            sUX3.resize(nu_1, 3);
        }
        
        if (chops.mNeedBufferP && sPXN1.rows() < nu_1) {
            sPXN1.resize(nu_1, spectral::nPEM * 1);
            sPX1.resize(nu_1, 1);
        }
        
        if (chops.mNeedBufferD && sDXN1.rows() < nu_1) {
            sDXN1.resize(nu_1, spectral::nPEM * 1);
            sDX1.resize(nu_1, 1);
        }
    }
    
    // workspace for record
    inline static eigen::RRow1 sX1;
    inline static eigen::RRow3 sU3;
    inline static eigen::RRow1 sP1;
    inline static eigen::RRow1 sD1;
    inline static eigen::CMatX1 sXX1 = eigen::CMatX1(0, 1);
    inline static eigen::CMatX3 sUX3 = eigen::CMatX3(0, 3);
    inline static eigen::CMatX1 sPX1 = eigen::CMatX1(0, 1);
    inline static eigen::CMatX1 sDX1 = eigen::CMatX1(0, 1);
    inline static
    eigen::CMatXN sXXN1 = eigen::CMatXN(0, spectral::nPEM * 1);
    inline static
    eigen::CMatXN3 sUXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static
    eigen::CMatXN sPXN1 = eigen::CMatXN(0, spectral::nPEM * 1);
    inline static
    eigen::CMatXN sDXN1 = eigen::CMatXN(0, spectral::nPEM * 1);
};

#endif /* StationFluid_hpp */
