//
//  StationSolid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/3/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station in solid

#ifndef StationSolid_hpp
#define StationSolid_hpp

#include "Station.hpp"
class SolidElement;

class StationSolid: public Station {
public:
    // constructor
    StationSolid(const std::string &key, double phi,
                 double theta, double backAzimuth):
    Station(key, phi, theta, backAzimuth) {
        // nothing
    }
    
    
    /////////////////////////// setup ///////////////////////////
    // set element
    void setElement(const std::shared_ptr<SolidElement> &element,
                    const eigen::DRowN &weights);
    
    // set in group
    void setInGroup(int dumpIntv, const channel::solid::ChannelOptions &chops);
    
    
    /////////////////////////// record ///////////////////////////
public:
    // record
    void record(int bufferLine, const channel::solid::ChannelOptions &chops);
    
private:
    // element
    std::shared_ptr<SolidElement> mElement = nullptr;
    
    // buffer
    eigen::RMatX3_RM mBufferU = eigen::RMatX3_RM(0, 3);
    eigen::RMatX9_RM mBufferG = eigen::RMatX9_RM(0, 9);
    eigen::RMatX6_RM mBufferE = eigen::RMatX6_RM(0, 6);
    eigen::RMatX3_RM mBufferR = eigen::RMatX3_RM(0, 3);
    eigen::RMatX6_RM mBufferS = eigen::RMatX6_RM(0, 6);
    
    
    /////////////////////////// process ///////////////////////////
public:
    // process and report to group
    void processReport(int bufferLine,
                       const channel::solid::ChannelOptions &chops,
                       int stationIndex, eigen::RTensor3 &bufferFields);
    
private:
    // process 1: rotate
    void rotate(int bufferLine, const channel::solid::ChannelOptions &chops);
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // expand workspace for record
    static void
    expandWorkspaceRecord(int nu_1,
                          const channel::solid::ChannelOptions &chops) {
        // for element
        if (chops.mNeedBufferU && sUXN3.rows() < nu_1) {
            sUXN3.resize(nu_1, spectral::nPEM * 3);
            sUX3.resize(nu_1, 3);
        }
        
        if (chops.mNeedBufferG && sGXN9.rows() < nu_1) {
            sGXN9.resize(nu_1, spectral::nPEM * 9);
            sGX9.resize(nu_1, 9);
        }
        
        if (chops.mNeedBufferE && sEXN6.rows() < nu_1) {
            sEXN6.resize(nu_1, spectral::nPEM * 6);
            sEX6.resize(nu_1, 6);
        }
        
        if (chops.mNeedBufferR && sRXN3.rows() < nu_1) {
            sRXN3.resize(nu_1, spectral::nPEM * 3);
            sRX3.resize(nu_1, 3);
        }
        
        if (chops.mNeedBufferS && sSXN6.rows() < nu_1) {
            sSXN6.resize(nu_1, spectral::nPEM * 6);
            sSX6.resize(nu_1, 6);
        }
    }
    
    // workspace for record
    inline static eigen::RRow3 sU3;
    inline static eigen::RRow9 sG9;
    inline static eigen::RRow6 sE6;
    inline static eigen::RRow3 sR3;
    inline static eigen::RRow6 sS6;
    inline static eigen::CMatX3 sUX3 = eigen::CMatX3(0, 3);
    inline static eigen::CMatX9 sGX9 = eigen::CMatX9(0, 9);
    inline static eigen::CMatX6 sEX6 = eigen::CMatX6(0, 6);
    inline static eigen::CMatX3 sRX3 = eigen::CMatX3(0, 3);
    inline static eigen::CMatX6 sSX6 = eigen::CMatX6(0, 6);
    inline static
    eigen::CMatXN3 sUXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static
    eigen::CMatXN9 sGXN9 = eigen::CMatXN9(0, spectral::nPEM * 9);
    inline static
    eigen::CMatXN6 sEXN6 = eigen::CMatXN6(0, spectral::nPEM * 6);
    inline static
    eigen::CMatXN3 sRXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static
    eigen::CMatXN6 sSXN6 = eigen::CMatXN6(0, spectral::nPEM * 6);
};

#endif /* StationSolid_hpp */
