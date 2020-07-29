//
//  ElementOpFluid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output in fluid

#ifndef ElementOpFluid_hpp
#define ElementOpFluid_hpp

#include "ElementOp.hpp"
#include "FluidElement.hpp"
#include "Point.hpp"

class ElementOpFluid: public ElementOp {
public:
    /////////////////////////// setup ///////////////////////////
    // set element
    void setElement(const std::shared_ptr<FluidElement> &element) {
        mElement = element;
    }
    
    // set in group
    void setInGroup(int dumpIntv, const channel::fluid::ChannelOptions &chops,
                    int npnts, int nphis);
    
    
    /////////////////////////// get sizes ///////////////////////////
    // get na
    int getNa(int nphis) const {
        if (nphis == 0) {
            return mElement->getNr();
        } else {
            return nphis;
        }
    }
    
    // get nu_1
    int getNu_1() const {
        return mElement->getNu_1();
    }
    
    // get element tag
    int getElementTag() const {
        return mElement->getQuadTag();
    }
    
    // get coords
    eigen::DRowX getCoords(const std::vector<int> &ipnts) const {
        eigen::DRowX sz(ipnts.size() * 2);
        for (int ip = 0; ip < ipnts.size(); ip++) {
            sz.block(0, ip * 2, 1, 2) =
            mElement->getPoint(ipnts[ip]).getCoords();
        }
        return sz;
    }
    
    
    /////////////////////////// record ///////////////////////////
public:
    // record
    void
    record(int bufferLine, const channel::fluid::ChannelOptions &chops,
           const std::vector<int> &ipnts, const eigen::CMatXX &expIAlphaPhi);
    
private:
    // element
    std::shared_ptr<FluidElement> mElement = nullptr;
    
    // buffer
    eigen::RTensor4 mBufferX;
    eigen::RTensor4 mBufferU;
    eigen::RTensor4 mBufferP;
    eigen::RTensor4 mBufferD;
    
    
    /////////////////////////// process ///////////////////////////
public:
    // process and report to group
    void
    processReport(int bufferLine, const channel::fluid::ChannelOptions &chops,
                  int elemIndexNaGrid, int naGridIndex,
                  std::vector<eigen::RTensor5> &ioBuffers);
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // expand workspace for record
    static void
    expandWorkspaceRecord(int nu_1, int na,
                          const channel::fluid::ChannelOptions &chops) {
        // use (na + 1) to handle Nyquist
        int na_1 = na + 1;
        
        if (chops.mNeedBufferX) {
            if (sCXXN1.rows() < nu_1) {
                sCXXN1.resize(nu_1, spectral::nPEM * 1);
                sCXXX.resize(nu_1, spectral::nPEM * 1);
            }
            if (sRXXX.rows() < na_1) {
                sRXXX.resize(na_1, spectral::nPEM * 1);
            }
        }
        if (chops.mNeedBufferU) {
            if (sCUXN3.rows() < nu_1) {
                sCUXN3.resize(nu_1, spectral::nPEM * 3);
                sCUXX.resize(nu_1, spectral::nPEM * 3);
            }
            if (sRUXX.rows() < na_1) {
                sRUXX.resize(na_1, spectral::nPEM * 3);
            }
        }
        if (chops.mNeedBufferP) {
            if (sCPXN1.rows() < nu_1) {
                sCPXN1.resize(nu_1, spectral::nPEM * 1);
                sCPXX.resize(nu_1, spectral::nPEM * 1);
            }
            if (sRPXX.rows() < na_1) {
                sRPXX.resize(na_1, spectral::nPEM * 1);
            }
        }
        if (chops.mNeedBufferD) {
            if (sCDXN1.rows() < nu_1) {
                sCDXN1.resize(nu_1, spectral::nPEM * 1);
                sCDXX.resize(nu_1, spectral::nPEM * 1);
            }
            if (sRDXX.rows() < na_1) {
                sRDXX.resize(na_1, spectral::nPEM * 1);
            }
        }
    }
    
    // workspace for record
    // get response from elememt
    inline static eigen::CMatXN  sCXXN1 = eigen::CMatXN (0, spectral::nPEM * 1);
    inline static eigen::CMatXN3 sCUXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static eigen::CMatXN  sCPXN1 = eigen::CMatXN (0, spectral::nPEM * 1);
    inline static eigen::CMatXN  sCDXN1 = eigen::CMatXN (0, spectral::nPEM * 1);
    
    // in-plane downsampling
    inline static eigen::CMatXX sCXXX = eigen::CMatXX(0, spectral::nPEM * 1);
    inline static eigen::CMatXX sCUXX = eigen::CMatXX(0, spectral::nPEM * 3);
    inline static eigen::CMatXX sCPXX = eigen::CMatXX(0, spectral::nPEM * 1);
    inline static eigen::CMatXX sCDXX = eigen::CMatXX(0, spectral::nPEM * 1);
    
    // making real
    inline static eigen::RMatXX_RM sRXXX = eigen::RMatXX(0, spectral::nPEM * 1);
    inline static eigen::RMatXX_RM sRUXX = eigen::RMatXX(0, spectral::nPEM * 3);
    inline static eigen::RMatXX_RM sRPXX = eigen::RMatXX(0, spectral::nPEM * 1);
    inline static eigen::RMatXX_RM sRDXX = eigen::RMatXX(0, spectral::nPEM * 1);
};

#endif /* ElementOpFluid_hpp */
