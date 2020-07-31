//
//  ElementOpSolid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output in solid

#ifndef ElementOpSolid_hpp
#define ElementOpSolid_hpp

#include "ElementOp.hpp"
#include "SolidElement.hpp"
#include "Point.hpp"

class ElementOpSolid: public ElementOp {
public:
    // constructor
    ElementOpSolid(const std::vector<int> &ipnts): ElementOp(ipnts) {
        // nothing
    }
    
    
    /////////////////////////// setup ///////////////////////////
    // set element
    void setElement(const std::shared_ptr<SolidElement> &element) {
        mElement = element;
    }
    
    // set in group
    void setInGroup(int dumpIntv, const channel::solid::ChannelOptions &chops,
                    int nphis);
    
    
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
    eigen::DRowX getCoords() const {
        eigen::DRowX sz(mIPnts.size() * 2);
        for (int ip = 0; ip < mIPnts.size(); ip++) {
            sz.block(0, ip * 2, 1, 2) =
            mElement->getPoint(mIPnts[ip]).getCoords();
        }
        return sz;
    }
    
    
    /////////////////////////// record ///////////////////////////
public:
    // record
    void
    record(int bufferLine, const channel::solid::ChannelOptions &chops,
           const eigen::CMatXX &expIAlphaPhi);
    
private:
    // element
    std::shared_ptr<SolidElement> mElement = nullptr;
    
    // buffer
    eigen::RTensor4 mBufferU;
    eigen::RTensor4 mBufferG;
    eigen::RTensor4 mBufferE;
    eigen::RTensor4 mBufferR;
    eigen::RTensor4 mBufferS;
    
    
    /////////////////////////// process ///////////////////////////
public:
    // process and report to group
    void
    processReport(int bufferLine, const channel::solid::ChannelOptions &chops,
                  int elemIndexNaGrid, int naGridIndex,
                  std::vector<eigen::RTensor5> &ioBuffers);
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
private:
    // expand workspace for record
    static void
    expandWorkspaceRecord(int nu_1, int na,
                          const channel::solid::ChannelOptions &chops) {
        // use (na + 1) to handle Nyquist
        int na_1 = na + 1;
        
        if (chops.mNeedBufferU) {
            if (sCUXN3.rows() < nu_1) {
                sCUXN3.resize(nu_1, spectral::nPEM * 3);
                sCUXX.resize(nu_1, spectral::nPEM * 3);
            }
            if (sRUXX.rows() < na_1) {
                sRUXX.resize(na_1, spectral::nPEM * 3);
            }
        }
        if (chops.mNeedBufferG) {
            if (sCGXN9.rows() < nu_1) {
                sCGXN9.resize(nu_1, spectral::nPEM * 9);
                sCGXX.resize(nu_1, spectral::nPEM * 9);
            }
            if (sRGXX.rows() < na_1) {
                sRGXX.resize(na_1, spectral::nPEM * 9);
            }
        }
        if (chops.mNeedBufferE) {
            if (sCEXN6.rows() < nu_1) {
                sCEXN6.resize(nu_1, spectral::nPEM * 6);
                sCEXX.resize(nu_1, spectral::nPEM * 6);
            }
            if (sREXX.rows() < na_1) {
                sREXX.resize(na_1, spectral::nPEM * 6);
            }
        }
        if (chops.mNeedBufferR) {
            if (sCRXN3.rows() < nu_1) {
                sCRXN3.resize(nu_1, spectral::nPEM * 3);
                sCRXX.resize(nu_1, spectral::nPEM * 3);
            }
            if (sRRXX.rows() < na_1) {
                sRRXX.resize(na_1, spectral::nPEM * 3);
            }
        }
        if (chops.mNeedBufferS) {
            if (sCSXN6.rows() < nu_1) {
                sCSXN6.resize(nu_1, spectral::nPEM * 6);
                sCSXX.resize(nu_1, spectral::nPEM * 6);
            }
            if (sRSXX.rows() < na_1) {
                sRSXX.resize(na_1, spectral::nPEM * 6);
            }
        }
    }
    
    // workspace for record
    // get response from elememt
    inline static eigen::CMatXN3 sCUXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static eigen::CMatXN9 sCGXN9 = eigen::CMatXN9(0, spectral::nPEM * 9);
    inline static eigen::CMatXN6 sCEXN6 = eigen::CMatXN6(0, spectral::nPEM * 6);
    inline static eigen::CMatXN3 sCRXN3 = eigen::CMatXN3(0, spectral::nPEM * 3);
    inline static eigen::CMatXN6 sCSXN6 = eigen::CMatXN6(0, spectral::nPEM * 6);
    
    // in-plane downsampling
    inline static eigen::CMatXX sCUXX = eigen::CMatXX(0, spectral::nPEM * 3);
    inline static eigen::CMatXX sCGXX = eigen::CMatXX(0, spectral::nPEM * 9);
    inline static eigen::CMatXX sCEXX = eigen::CMatXX(0, spectral::nPEM * 6);
    inline static eigen::CMatXX sCRXX = eigen::CMatXX(0, spectral::nPEM * 3);
    inline static eigen::CMatXX sCSXX = eigen::CMatXX(0, spectral::nPEM * 6);
    
    // making real
    inline static eigen::RMatXX_RM sRUXX = eigen::RMatXX(0, spectral::nPEM * 3);
    inline static eigen::RMatXX_RM sRGXX = eigen::RMatXX(0, spectral::nPEM * 9);
    inline static eigen::RMatXX_RM sREXX = eigen::RMatXX(0, spectral::nPEM * 6);
    inline static eigen::RMatXX_RM sRRXX = eigen::RMatXX(0, spectral::nPEM * 3);
    inline static eigen::RMatXX_RM sRSXX = eigen::RMatXX(0, spectral::nPEM * 6);
};

#endif /* ElementOpSolid_hpp */
