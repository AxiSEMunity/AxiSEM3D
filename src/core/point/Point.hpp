//
//  Point.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/24/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  GLL point

#ifndef Point_hpp
#define Point_hpp

#include "eigen_point.hpp"
#include "Mass.hpp"
#include "bstring.hpp"

class Point {
public:
    // constructor
    Point(int nr, const eigen::DRow2 &crds, int meshTag,
          std::unique_ptr<const Mass> &mass):
    mNr(nr), mNu_1(mNr / 2 + 1), mCoords(crds),
    mMeshTag(meshTag), mMass(mass.release()) {
        // nothing
    }
    
    // destructor
    virtual ~Point() = default;
    
    // type info
    std::string typeInfo() const {
        return bstring::typeName(*this) + "$" + bstring::typeName(*mMass);
    }
    
    // cost signature
    std::string costSignature() const {
        std::stringstream ss;
        ss << typeInfo() << "$" << mNr;
        return ss.str();
    }
    
    
    /////////////////////////// properties ///////////////////////////
    // nr
    int getNr() const {
        return mNr;
    }
    
    // nu
    int getNu_1() const {
        return mNu_1;
    }
    
    // location
    const eigen::DRow2 &getCoords() const {
        return mCoords;
    }
    
    // tag
    int getMeshTag() const {
        return mMeshTag;
    }
    
    
    /////////////////////////// domain ///////////////////////////
    // set domain tag
    void setDomainTag(int domainTag) {
        mDomainTag = domainTag;
    }
    
    // get domain tag
    int getDomainTag() const {
        return mDomainTag;
    }
    
    
    /////////////////////// wavefield scanning ///////////////////////
    // enable scanning
    // need this abstract for vertex-only option
    virtual void enableScanning() = 0;
    
protected:
    // order
    const int mNr;
    const int mNu_1;
    
    // location (s, z)
    const eigen::DRow2 mCoords;
    
    // tag in spectral-element mesh
    // * determined by location, mpi-independent
    // * a pair of solid and fluid points at the same location
    //   have the same mesh tag
    const int mMeshTag;
    
    // mass
    const std::unique_ptr<const Mass> mMass;
    
    // tag (index) in computational domain, mpi-dependent
    int mDomainTag = -1;
};


/////////////////////// wavefield scanning ///////////////////////
class Scanning1D {
    // typedef
    typedef numerical::Real Real;
    // use H2 as the key for fast insertion and deletion
    typedef std::map<Real, int> ScanMap;
    
public:
    // do scanning
    template <class CColX>
    void doScanning(Real relTolFourierH2, Real relTolH2, Real absTolH2,
                    int maxNumPeaks, const CColX &fseries) {
        // l2 and h2 norm
        Real L2 = fseries.squaredNorm();
        Real H2 = L2 - .5 * fseries.row(0).squaredNorm();
        
        // current max
        Real maxH2 = 0.;
        if (mCandidates.size() > 0) {
            maxH2 = mCandidates.rbegin()->first;
        }
        
        // determine whether this can be a candidate
        if (!(H2 > absTolH2 && H2 > relTolH2 * maxH2 &&
              mH2Prev > H2 && mH2Prev > mH2PrevPrev)) {
            // disqualified, update H2 series
            mH2PrevPrev = mH2Prev;
            mH2Prev = H2;
            return;
        }
        
        // test a smaller Nu, starting from 0 (most aggressive)
        Real absFourier = H2 * relTolFourierH2;
        int newNu_1 = 1;
        for (; newNu_1 < fseries.size(); newNu_1++) {
            // no need to differentiate L2 and H2 for difference
            if (L2 - fseries.topRows(newNu_1).squaredNorm() < absFourier) {
                // if this never happens, newNu_1 will become fseries.size()
                break;
            }
        }
        
        // add candidate
        auto itFind = mCandidates.find(H2);
        if (itFind != mCandidates.end()) {
            // replace
            itFind->second = std::max(itFind->second, newNu_1);
        } else {
            // remove tiny values due to new maximum
            if (maxH2 < H2) {
                // find deletion point
                auto deleteIter = mCandidates.lower_bound(relTolH2 * H2);
                // delete from begin to deletion point
                mCandidates.erase(mCandidates.cbegin(), deleteIter);
            }
            // insert
            mCandidates.insert({H2, newNu_1});
        }
        
        // check number of peaks
        if (mCandidates.size() > maxNumPeaks) {
            // remove the first (smallest H2)
            mCandidates.erase(mCandidates.cbegin());
        }
        
        // update H2 series
        mH2PrevPrev = mH2Prev;
        mH2Prev = H2;
    }
    
    // report scanning Nr
    int reportScanningNr(int originalNr) const {
        if (mCandidates.size() == 0) {
            // use 1 for empty (wave has not arrived at this point)
            return 1;
        } else {
            // maximum value in map
            int nu_1 = std::max_element(mCandidates.begin(), mCandidates.end(),
                                        [](const ScanMap::value_type &a,
                                           const ScanMap::value_type &b) {
                return a.second < b.second;})->second;
            // return nr
            return std::min(nu_1 * 2 + 1, originalNr);
        }
    }
    
private:
    // previous H2 values for peak detection
    Real mH2PrevPrev = 0.;
    Real mH2Prev = 0.;
    
    // candidates
    ScanMap mCandidates;
};

#endif /* Point_hpp */
