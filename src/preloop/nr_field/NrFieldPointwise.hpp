//
//  NrFieldPointwise.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  pointwise Nr(s,z)

#ifndef NrFieldPointwise_hpp
#define NrFieldPointwise_hpp

#include "NrField.hpp"
#include "RTreeND.hpp"

class NrFieldPointwise: public NrField {
public:
    // constructor
    NrFieldPointwise(const std::string &fname, double factor,
                     double distTolExact);
    
    // get nr by (s, z)
    eigen::IColX getNrAtPoints(const eigen::DMatX2_RM &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // file name
    const std::string mFilename;
    
    // factor
    const double mFactor;
    
    // dist tolerance for exact match
    const double mDistTolExact;
    
    // rtree
    std::unique_ptr<const RTreeND<2, 1, int>> mRTree = nullptr;
    
    // scanning
    long mSumNrStart = -1;
};

#endif /* NrFieldPointwise_hpp */
