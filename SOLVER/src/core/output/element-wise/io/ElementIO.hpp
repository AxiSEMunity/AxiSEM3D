//
//  ElementIO.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/26/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  IO for element output

#ifndef ElementIO_hpp
#define ElementIO_hpp

#include "eigen_element_op.hpp"

class ElementIO {
public:
    // destructor
    virtual ~ElementIO() = default;
    
    // initialize
    virtual void initialize(const std::string &groupName,
                            int numRecordSteps,
                            const std::vector<std::string> &channels,
                            int npnts, const std::vector<int> &naGrid,
                            const eigen::IMatX4_RM &elemNaInfo,
                            const eigen::DMatXX_RM &elemCoords);
    
    // finalize
    virtual void finalize() = 0;
    
    // dump to file
    virtual void
    dumpToFile(const eigen::DColX &bufferTime,
               const std::vector<eigen::RTensor5> &bufferFields,
               int bufferLine) = 0;
    
    // set flush
    void setFlush(bool flush) {
        mFlush = flush;
    }
    
    //////////////// mpi global properties ////////////////
protected:
    // total number of elements
    int mNumElementsGlobal = 0;
    
    // mpi rank with maximum number of elements
    int mRankWithMaxNumElements = -1;
    
    // flush
    bool mFlush = false;
};

#endif /* ElementIO_hpp */
