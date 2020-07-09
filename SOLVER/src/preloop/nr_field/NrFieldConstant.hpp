//
//  NrFieldConstant.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  constant Nr(s,z)

#ifndef NrFieldConstant_hpp
#define NrFieldConstant_hpp

#include "NrField.hpp"

class NrFieldConstant: public NrField {
public:
    // constructor
    NrFieldConstant(int nr): mNr(nr) {
        // nothing
    }
    
    // get nr by (s, z)
    eigen::IColX getNrAtPoints(const eigen::DMatX2_RM &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // the constant nr value
    const int mNr;
};

#endif /* NrFieldConstant_hpp */
