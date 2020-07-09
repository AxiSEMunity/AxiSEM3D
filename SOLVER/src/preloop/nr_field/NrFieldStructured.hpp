//
//  NrFieldStructured.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  structured Nr(s,z)

#ifndef NrFieldStructured_hpp
#define NrFieldStructured_hpp

#include "NrField.hpp"
#include "StructuredGrid.hpp"

class NrFieldStructured: public NrField {
public:
    // constructor
    NrFieldStructured(const std::string &fname, int valOutOfRange);
    
    // get nr by (s, z)
    eigen::IColX getNrAtPoints(const eigen::DMatX2_RM &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // file name
    const std::string mFilename;
    
    // factor
    const int mValueOutOfRange;
    
    // grid
    std::unique_ptr<const StructuredGrid<2, int>> mGrid = nullptr;
};

#endif /* NrFieldStructured_hpp */
