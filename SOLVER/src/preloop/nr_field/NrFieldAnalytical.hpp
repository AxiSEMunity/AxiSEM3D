//
//  NrFieldAnalytical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  analytical Nr(s,z)

#ifndef NrFieldAnalytical_hpp
#define NrFieldAnalytical_hpp

#include "NrField.hpp"
#include <vector>

class NrFieldAnalytical: public NrField {
public:
    // constructor
    NrFieldAnalytical();
    
    // get nr by (s, z)
    eigen::IColX getNrAtPoints(const eigen::DMatX2_RM &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // code ID
    static std::string sCodeID;
    
    // TODO: add your data here
    // below are data for
    // sCodeID = "depth-dependent (AxiSEM3D default)"
    std::vector<double> mControlDepths;
    std::vector<double> mControlNrs;
};

#endif /* NrFieldAnalytical_hpp */
