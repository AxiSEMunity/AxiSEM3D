//
//  NrField.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#ifndef NrField_hpp
#define NrField_hpp

#include "eigen_generic.hpp"
#include <memory>

// need distance tolerance for NrFieldPointwise
class ExodusMesh;

namespace eigen {
    // coords
    typedef Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> DMatX2_RM;
}

class NrField {
public:
    // destructor
    virtual ~NrField() = default;
    
    // get nr by (s, z)
    virtual eigen::IColX getNrAtPoints(const eigen::DMatX2_RM &sz) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    
    ////////////////////////////// static //////////////////////////////
public:
    // build from inparam
    static std::unique_ptr<const NrField>
    buildInparam(const ExodusMesh &exodusMesh);
    
    // is lucky number
    static bool isLuckyNumber(int n);
    
    // next lucky number
    static int nextLuckyNumber(int n);
    
private:
    // lucky numbers
    static std::vector<int> sLuckyNumbers;
};

#endif /* NrField_hpp */
