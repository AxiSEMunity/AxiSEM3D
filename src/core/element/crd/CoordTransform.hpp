//
//  CoordTransform.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  coordinate transform between (s,phi,z) and (R,T,Z)
//  for vector and 2nd-order tensor fields

#ifndef CoordTransform_hpp
#define CoordTransform_hpp

#include "eigen_element.hpp"

class CoordTransform {
public:
    // destructor
    virtual ~CoordTransform() = default;
    
    // (s,phi,z) -> (R,T,Z)
    virtual
    void transformSPZ_RTZ3(eigen::vec_ar3_CMatPP_RM &ui, int nu_1) const = 0;
    
    // (R,T,Z) -> (s,phi,z)
    virtual
    void transformRTZ_SPZ3(eigen::vec_ar3_CMatPP_RM &ui, int nu_1) const = 0;
    
    // (s,phi,z) -> (R,T,Z) for nabla
    virtual
    void transformSPZ_RTZ9(eigen::vec_ar9_CMatPP_RM &nij, int nu_1) const = 0;
    
    // (R,T,Z) -> (s,phi,z) for nabla
    virtual
    void transformRTZ_SPZ9(eigen::vec_ar9_CMatPP_RM &nij, int nu_1) const = 0;
    
    // (s,phi,z) -> (R,T,Z) for Voigt
    virtual
    void transformSPZ_RTZ6(eigen::vec_ar6_CMatPP_RM &eij, int nu_1) const = 0;
    
    // (R,T,Z) -> (s,phi,z) for Voigt
    virtual
    void transformRTZ_SPZ6(eigen::vec_ar6_CMatPP_RM &sij, int nu_1) const = 0;
};

#endif /* CoordTransform_hpp */
