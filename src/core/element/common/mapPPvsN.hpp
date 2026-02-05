//
//  mapPPvsN.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/9/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  conversion between vec_arD_MatPP_RM and MatXND
//  where D = dimension of data

#ifndef mapPPvsN_hpp
#define mapPPvsN_hpp

#include "eigen.hpp"
#include "spectral.hpp"

namespace mapPPvsN {
    using spectral::nPEM;
    
    // PP -> N
    template <typename vec_arD_MatPP_RM, typename MatXND,
    typename RowN = Eigen::Matrix<typename MatXND::Scalar, 1, nPEM>,
    int D = MatXND::ColsAtCompileTime / nPEM>
    void PP2N(const vec_arD_MatPP_RM &inPP, MatXND &outN, int nu_1) {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            for (int idim = 0; idim < D; idim++) {
                outN.block(alpha, idim * nPEM, 1, nPEM) =
                Eigen::Map<const RowN>(inPP[alpha][idim].data());
            }
        }
    }
    
    // N -> PP
    template <typename vec_arD_MatPP_RM, typename MatXND,
    typename RowN = Eigen::Matrix<typename MatXND::Scalar, 1, nPEM>,
    int D = MatXND::ColsAtCompileTime / nPEM>
    void N2PP(const MatXND &inN, vec_arD_MatPP_RM &outPP, int nu_1) {
        for (int alpha = 0; alpha < nu_1; alpha++) {
            for (int idim = 0; idim < D; idim++) {
                Eigen::Map<RowN>(outPP[alpha][idim].data()) =
                inN.block(alpha, idim * nPEM, 1, nPEM);
            }
        }
    }
}

#endif /* mapPPvsN_hpp */
