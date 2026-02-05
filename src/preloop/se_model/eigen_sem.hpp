//
//  eigen_sem.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/25/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  eigen typedef for SEM

#ifndef eigen_sem_hpp
#define eigen_sem_hpp

#include "eigen_generic.hpp"
#include "spectral.hpp"

namespace eigen {
    using Eigen::Dynamic;
    using Eigen::RowMajor;
    using spectral::nPEM;
    using spectral::nPED;
    
    ///////// nodal /////////
    // 2D coords
    typedef Eigen::Matrix<double, 2, 1> DCol2;
    typedef Eigen::Matrix<double, 1, 2> DRow2;
    // Jacobian
    typedef Eigen::Matrix<double, 2, 2> DMat22;
    // nodal field
    typedef Eigen::Matrix<double, 1, 4> DRow4;
    // nodal coords
    typedef Eigen::Matrix<double, 2, 4> DMat24;
    
    ///////// GLL coordinates /////////
    // 2D coordinates on an element
    typedef Eigen::Matrix<double, 2, nPEM> DMat2N;
    // 3D coordinates flattened
    typedef Eigen::Matrix<double, Dynamic, 3> DMatX3;
    
    ///////// GLL fields /////////
    // 1D scalar field on an element
    typedef Eigen::Matrix<int, 1, nPEM> IRowN;
    typedef Eigen::Matrix<double, 1, nPEM> DRowN;
    // 3D scalar field on an element with max Nr
    typedef Eigen::Matrix<double, Dynamic, nPEM> DMatXN;
    // 3D scalar field on an element with different Nr's
    typedef std::array<eigen::DColX, nPEM> arN_DColX;
    typedef std::array<eigen::IColX, nPEM> arN_IColX;
    // 3D scalar field on an edge with different Nr's
    typedef std::array<eigen::DColX, nPED> arP_DColX;
    // 3D vector field on an edge with different Nr's
    typedef std::array<eigen::DMatX3, nPED> arP_DMatX3;
    
    ///////// GLL structured /////////
    typedef Eigen::Matrix<double, nPED, 1> DColP;
    typedef Eigen::Matrix<double, nPED, nPED, Eigen::RowMajor> DMatPP_RM;
}


/////////////////////////// 1D/3D operations ///////////////////////////
namespace op1D_3D {
    // cwiseProduct
    template <class MatA, class MatB, class MatR>
    void times(const MatA &A, const MatB &B, MatR &R) {
        if (A.rows() == B.rows()) {
            // 1Dx1D or 3Dx3D
            R = A.cwiseProduct(B);
        } else if (A.rows() == 1) {
            // 1Dx3D
            R = (B * A.asDiagonal()).eval(); // gcc needs eval()
        } else if (B.rows() == 1) {
            // 3Dx1D
            R = (A * B.asDiagonal()).eval(); // gcc needs eval()
        } else {
            throw std::runtime_error("op1D_3D::times || "
                                     "Incompatible dimemsions.");
        }
    }
    
    // addTo
    template <class MatIn, class MatMe>
    void addTo(const MatIn &in, MatMe &me) {
        // input can be zero-sized with mpi
        if (in.rows() == 0) {
            return;
        }
        
        // first visit
        if (me.rows() == 0) {
            me = in;
            return;
        }
        
        // add to
        if (me.rows() == in.rows()) {
            // add 1D to 1D or add 3D to 3D
            me += in;
        } else if (me.rows() == 1) {
            // add 3D to 1D
            me = (in.rowwise() + me.row(0)).eval(); // gcc needs eval()
        } else if (in.rows() == 1) {
            // add 1D to 3D
            me.rowwise() += in.row(0);
        } else {
            throw std::runtime_error("op1D_3D::addTo || "
                                     "Incompatible dimemsions.");
        }
    }
    
    // regularize 1D or 3D
    template <class Mat>
    bool regularize1D(const std::vector<std::reference_wrapper<Mat>> &mats) {
        // find max nr
        int nrMax = -1;
        for (const std::reference_wrapper<Mat> &mat: mats) {
            int nr = (int)mat.get().rows();
            if (nr == 0) {
                throw std::runtime_error("op1D_3D::regularize1D || "
                                         "Uninitialized matrix.");
            }
            nrMax = std::max(nr, nrMax);
        }
        
        // 1D
        if (nrMax == 1) {
            return true;
        }
        
        // all extend to 3D
        for (const std::reference_wrapper<Mat> &mat: mats) {
            int nr = (int)mat.get().rows();
            if (nr == 1) {
                // extend 1D to 3D
                // gcc needs eval()
                mat.get() = mat.get().replicate(nrMax, 1).eval();
            } else if (nr == nrMax) {
                // nothing
            } else {
                throw std::runtime_error("op1D_3D::regularize1D || "
                                         "Incompatible dimemsions.");
            }
        }
        return false;
    }
    
    // extend 1D to 3D
    template <class Mat>
    Mat to3D(const Mat &m, int nr) {
        if (m.rows() == nr) {
            // already 3D, do nothing
            return m;
        }
        // extend ONLY 1D to 3D
        if (m.rows() == 1) {
            return m.replicate(nr, 1);
        } else {
            throw std::runtime_error("op1D_3D::to3D || "
                                     "Incompatible dimemsions.");
        }
    }
    
    // reduce 3D to 1D if possible
    template <class Mat, typename T = typename Mat::Scalar>
    void tryReduceTo1D(Mat &m) {
        // null or 1D
        if (m.rows() == 0 || m.rows() == 1) {
            return;
        }
        // reduce to 1D if all rows are the same
        if ((m.rowwise() - m.row(0)).norm() <
            m.norm() * numerical::epsilon<T>()) {
            m.conservativeResize(1, m.cols());
        }
    }
    
    // 1D structured
    inline eigen::DMatPP_RM toPP(const eigen::DMatXN &mat) {
        return Eigen::Map<const eigen::DMatPP_RM>(mat.data());
    }
}

#endif /* eigen_sem_hpp */
