//
//  eigen_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/11/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  eigen tools

#ifndef eigen_tools_hpp
#define eigen_tools_hpp

#include "eigen.hpp"
#include "bstring.hpp"

namespace eigen_tools {
    /////////////////////////////// memory ///////////////////////////////
    // memory info
    static double useless = 0.;
    template <class EigenMat>
    std::string memoryInfo(const EigenMat &mat, const std::string &title,
                           double &addToMemGB = useless, int count = 1) {
        // type name
        const std::string &tname =
        bstring::typeName<typename EigenMat::Scalar>();
        // scalar size
        size_t scalar = sizeof(typename EigenMat::Scalar);
        // total size
        double memGB = count * (mat.size() * scalar) / 1e9;
        // format
        std::stringstream ss;
        ss << title << ": ";
        if (count > 1) {
            ss << "matrix count = " << count << "; ";
        }
        ss << "dimensions = " << mat.rows() << " x " << mat.cols() << "; ";
        ss << "scalar type = " << tname << " (" << scalar << " bytes); ";
        ss << "memory = " << memGB << " GB";
        // add to
        addToMemGB += memGB;
        return ss.str();
    }
    
    // memory info for tensor
    template <class EigenTensor>
    std::string memoryInfoTensor(const EigenTensor &tensor,
                                 const std::string &title,
                                 double &addToMemGB = useless, int count = 1) {
        // type name
        const std::string &tname =
        bstring::typeName<typename EigenTensor::Scalar>();
        // scalar size
        size_t scalar = sizeof(typename EigenTensor::Scalar);
        // total size
        double memGB = count * (tensor.size() * scalar) / 1e9;
        // format
        std::stringstream ss;
        ss << title << ": ";
        if (count > 1) {
            ss << "tensor count = " << count << "; ";
        }
        ss << "dimensions = ";
        for (int ir = 0; ir < tensor.dimensions().size() - 1; ir++) {
            ss << tensor.dimension(ir) << " x ";
        }
        ss << tensor.dimension(tensor.dimensions().size() - 1) << "; ";
        ss << "scalar type = " << tname << " (" << scalar << " bytes); ";
        ss << "memory = " << memGB << " GB";
        // add to
        addToMemGB += memGB;
        return ss.str();
    }
    
    
    /////////////////////////////// Fourier ///////////////////////////////
    // compute 2 * exp(I * alpha * phi)
    template <typename CColX, typename T = typename CColX::Scalar::value_type>
    void computeTwoExpIAlphaPhi(int nu_1, double phi, CColX &twoExpIAlphaPhi) {
        if (nu_1 <= 0) {
            // nothing
            return;
        }
        // no factor two on order zero
        twoExpIAlphaPhi(0) = 1.;
        // non-zero orders
        static const T two = 2.;
        std::complex<T> i_phi(0., phi);
        for (int alpha = 1; alpha < nu_1; alpha++) {
            twoExpIAlphaPhi(alpha) = two * exp((T)alpha * i_phi);
        }
    }
    
    // compute real value of a Fourier series at phi with precomputed exp
    template <typename RRowM, typename CColX, typename CMatXM>
    void computeFourierAtPhiExp(const CMatXM &coeffs, int nu_1,
                                const CColX &twoExpIAlphaPhi, RRowM &reals,
                                int realRow = 0, int ncol = -1) {
        if (ncol < 0) {
            ncol = (int)coeffs.cols();
        }
        // sum over Fourier orders
        reals.block(realRow, 0, 1, ncol) =
        (twoExpIAlphaPhi.topRows(nu_1).asDiagonal() *
         coeffs.block(0, 0, nu_1, ncol)).real().colwise().sum();
    }
    
    // compute real value of a Fourier series at phi
    template <typename RRowM, typename CColX, typename CMatXM>
    void computeFourierAtPhi(const CMatXM &coeffs, int nu_1, double phi,
                             RRowM &reals, int realRow = 0, int ncol = -1) {
        // 2 * exp(I * alpha * phi)
        CColX twoExpIAlphaPhi(nu_1);
        computeTwoExpIAlphaPhi(nu_1, phi, twoExpIAlphaPhi);
        // sum over Fourier orders
        computeFourierAtPhiExp(coeffs, nu_1, twoExpIAlphaPhi, reals,
                               realRow, ncol);
    }
}

#endif /* eigen_tools_hpp */
