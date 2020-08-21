//
//  ElementOp.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/27/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  element output

#ifndef ElementOp_hpp
#define ElementOp_hpp

#include "eigen_element_op.hpp"
#include "eigen_tools.hpp"
#include "channel.hpp"

class ElementOp {
public:
    // constructor
    ElementOp(const std::vector<int> &ipnts): mIPnts(ipnts) {
        // nothing
    }
    
    // destructor
    virtual ~ElementOp() = default;
    
    ///////////////////////// template functions /////////////////////////
    // record: inplane downsampling and making real
    template <int D,
    typename CMatXND =
    Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, spectral::nPEM * D>,
    typename RMatXND_RM =
    Eigen::Matrix<numerical::Real, Eigen::Dynamic, spectral::nPEM * D,
    Eigen::RowMajor>>
    void recordToElem(CMatXND &cxnd, int nu_1, RMatXND_RM &rxnd,
                      const eigen::CMatXX &expIAlphaPhi,
                      eigen::RTensor4 &field, int bufferLine) const {
        // inplane downsampling
        int npnts = (int)mIPnts.size();
        if (npnts != spectral::nPEM) {
            for (int dim = 0; dim < D; dim++) {
                cxnd.block(0, npnts * dim, nu_1, npnts) =
                cxnd.block(0, spectral::nPEM * dim, nu_1, spectral::nPEM)
                (Eigen::all, mIPnts);
            }
        }
        
        // making real
        int nphis = (int)expIAlphaPhi.cols();
        int npntsD = npnts * D;
        if (nphis == 0) {
            // no Fourier interpolation, just reform complex to real
            // zeroth
            rxnd.block(0, 0, 1, npntsD) = cxnd.block(0, 0, 1, npntsD).real();
            // higher (Nyquist truncated later)
            for (int alpha = 1; alpha < nu_1; alpha++) {
                rxnd.block(alpha * 2 - 1, 0, 1, npntsD) =
                cxnd.block(alpha, 0, 1, npntsD).real();
                rxnd.block(alpha * 2 - 0, 0, 1, npntsD) =
                cxnd.block(alpha, 0, 1, npntsD).imag();
            }
        } else {
            // Fourier interpolation
            for (int iphi = 0; iphi < nphis; iphi++) {
                eigen_tools::
                computeFourierAtPhiExp(cxnd, nu_1, expIAlphaPhi.col(iphi),
                                       rxnd, iphi, npntsD);
            }
        }
        
        // write to buffer (Nyquist truncated here)
        static const eigen::IArray3 shuffle = {0, 2, 1};
        static const eigen::IArray3 zero = {0, 0, 0};
        int na = (int)field.dimension(0);
        field.slice(eigen::IArray4({0, 0, 0, bufferLine}),
                    eigen::IArray4({na, npnts, D, 1})).
        reshape(eigen::IArray3{na, npnts, D})
        = // the longest C++ statement I have written
        Eigen::TensorMap<eigen::RTensor3>
        (rxnd.data(), eigen::IArray3({rxnd.rows(), rxnd.cols(), 1})).
        slice(zero, eigen::IArray3({na, npntsD, 1})).
        reshape(eigen::IArray3{na, D, npnts}).shuffle(shuffle);
    }
    
    // dump: compute channel and feed IO buffer
    template <int D> static
    void dumpToIO(const eigen::RTensor4 &field, int fieldIndex,
                  int bufferLine, int channelIndex, int elemIndexNaGrid,
                  int naGridIndex, std::vector<eigen::RTensor5> &ioBuffers) {
        // result reference
        int na = (int)field.dimension(0);
        int npnts = (int)field.dimension(1);
        eigen::IArray5 loc5 = {elemIndexNaGrid, 0, 0, channelIndex, 0};
        eigen::IArray5 len5 = {1, na, npnts, 1, bufferLine};
        eigen::IArray3 copy = {na, npnts, bufferLine};
        auto res = ioBuffers[naGridIndex].slice(loc5, len5).reshape(copy);
        
        // indexing of input elemBuffer
        eigen::IArray4 loc4 = {0, 0, 0, 0};
        eigen::IArray4 len4 = {na, npnts, 1, bufferLine};
        
        // channel index
        static const int chrank = 2;
        static const eigen::IArray1 chDim = {2};
        
        // channel
        if (fieldIndex >= 0) {
            loc4[chrank] = fieldIndex;
            len4[chrank] = 1;
            res = field.slice(loc4, len4).reshape(copy);
        } else {
            if (D == 3) {
                if (fieldIndex == -1) {
                    res.setZero();
                    len4[chrank] = 1;
                    for (int dim = 0; dim < 3; dim++) {
                        loc4[chrank] = dim;
                        res += field.slice(loc4, len4).reshape(copy).square();
                    }
                    res = res.sqrt();
                } else {
                    throw std::runtime_error("ElementOp::dumpToIO || "
                                             "Invalid channel setting.");
                }
            } else if (D == 6) {
                if (fieldIndex == -1) {
                    // trace
                    loc4[chrank] = 0;
                    len4[chrank] = 3;
                    res = field.slice(loc4, len4).sum(chDim);
                } else if (fieldIndex == -2) {
                    // J2 = I1 ^ 2 / 3 - I2
                    // I1
                    loc4[chrank] = 0;
                    len4[chrank] = 3;
                    res = (field.slice(loc4, len4).sum(chDim).
                           square() * numerical::Real(1./3));
                    // I2: permutation 0, 1, 2
                    len4[chrank] = 1;
                    eigen::IArray4 lox4 = loc4;
                    for (int dim = 0; dim < 3; dim++) {
                        loc4[chrank] = dim;
                        lox4[chrank] = (dim == 2) ? 0 : (dim + 1);
                        res -= (field.slice(loc4, len4).reshape(copy) *
                                field.slice(lox4, len4).reshape(copy));
                    }
                    // I2: permutation 3, 4, 5
                    for (int dim = 3; dim < 6; dim++) {
                        loc4[chrank] = dim;
                        lox4[chrank] = (dim == 5) ? 3 : (dim + 1);
                        res += (field.slice(loc4, len4).reshape(copy) *
                                field.slice(lox4, len4).reshape(copy));
                    }
                } else {
                    throw std::runtime_error("ElementOp::dumpToIO || "
                                             "Invalid channel setting.");
                }
            } else if (D == 9) {
                if (fieldIndex == -1) {
                    // trace
                    loc4[chrank] = 0;
                    len4[chrank] = 3;
                    res = field.slice(loc4, len4).sum(chDim);
                } else {
                    throw std::runtime_error("ElementOp::dumpToIO || "
                                             "Invalid channel setting.");
                }
            } else {
                throw std::runtime_error("ElementOp::dumpToIO || "
                                         "Invalid channel setting.");
            }
        }
    }
    
    ////////// data //////////
protected:
    // ipnts are different for elements if edge is specified
    const std::vector<int> mIPnts;
};

#endif /* ElementOp_hpp */
