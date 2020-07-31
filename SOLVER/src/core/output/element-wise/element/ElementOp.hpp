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
    // record: in-plane downsampling and making real
    template <int D, typename CMatXND =
    Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, spectral::nPEM * D>>
    void recordToElem(const CMatXND &cxnd, int nu_1, eigen::CMatXX &cxad,
                      eigen::RMatXX_RM &rxad, const eigen::CMatXX &expIAlphaPhi,
                      eigen::RTensor4 &field, int bufferLine) const {
        // in-plane downsampling
        int npnts = (int)mIPnts.size();
        if (npnts == spectral::nPEM) {
            // no downsampling
            cxad.topRows(nu_1) = cxnd.topRows(nu_1);
        } else {
            // downsampling for each dimension
            for (int dim = 0; dim < D; dim++) {
                cxad.block(0, npnts * dim, nu_1, npnts) =
                cxnd.block(0, spectral::nPEM * dim, nu_1, spectral::nPEM)
                (Eigen::all, mIPnts);
            }
        }
        
        // making real
        int nphis = (int)expIAlphaPhi.cols();
        if (nphis == 0) {
            // no Fourier interpolation, just reform complex to real
            // zeroth
            rxad.block(0, 0, 1, npnts) = cxad.block(0, 0, 1, npnts).real();
            // higher (Nyquist truncated later)
            for (int alpha = 1; alpha < nu_1; alpha++) {
                rxad.block(alpha * 2 - 1, 0, 1, npnts) =
                cxad.block(alpha, 0, 1, npnts).real();
                rxad.block(alpha * 2 - 0, 0, 1, npnts) =
                cxad.block(alpha, 0, 1, npnts).imag();
            }
        } else {
            // Fourier interpolation
            for (int iphi = 0; iphi < nphis; iphi++) {
                eigen_tools::
                computeFourierAtPhiExp(cxad, nu_1, expIAlphaPhi.col(iphi),
                                       rxad, iphi, npnts);
            }
        }
        
        // write to buffer (Nyquist truncated here)
        const static eigen::IArray3 shuffle = {0, 2, 1};
        int na = (int)field.dimension(0);
        eigen::IArray3 copy = {na, npnts, D};
        field.slice(eigen::IArray4({0, 0, 0, bufferLine}),
                    eigen::IArray4({na, npnts, D, 1})).reshape(copy) =
        Eigen::TensorMap<eigen::RTensor3>
        (rxad.data(), eigen::IArray3({rxad.rows(), D, rxad.cols() / D})).
        shuffle(shuffle).slice(eigen::IArray3({0, 0, 0}), copy);
    };
    
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
        eigen::IArray4 copy = {na, npnts, bufferLine};
        auto res = ioBuffers[naGridIndex].slice(loc5, len5).reshape(copy);
        
        // indexing of input elemBuffer
        eigen::IArray4 loc4 = {0, 0, 0, 0};
        eigen::IArray4 len4 = {na, npnts, 1, bufferLine};
        
        // channel
        const int chrank = 2;
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
                    res = field.slice(loc4, len4).sum(eigen::IArray1({2}));
                } else if (fieldIndex == -2) {
                    // J2 = I1 ^ 2 / 3 - I2
                    // I1
                    loc4[chrank] = 0;
                    len4[chrank] = 3;
                    res = (field.slice(loc4, len4).sum(eigen::IArray1({2})).
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
                    res = field.slice(loc4, len4).sum(eigen::IArray1({2}));
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
