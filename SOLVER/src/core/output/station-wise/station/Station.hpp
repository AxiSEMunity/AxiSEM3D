//
//  Station.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/9/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  station

#ifndef Station_hpp
#define Station_hpp

#include "eigen_station.hpp"
#include "eigen_tools.hpp"
#include "channel.hpp"

// short alias for column block operation
#define bblk(dim) block(0, dim, bufferLine, 1)

class Station {
public:
    // constructor
    Station(const std::string &key, double phi,
            double theta, double backAzimuth):
    mKey(key), mPhi(phi), mTheta(theta),
    mBAzPlus90(backAzimuth + numerical::dPi * 0.5) {
        // nothing
    }
    
    // destructor
    virtual ~Station() = default;
    
    // get key
    const std::string &getKey() const {
        return mKey;
    }
    
protected:
    // set element: inplane weights and nu
    void setElement(const eigen::DRowN &weights, int nu_1);
    
    
    ///////////////////////// template functions /////////////////////////
    // record: inplane and Fourier interpolation
    template <int D,
    typename CMatXND =
    Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, spectral::nPEM * D>,
    typename CMatXD =
    Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, D>,
    typename RRowD =
    Eigen::Matrix<numerical::Real, 1, D>>
    void interpolate(const CMatXND &xnd, CMatXD &xd, RRowD &d, int nu_1) const {
        // inplane interpolation
        for (int idim = 0; idim < D; idim++) {
            // Eigen 3.4 feature
            xd.block(0, idim, nu_1, 1) =
            (xnd.block(0, idim * spectral::nPEM, nu_1, spectral::nPEM)
             (Eigen::all, mNonZeroIndices) *
             mNonZeroWeights.asDiagonal()).rowwise().sum();
        }
        
        // azimuthal interpolation
#ifndef _SAVE_MEMORY
        // sum exp
        eigen_tools::computeFourierAtPhiExp(xd, nu_1, m2ExpIAlphaPhi, d);
#else
        // compute exp on-the-fly
        eigen_tools::computeTwoExpIAlphaPhi(nu_1, mPhi, s2ExpIAlphaPhi);
        // sum exp
        eigen_tools::computeFourierAtPhiExp(xd, nu_1, s2ExpIAlphaPhi, d);
#endif
    }
    
    // process: rotate field
    template <int D,
    typename RMatXD_RM =
    Eigen::Matrix<numerical::Real, Eigen::Dynamic, D, Eigen::RowMajor>>
    void rotateField(RMatXD_RM &field, int bufferLine, bool fieldInRTZ,
                     channel::WavefieldCS wcs, bool cartesian) const {
        if (cartesian) {
            if (wcs == channel::WavefieldCS::ENZ) {
                // RTZ (= SPZ) -> ENZ
                rotate(field, bufferLine, 2, mBAzPlus90);
            } else if (wcs == channel::WavefieldCS::xyz) {
                // RTZ (= SPZ) -> XYZ
                rotate(field, bufferLine, 2, -mPhi);
            }
        } else {
            if (fieldInRTZ) {
                if (wcs == channel::WavefieldCS::spz) {
                    // RTZ -> SPZ
                    rotate(field, bufferLine, 1, -mTheta);
                } else if (wcs == channel::WavefieldCS::ENZ) {
                    // RTZ -> ENZ
                    rotate(field, bufferLine, 2, mBAzPlus90);
                } else if (wcs == channel::WavefieldCS::xyz) {
                    // RTZ -> SPZ
                    rotate(field, bufferLine, 1, -mTheta);
                    // SPZ -> XYZ
                    rotate(field, bufferLine, 2, -mPhi);
                }
            } else {
                if (wcs == channel::WavefieldCS::RTZ) {
                    // SPZ -> RTZ
                    rotate(field, bufferLine, 1, mTheta);
                } else if (wcs == channel::WavefieldCS::ENZ) {
                    // SPZ -> RTZ
                    rotate(field, bufferLine, 1, mTheta);
                    // RTZ -> ENZ
                    rotate(field, bufferLine, 2, mBAzPlus90);
                } else if (wcs == channel::WavefieldCS::xyz) {
                    // SPZ -> XYZ
                    rotate(field, bufferLine, 2, -mPhi);
                }
            }
        }
    }
    
    // process: compute channel and feed buffer
    template <int D,
    typename RMatXD_RM =
    Eigen::Matrix<numerical::Real, Eigen::Dynamic, D, Eigen::RowMajor>>
    void computeFeedChannel(const RMatXD_RM &field, int fieldIndex,
                            int bufferLine, int channelIndex, int stationIndex,
                            eigen::RTensor3 &bufferFields) const {
        if (fieldIndex >= 0) {
            sColBuffer.topRows(bufferLine) = field.bblk(fieldIndex);
        } else {
            if (D == 3) {
                if (fieldIndex == -1) {
                    sColBuffer.topRows(bufferLine) =
                    field.topRows(bufferLine).rowwise().norm();
                } else {
                    throw std::runtime_error("Station::computeFeedChannel || "
                                             "Invalid channel setting.");
                }
            } else if (D == 6) {
                if (fieldIndex == -1) {
                    sColBuffer.topRows(bufferLine) =
                    field.bblk(0) + field.bblk(1) + field.bblk(2);
                } else if (fieldIndex == -2) {
                    // J2 = I1 ^ 2 / 3 - I2
                    sColBuffer.topRows(bufferLine) =
                    (field.bblk(0) + field.bblk(1) +
                     field.bblk(2)).array().square() * numerical::Real(1./3);
                    sColBuffer.topRows(bufferLine) -=
                    (field.bblk(0).cwiseProduct(field.bblk(1)) +
                     field.bblk(1).cwiseProduct(field.bblk(2)) +
                     field.bblk(2).cwiseProduct(field.bblk(0)) -
                     field.bblk(3).cwiseProduct(field.bblk(3)) -
                     field.bblk(4).cwiseProduct(field.bblk(4)) -
                     field.bblk(5).cwiseProduct(field.bblk(5)));
                } else {
                    throw std::runtime_error("Station::computeFeedChannel || "
                                             "Invalid channel setting.");
                }
            } else if (D == 9) {
                if (fieldIndex == -1) {
                    sColBuffer.topRows(bufferLine) =
                    field.bblk(0) + field.bblk(4) + field.bblk(8);
                } else {
                    throw std::runtime_error("Station::computeFeedChannel || "
                                             "Invalid channel setting.");
                }
            } else {
                throw std::runtime_error("Station::computeFeedChannel || "
                                         "Invalid channel setting.");
            }
        }
        
        // feed to tensor
        eigen::IArray3 loc = {stationIndex, channelIndex, 0};
        eigen::IArray3 len = {1, 1, bufferLine};
        bufferFields.slice(loc, len).reshape(eigen::IArray1{bufferLine}) =
        Eigen::TensorMap<eigen::RTensor1>(sColBuffer.data(), bufferLine);
    }
    
    
    ///////////////////////// data /////////////////////////
private:
    /////////// for io ///////////
    // key (network.name)
    const std::string mKey;
    
    /////////// for record ///////////
    // inplane interpolation weights
    // NOTE: stations are likely to be located on an element edge such as
    //       the surface, so we only store non-zero inplane weights
    eigen::RRowX mNonZeroWeights = eigen::RRowX(0);
    std::vector<int> mNonZeroIndices;
    
    // azimuth for Fourier interpolation
    const double mPhi;
    
    // 2 * exp(i * alpha * phi) for Fourier interpolation
#ifndef _SAVE_MEMORY
    // precompute
    eigen::CColX m2ExpIAlphaPhi = eigen::CColX(0);
#else
    // on-the-fly
    inline static eigen::CColX s2ExpIAlphaPhi = eigen::CColX(0);
#endif
    
    /////////// for process ///////////
    // angles for rotation to RTZ and ENZ
    const double mTheta;
    const double mBAzPlus90;
    
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    ///////////////////////// rotate /////////////////////////
private:
    // rotate vector
    static void rotate(eigen::RMatX3_RM &V, int nrow, int k, double angle);
    
    // rotate tensor
    static void rotate(eigen::RMatX9_RM &T, int nrow, int k, double angle);
    
    // rotate tensor in Voigt
    static void rotate(eigen::RMatX6_RM &S, int nrow, int k, double angle);
    
    
    ///////////////////////// workspace /////////////////////////
protected:
    // expand workspace for process
    static void expandWorkspaceProcess(int dumpIntv, bool needTensor33) {
        if (sColBuffer.rows() < dumpIntv) {
            sColBuffer.resize(dumpIntv, 1);
        }
        if (sTensor33.rows() < dumpIntv && needTensor33) {
            sTensor33.resize(dumpIntv, 9);
        }
    }
    
private:
    // workspace for process
    inline static eigen::RColX sColBuffer = eigen::RColX(0);
    inline static eigen::RMatX9_RM sTensor33 = eigen::RMatX9_RM(0, 9);
};

#endif /* Station_hpp */
