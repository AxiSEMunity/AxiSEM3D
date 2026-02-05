//
//  PhysicalProperty.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 4/4/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  physical property on element
//  NOTE: 1) data are stored pointwise by std::array<eigen::DColX, nPEM>
//        2) data can only be empty (rows=0), 1D (rows=1) or 3D (rows=nr)
//        3) computed elemental data are stored in eigen::DMatXN,
//           either 1D (row=1) or 3D (row=nr_max)
//        4) getPointwise() and getElemental() return 1D zero if data is empty

#ifndef PhysicalProperty_hpp
#define PhysicalProperty_hpp

#include "spectrals.hpp"

// property without nodal
template <int M>
class PhysicalProperty {
protected:
    typedef std::array<eigen::DColX, M> arM_DColX;
    typedef Eigen::Matrix<double, Eigen::Dynamic, M> DMatXM;
    
public:
    // constructor
    PhysicalProperty(): mGLL() {
        mGLL.fill(eigen::DColX(0));
    }
    
    // constructor
    PhysicalProperty(const arM_DColX &gll): mGLL(gll) {
        // nothing
    }
    
    // copy constructor
    PhysicalProperty(const PhysicalProperty &other): mGLL(other.mGLL) {
        // nothing
    }
    
    // bool
    operator bool() const {
        return mGLL[0].rows() != 0;
    }
    
    // get pointwise
    virtual arM_DColX getPointwise() const {
        if (*this) {
            // GLL
            return mGLL;
        } else {
            // return zero array with 1D size
            arM_DColX zero;
            zero.fill(eigen::DColX::Zero(1));
            return zero;
        }
    }
    
    // get elemental
    virtual DMatXM getElemental() const {
        if (*this) {
            // max nr on element
            int maxNr = 0;
            for (int ipnt = 0; ipnt < M; ipnt++) {
                maxNr = std::max(maxNr, (int)mGLL[ipnt].rows());
            }
            // allocate
            DMatXM mat(maxNr, M);
            for (int ipnt = 0; ipnt < M; ipnt++) {
                // use linear interpolation to expand dimension
                // because using fft may cause over-shooting
                int nr = (int)mGLL[ipnt].rows();
                linInterpPhi(mGLL[ipnt], mat, 0, ipnt, nr, maxNr);
            }
            return mat;
        } else {
            // return zero matrix with 1D size
            return DMatXM::Zero(1, M);
        }
    }
    
    // set GLL
    void setGLL(const arM_DColX &gll) {
        mGLL = gll;
        tryReduceTo1D();
    }
    
    // add GLL
    void addGLL(const arM_DColX &gll) {
        for (int ipnt = 0; ipnt < M; ipnt++) {
            op1D_3D::addTo(gll[ipnt], mGLL[ipnt]);
        }
        tryReduceTo1D();
    }
    
    // linear interpolation over phi
    template <class MatIn, class MatOut>
    static void linInterpPhi(const MatIn &dataIn, MatOut &dataOut,
                             int colIn, int colOut, int nrIn, int nrOut) {
        if (nrIn == nrOut) {
            // same size, direct copy
            dataOut.col(colOut) = dataIn.col(colIn);
        } else {
            // linear interpolation
            double dPhiIn = 1. / nrIn;
            double dPhiOut = 1. / nrOut;
            for (int alpha = 0; alpha < nrOut; alpha++) {
                double phi = alpha * dPhiOut;
                int loc0 = (int)(phi / dPhiIn);
                int loc1 = (loc0 == nrIn - 1) ? 0 : loc0 + 1;
                double val0 = dataIn(loc0, colIn);
                double val1 = dataIn(loc1, colIn);
                dataOut(alpha, colOut) = (val0 + (val1 - val0) / dPhiIn *
                                          (phi - loc0 * dPhiIn));
            }
        }
    }
    
private:
    // reduce 3D to 1D if possible
    void tryReduceTo1D() {
        // check 1D or 3D
        bool all1D = true;
        for (int ipnt = 0; ipnt < M; ipnt++) {
            if ((mGLL[ipnt].array() - mGLL[ipnt](0)).matrix().norm() >
                mGLL[ipnt].norm() * numerical::dEpsilon) {
                all1D = false;
                break;
            }
        }
        // reduce to 1D if all are 1D
        if (all1D) {
            for (int ipnt = 0; ipnt < M; ipnt++) {
                mGLL[ipnt].conservativeResize(1);
            }
        }
    }
    
protected:
    // GLL
    arM_DColX mGLL;
};


// property with nodal, M = spectral::nPEM
class NodalPhysicalProperty: public PhysicalProperty<spectral::nPEM> {
public:
    // constructor with nodal
    NodalPhysicalProperty(const eigen::DRow4 &nodal, bool axial):
    PhysicalProperty<spectral::nPEM>(),
    mNodal(nodal), mAxial(axial) {
        // nothing
    }
    
    // copy constructor
    NodalPhysicalProperty(const NodalPhysicalProperty &other):
    PhysicalProperty<spectral::nPEM>(other),
    mNodal(other.mNodal), mAxial(other.mAxial) {
        // nothing
    }
    
    // get pointwise from nodal
    arM_DColX getPointwiseNodal() const {
        const eigen::DRowN &mat = spectrals::interpolateGLL(mNodal, mAxial);
        eigen::arN_DColX ar;
        for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
            ar[ipnt] = (eigen::DColX(1) << mat(ipnt)).finished();
        }
        return ar;
    }
    
    // get elemental from nodal
    eigen::DMatXN getElementalNodal() const {
        return spectrals::interpolateGLL(mNodal, mAxial);
    }
    
    // get pointwise
    arM_DColX getPointwise() const {
        if (*this) {
            // preferentially use GLL
            return PhysicalProperty<spectral::nPEM>::getPointwise();
        } else {
            // use nodal without GLL
            return getPointwiseNodal();
        }
    }
    
    // get elemental
    eigen::DMatXN getElemental() const {
        if (*this) {
            // preferentially use GLL
            return PhysicalProperty<spectral::nPEM>::getElemental();
        } else {
            // use nodal without GLL
            return getElementalNodal();
        }
    }
    
    // get nodal
    const eigen::DRow4 getNodalData() const {
        return mNodal;
    }
    
    // get axial
    bool axial() const {
        return mAxial;
    }
    
    
    /////////////////////////// operators ///////////////////////////
    // operator * scalar
    NodalPhysicalProperty operator *(double factor) const {
        // copy this
        NodalPhysicalProperty result(*this);
        // nodal
        result.mNodal *= factor;
        // GLL
        if (result) {
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                result.mGLL[ipnt] *= factor;
            }
        }
        return result;
    }
    
    // operator / scalar
    NodalPhysicalProperty operator /(double factor) const {
        return *this * (1. / factor);
    }
    
    // pow
    NodalPhysicalProperty pow(double p) const {
        // copy this
        NodalPhysicalProperty result(*this);
        // nodal
        result.mNodal = result.mNodal.array().pow(p);
        // GLL
        if (result) {
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                result.mGLL[ipnt] = result.mGLL[ipnt].array().pow(p);
            }
        }
        return result;
    }
    
    // operator +
    NodalPhysicalProperty operator +(const NodalPhysicalProperty &other) const {
        // copy this
        NodalPhysicalProperty result(*this);
        // nodal
        result.mNodal += other.mNodal;
        // GLL: no need to compute GLL if both are nodal only
        if (*this || other) {
            const eigen::arN_DColX &gllThis = this->getPointwise();
            const eigen::arN_DColX &gllOther = other.getPointwise();
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                result.mGLL[ipnt] = gllThis[ipnt];
                op1D_3D::addTo(gllOther[ipnt], result.mGLL[ipnt]);
            }
        }
        return result;
    }
    
    // operator -
    NodalPhysicalProperty operator -(const NodalPhysicalProperty &other) const {
        return *this + other * (-1.);
    }
    
    // operator *
    NodalPhysicalProperty operator *(const NodalPhysicalProperty &other) const {
        // copy this
        NodalPhysicalProperty result(*this);
        // nodal
        result.mNodal.array() *= other.mNodal.array();
        // GLL: no need to compute GLL if both are nodal only
        if (*this || other) {
            const eigen::arN_DColX &gllThis = this->getPointwise();
            const eigen::arN_DColX &gllOther = other.getPointwise();
            for (int ipnt = 0; ipnt < spectral::nPEM; ipnt++) {
                op1D_3D::times(gllThis[ipnt], gllOther[ipnt],
                               result.mGLL[ipnt]);
            }
        }
        return result;
    }
    
    // operator /
    NodalPhysicalProperty operator /(const NodalPhysicalProperty &other) const {
        return *this * other.pow(-1.);
    }
    
private:
    // axial flag
    bool mAxial;
    
    // nodal
    eigen::DRow4 mNodal;
};

#endif /* PhysicalProperty_hpp */
