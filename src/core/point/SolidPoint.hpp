//
//  SolidPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid GLL point

#ifndef SolidPoint_hpp
#define SolidPoint_hpp

#include "Point.hpp"
#include "eigen_element.hpp"
class TimeScheme;

class SolidPoint: public Point {
public:
    // constructor
    SolidPoint(int nr, const eigen::DRow2 &crds, int meshTag,
               std::unique_ptr<const Mass> &mass,
               const TimeScheme &timeScheme);
    
public:
    /////////////////////////// measure ///////////////////////////
    // random displ
    void randomDispl() {
        mFields.mDispl.setRandom();
        mFields.mDispl.row(0).imag().setZero();
    }
    
    // random stiff
    void randomStiff() {
        mFields.mStiff.setRandom();
        mFields.mStiff.row(0).imag().setZero();
    }
    
    // reset to zero
    void resetToZero() {
        mFields.mStiff.setZero();
        mFields.mDispl.setZero();
        mFields.mVeloc.setZero();
        mFields.mAccel.setZero();
    }
    
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    int sizeComm() const {
        return (int)mFields.mStiff.size();
    }
    
    // feed to mpi buffer
    void feedComm(eigen::CColX &buffer, int &row) const {
        buffer.block(row, 0, mFields.mStiff.size(), 1) =
        Eigen::Map<const eigen::CColX>(mFields.mStiff.data(),
                                       mFields.mStiff.size());
        row += mFields.mStiff.size();
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::CColX &buffer, int &row) {
        mFields.mStiff +=
        Eigen::Map<const eigen::CMatX3>(&buffer(row), mFields.mStiff.rows(), 3);
        row += mFields.mStiff.size();
    }
    
    
    /////////////////////////// time loop ///////////////////////////
    // check stability
    bool stable() const {
        return mFields.mDispl.allFinite();
    }
    
    // stiff to accel
    void computeStiffToAccel();
    
    
    /////////////////////////// element ///////////////////////////
    // scatter displ to element
    void scatterDisplToElement(eigen::vec_ar3_CMatPP_RM &displ,
                               int nu_1_element, int ipol, int jpol) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            displ[alpha][0](ipol, jpol) = mFields.mDispl(alpha, 0);
            displ[alpha][1](ipol, jpol) = mFields.mDispl(alpha, 1);
            displ[alpha][2](ipol, jpol) = mFields.mDispl(alpha, 2);
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = mNu_1; alpha < nu_1_element; alpha++) {
            displ[alpha][0](ipol, jpol) = czero;
            displ[alpha][1](ipol, jpol) = czero;
            displ[alpha][2](ipol, jpol) = czero;
        }
    }
    
    // gather stiff from element
    void gatherStiffFromElement(const eigen::vec_ar3_CMatPP_RM &stiff,
                                int ipol, int jpol) {
        // add lower orders only
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            mFields.mStiff(alpha, 0) -= stiff[alpha][0](ipol, jpol);
            mFields.mStiff(alpha, 1) -= stiff[alpha][1](ipol, jpol);
            mFields.mStiff(alpha, 2) -= stiff[alpha][2](ipol, jpol);
        }
    }
    
    
    /////////////////////////// source ///////////////////////////
    // add force source (external)
    void addForceSource(const eigen::CMatXN3 &force,
                        int nu_1_force, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(mNu_1, nu_1_force);
        mFields.mStiff.topRows(nu_1_min) +=
        force(Eigen::seqN(Eigen::fix<0>, nu_1_min),
              Eigen::seqN(ipnt, Eigen::fix<3>, Eigen::fix<spectral::nPEM>));
    }
    
    
    /////////////////////////// fields ///////////////////////////
    // fields on a solid point
    struct Fields {
        eigen::CMatX3 mStiff = eigen::CMatX3(0, 3);
        eigen::CMatX3 mDispl = eigen::CMatX3(0, 3);
        eigen::CMatX3 mVeloc = eigen::CMatX3(0, 3);
        eigen::CMatX3 mAccel = eigen::CMatX3(0, 3);
    };
    
    // get
    const Fields &getFields() const {
        return mFields;
    }
    
    // set
    Fields &getFields() {
        return mFields;
    }
    
private:
    // fields on a solid point
    Fields mFields;
    
    
    /////////////////////////// wavefield scanning ///////////////////////////
public:
    // enable scanning
    void enableScanning() {
        if (!mScanningS) {
            mScanningS = std::make_unique<Scanning1D>();
            mScanningP = std::make_unique<Scanning1D>();
            mScanningZ = std::make_unique<Scanning1D>();
        }
    }
    
    // do scanning
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) const {
        if (mScanningS) {
            mScanningS->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(0));
            mScanningP->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(1));
            mScanningZ->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(2));
        }
    }
    
    // report scanning Nr
    int reportScanningNr() const {
        if (mScanningS) {
            int nrS = mScanningS->reportScanningNr(mNr);
            int nrP = mScanningP->reportScanningNr(mNr);
            int nrZ = mScanningZ->reportScanningNr(mNr);
            return std::max(std::max(nrS, nrP), nrZ);
        } else {
            return -1;
        }
    }
    
private:
    // scanning
    std::unique_ptr<Scanning1D> mScanningS = nullptr;
    std::unique_ptr<Scanning1D> mScanningP = nullptr;
    std::unique_ptr<Scanning1D> mScanningZ = nullptr;
};

#endif /* SolidPoint_hpp */
