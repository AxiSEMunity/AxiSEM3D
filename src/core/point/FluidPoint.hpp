//
//  FluidPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  fluid GLL point

#ifndef FluidPoint_hpp
#define FluidPoint_hpp

#include "Point.hpp"
#include "eigen_element.hpp"
class TimeScheme;

class FluidPoint: public Point {
public:
    // constructor
    FluidPoint(int nr, const eigen::DRow2 &crds, int meshTag,
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
        mFields.mPressureSource.setZero();
        mFields.mPressureStore.setZero();
        mFields.mDeltaStore.setZero();
    }
    
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    int sizeComm() const {
        return (int)mFields.mStiff.size();
    }
    
    // feed to mpi buffer
    void feedComm(eigen::CColX &buffer, int &row) const {
        buffer.block(row, 0, mFields.mStiff.rows(), 1) = mFields.mStiff;
        row += mFields.mStiff.rows();
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::CColX &buffer, int &row) {
        mFields.mStiff += buffer.block(row, 0, mFields.mStiff.rows(), 1);
        row += mFields.mStiff.rows();
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
    void scatterDisplToElement(eigen::vec_ar1_CMatPP_RM &displ,
                               int nu_1_element, int ipol, int jpol) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            displ[alpha][0](ipol, jpol) = mFields.mDispl(alpha);
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = mNu_1; alpha < nu_1_element; alpha++) {
            displ[alpha][0](ipol, jpol) = czero;
        }
    }
    
    // gather stiff from element
    void gatherStiffFromElement(const eigen::vec_ar1_CMatPP_RM &stiff,
                                int ipol, int jpol) {
        // add lower orders only
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            mFields.mStiff(alpha) -= stiff[alpha][0](ipol, jpol);
        }
    }
    
    
    /////////////////////////// source ///////////////////////////
    // prepare pressure source
    void preparePressureSource() {
        mFields.mPressureSource = eigen::CColX::Zero(mNu_1, 1);
    }
    
    // add pressure source
    void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(mNu_1, nu_1_pressure);
        mFields.mPressureSource.topRows(nu_1_min) +=
        pressure.block(0, ipnt, nu_1_min, 1);
    }
    
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare pressure output
    void preparePressureOutput() {
        // pressure is mAccel, which may not be allocated by the time scheme
        if (mFields.mPressureStore.rows() == 0) {
            mFields.mPressureStore = eigen::CColX::Zero(mNu_1, 1);
        }
    }
    
    // prepare delta output
    void prepareDeltaOutput() {
        // delta is mStiff, but we need to store it because mStiff is set
        // to zero after dividing by mass
        if (mFields.mDeltaStore.rows() == 0) {
            mFields.mDeltaStore = eigen::CColX::Zero(mNu_1, 1);
        }
    }
    
    // scatter pressure to element
    void scatterPressureToElement(eigen::CMatXN &pressure,
                                  int nu_1_element, int ipnt) const {
        // copy lower orders
        pressure.block(0, ipnt, mNu_1, 1) = mFields.mPressureStore;
        
        // mask higher orders
        pressure.block(mNu_1, ipnt, nu_1_element - mNu_1, 1).setZero();
    }
    
    // scatter delta to element
    void scatterDeltaToElement(eigen::CMatXN &delta,
                               int nu_1_element, int ipnt) const {
        // copy lower orders
        delta.block(0, ipnt, mNu_1, 1) = mFields.mDeltaStore;
        
        // mask higher orders
        delta.block(mNu_1, ipnt, nu_1_element - mNu_1, 1).setZero();
    }
    
    
    /////////////////////////// fields ///////////////////////////
    // fields on a fluid point
    struct Fields {
        eigen::CColX mStiff = eigen::CColX(0, 1);
        eigen::CColX mDispl = eigen::CColX(0, 1);
        eigen::CColX mVeloc = eigen::CColX(0, 1);
        eigen::CColX mAccel = eigen::CColX(0, 1);
        
        // pressure source (to be added to accel)
        eigen::CColX mPressureSource = eigen::CColX(0, 1);
        
        // acceleration storage for pressure output
        // NOTE: duplicated in Newmark
        eigen::CColX mPressureStore = eigen::CColX(0, 1);
        
        // stiffness storage for delta output
        eigen::CColX mDeltaStore = eigen::CColX(0, 1);
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
    // fields on a fluid point
    Fields mFields;
    
    
    /////////////////////////// wavefield scanning ///////////////////////////
public:
    // enable scanning
    void enableScanning() {
        if (!mScanningChi) {
            mScanningChi = std::make_unique<Scanning1D>();
        }
    }
    
    // disable scanning
    void disableScanning() {
        if (mScanningChi) {
            mScanningChi.reset();
            mScanningChi = nullptr;
        }
    }
    
    // do scanning
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) const {
        if (mScanningChi) {
            mScanningChi->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                     maxNumPeaks, mFields.mDispl);
        }
    }
    
    // report scanning Nr
    int reportScanningNr() const {
        if (mScanningChi) {
            return mScanningChi->reportScanningNr(mNr);
        } else {
            return -1;
        }
    }
    
private:
    // scanning
    std::unique_ptr<Scanning1D> mScanningChi = nullptr;
};

#endif /* FluidPoint_hpp */
