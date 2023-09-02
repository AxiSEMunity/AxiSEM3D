//
//  Sponge.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  sponge ABC

#ifndef Sponge_hpp
#define Sponge_hpp

#include "numerical.hpp"
#include "eigen_generic.hpp"
#include "fft.hpp"

//////////////// solid/fluid interface for Sponge ////////////////
template <class PointSF>
struct InterfaceSF {
};

template <>
struct InterfaceSF<SolidPoint> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 3> RMatSF;
    inline static SolverFFTW<numerical::Real, 3> &fftSF = fft::gFFT_3;
};

template <>
struct InterfaceSF<FluidPoint> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RMatSF;
    inline static SolverFFTW<numerical::Real, 1> &fftSF = fft::gFFT_1;
};


///////////////////// Sponge class /////////////////////
template <class PointSF>
class Sponge {
public:
    // 1D constructor
    Sponge(const std::shared_ptr<PointSF> &point, double gamma):
    mPoint(point), mGamma((eigen::RColX(1) << gamma).finished()) {
        // nothing
    }
    
    // 3D constructor
    Sponge(const std::shared_ptr<PointSF> &point, const eigen::DColX &gamma):
    mPoint(point), mGamma(gamma.cast<numerical::Real>()) {
        // check size
        int nr = (int)gamma.rows();
        if (nr != point->getNr()) {
            throw std::runtime_error("Sponge::Sponge || Incompatible sizes.");
        }
        
        // fft
        InterfaceSF<PointSF>::fftSF.addNR(nr);
        
        // workspace
        if (sStiff.rows() < nr) {
            sStiff.resize(nr, sStiff.cols());
            sDispl.resize(nr, sDispl.cols());
            sVeloc.resize(nr, sVeloc.cols());
        }
    }
    
    // get point
    const std::shared_ptr<PointSF> &getPoint() const {
        return mPoint;
    }
    
    // apply ABC
    // must be called after point->computeStiffToAccel(),
    // so here "stiff" has been converted to acceleration
    void apply() const {
        static const numerical::Real two = 2.;
        
        // get fields from point
        auto &stiff = mPoint->getFields().mStiff;
        const auto &veloc = mPoint->getFields().mVeloc;
        const auto &displ = mPoint->getFields().mDispl;
        
        // update acceleration
        if (mGamma.rows() == 1) {
            // 1D
            numerical::Real gamma = mGamma(0);
            stiff -= (two * gamma) * veloc + (gamma * gamma) * displ;
        } else {
            // 3D
            int nr = (int)mGamma.rows();
            InterfaceSF<PointSF>::fftSF.computeC2R(stiff, sStiff, nr);
            InterfaceSF<PointSF>::fftSF.computeC2R(displ, sDispl, nr);
            InterfaceSF<PointSF>::fftSF.computeC2R(veloc, sVeloc, nr);
            sStiff.topRows(nr) -=
            (two * mGamma).asDiagonal() * sVeloc.topRows(nr) +
            mGamma.array().square().matrix().asDiagonal() * sDispl.topRows(nr);
            InterfaceSF<PointSF>::fftSF.computeR2C(sStiff, stiff, nr);
        }
    }
    
private:
    // point
    const std::shared_ptr<PointSF> mPoint;
    
    // gamma
    const eigen::RColX mGamma;
    
    //////////////////// static ////////////////////
    // 3D workspace
    inline static typename
    InterfaceSF<PointSF>::RMatSF sStiff, sVeloc, sDispl;
};

#endif /* Sponge_hpp */
