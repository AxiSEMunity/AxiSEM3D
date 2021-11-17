//
//  Mechanism.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of source mechanism

#include "Mechanism.hpp"
#include "inparam.hpp"

#include "Quad.hpp"
#include "FluidPressure.hpp"
#include "SolidForce.hpp"
#include "SolidMoment.hpp"
#include "Domain.hpp"
#include "Element.hpp"

// build from inparam
std::unique_ptr<const Mechanism> Mechanism::
buildInparam(int sindex, const std::string &sourceName) {
    // short alias
    const InparamYAML &gm = inparam::gInparamSource;
    const std::string &root = ("list_of_sources:[" + bstring::toString(sindex)
                               + "]:" + sourceName + ":mechanism");
    
    // type
    const SM_Type &type =
    gm.getWithLimits<SM_Type>(root + ":type", {
        {"MOMENT_TENSOR", SM_Type::MomentTensor},
        {"FORCE_VECTOR", SM_Type::ForceVector},
        {"FLUID_PRESSURE", SM_Type::FluidPressure}});
    
    // data
    const std::vector<double> &vec = gm.getVector<double>(root + ":data");
    eigen::DColX data = Eigen::Map<const eigen::DColX>(vec.data(), vec.size());
    data *= gm.get<double>(root + ":unit");
    return std::make_unique<Mechanism>(type, data);
}

// release element source
void Mechanism::
release(const eigen::DMat33 &Qzsp, bool sourceOnAxis,
        const eigen::DRowN &inplaneFactor, double phi,
        const Quad &quad, std::unique_ptr<STF> &stf, Domain &domain) const {
    // constants
    using spectral::nPEM;
    using spectral::nPED;
    static const numerical::ComplexD i(0., 1.);
    const numerical::ComplexD i_phi(0., -phi);
    
    // nu + 1
    int nu_1 = quad.getElement()->getNu_1();
    
    // create source
    std::unique_ptr<const ElementSource> src = nullptr;
    if (mType == SM_Type::FluidPressure) {
        // interpolated data
        const eigen::DRowN &p = mData(0) * inplaneFactor;
        // allocate
        if (sourceOnAxis) {
            nu_1 = std::min(nu_1, 1);
        }
        eigen::ZMatXN pattern(nu_1, nPEM);
        // non-axial
        if (sourceOnAxis) {
            pattern.setZero();
        } else {
            for (int alpha = 0; alpha < nu_1; alpha++) {
                pattern.row(alpha) = exp(1. * alpha * i_phi) * p;
            }
        }
        // axial
        if (quad.axial()) {
            // mask
            pattern.block(0, 0, nu_1, nPED).setZero();
            // monopole
            pattern.block(0, 0, 1, nPED) = p.block(0, 0, 1, nPED);
            // axial scaling
            pattern /= (2. * numerical::dPi);
        }
        // element source
        src = std::make_unique<const FluidPressure>
        (stf, quad.getFluidElement(), pattern.cast<numerical::ComplexR>());
    } else if (mType == SM_Type::ForceVector) {
        // rotate data
        const eigen::DColX &fzsp = Qzsp * mData;
        // interpolated data
        const eigen::DRowN &fs = fzsp(1) * inplaneFactor;
        const eigen::DRowN &fp = fzsp(2) * inplaneFactor;
        const eigen::DRowN &fz = fzsp(0) * inplaneFactor;
        // allocate
        if (sourceOnAxis) {
            nu_1 = std::min(nu_1, 2);
        }
        eigen::ZMatXN3 pattern(nu_1, nPEM * 3);
        // non-axial
        if (sourceOnAxis) {
            pattern.setZero();
        } else {
            for (int alpha = 0; alpha < nu_1; alpha++) {
                numerical::ComplexD exp_iAlphaPhi = exp(1. * alpha * i_phi);
                pattern.block(alpha, nPEM * 0, 1, nPEM) = exp_iAlphaPhi * fs;
                pattern.block(alpha, nPEM * 1, 1, nPEM) = exp_iAlphaPhi * fp;
                pattern.block(alpha, nPEM * 2, 1, nPEM) = exp_iAlphaPhi * fz;
            }
        }
        // axial
        if (quad.axial()) {
            // mask
            pattern.block(0, nPEM * 0, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 1, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 2, nu_1, nPED).setZero();
            // monopole
            pattern.block(0, nPEM * 2, 1, nPED) = fz.block(0, 0, 1, nPED);
            // dipole
            if (nu_1 >= 2) {
                pattern.block(1, nPEM * 0, 1, nPED) =
                (fs.block(0, 0, 1, nPED) - i * fp.block(0, 0, 1, nPED)) / 2.;
                pattern.block(1, nPEM * 1, 1, nPED) =
                i * pattern.block(1, nPEM * 0, 1, nPED);
            }
            // axial scaling
            pattern /= (2. * numerical::dPi);
        }
        // element source
        src = std::make_unique<const SolidForce>
        (stf, quad.getSolidElement(), pattern.cast<numerical::ComplexR>());
    } else {
        // rotate data
        static eigen::DMat33 m;
        m(0, 0) = mData(0);
        m(1, 1) = mData(1);
        m(2, 2) = mData(2);
        m(0, 1) = mData(3);
        m(1, 0) = mData(3);
        m(0, 2) = mData(4);
        m(2, 0) = mData(4);
        m(1, 2) = mData(5);
        m(2, 1) = mData(5);
        m = (Qzsp * m * Qzsp.transpose()).eval();
        // interpolated data
        // wanted
        // mss msp msz
        // msp mpp mpz
        // msz mpz mzz
        // stored
        // mzz msz mpz
        // msz mss msp
        // mpz msp mpp
        const eigen::DRowN &mss = m(1, 1) * inplaneFactor;
        const eigen::DRowN &mpp = m(2, 2) * inplaneFactor;
        const eigen::DRowN &mzz = m(0, 0) * inplaneFactor;
        const eigen::DRowN &mpz = m(0, 2) * inplaneFactor;
        const eigen::DRowN &msz = m(0, 1) * inplaneFactor;
        const eigen::DRowN &msp = m(1, 2) * inplaneFactor;
        // allocate
        if (sourceOnAxis) {
            nu_1 = std::min(nu_1, 3);
        }
        eigen::ZMatXN6 pattern(nu_1, nPEM * 6);
        // non-axial
        if (sourceOnAxis) {
            pattern.setZero();
        } else {
            for (int alpha = 0; alpha < nu_1; alpha++) {
                // Voigt notation
                // 0   5   4
                //     1   3
                //         2
                // mss msp msz
                // msp mpp mpz
                // msz mpz mzz
                numerical::ComplexD exp_iAlphaPhi = exp(1. * alpha * i_phi);
                pattern.block(alpha, nPEM * 0, 1, nPEM) = exp_iAlphaPhi * mss;
                pattern.block(alpha, nPEM * 1, 1, nPEM) = exp_iAlphaPhi * mpp;
                pattern.block(alpha, nPEM * 2, 1, nPEM) = exp_iAlphaPhi * mzz;
                pattern.block(alpha, nPEM * 3, 1, nPEM) = exp_iAlphaPhi * mpz;
                pattern.block(alpha, nPEM * 4, 1, nPEM) = exp_iAlphaPhi * msz;
                pattern.block(alpha, nPEM * 5, 1, nPEM) = exp_iAlphaPhi * msp;
            }
        }
        // axial
        if (quad.axial()) {
            // mask
            pattern.block(0, nPEM * 0, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 1, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 2, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 3, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 4, nu_1, nPED).setZero();
            pattern.block(0, nPEM * 5, nu_1, nPED).setZero();
            // monopole
            pattern.block(0, nPEM * 0, 1, nPED) =
            pattern.block(0, nPEM * 1, 1, nPED) =
            (mss.block(0, 0, 1, nPED) + mpp.block(0, 0, 1, nPED)) / 2.;
            pattern.block(0, nPEM * 2, 1, nPED) = mzz.block(0, 0, 1, nPED);
            // dipole
            if (nu_1 >= 2) {
                pattern.block(1, nPEM * 4, 1, nPED) =
                (msz.block(0, 0, 1, nPED) - i * mpz.block(0, 0, 1, nPED)) / 2.;
                pattern.block(1, nPEM * 3, 1, nPED) =
                i * pattern.block(1, nPEM * 4, 1, nPED);
            }
            // quadrupole
            if (nu_1 >= 3) {
                pattern.block(2, nPEM * 0, 1, nPED) =
                (mss.block(0, 0, 1, nPED) - mpp.block(0, 0, 1, nPED)
                 - 2. * i * msp.block(0, 0, 1, nPED)) / 4.;
                pattern.block(2, nPEM * 1, 1, nPED) =
                -pattern.block(2, nPEM * 0, 1, nPED);
                pattern.block(2, nPEM * 5, 1, nPED) =
                i * pattern.block(2, nPEM * 0, 1, nPED);
            }
            // axial scaling
            pattern /= (2. * numerical::dPi);
        }
        // element source
        src = std::make_unique<const SolidMoment>
        (stf, quad.getSolidElement(), pattern.cast<numerical::ComplexR>());
    }
    
    // finally add to domain
    domain.addElementSource(src);
}
