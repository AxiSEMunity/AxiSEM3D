//
//  GaussianSTF.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  Gaussian source-time function

#ifndef GaussianSTF_hpp
#define GaussianSTF_hpp

#include "STF.hpp"
#include "numerical.hpp"
#include "bstring.hpp"
#include <cmath>

class GaussianSTF: public STF {
public:
    // constructor
    GaussianSTF(double halfDuration, double decayFactor,
                double timeShift, int derivativeOrder):
    mHalfDuration(halfDuration), mDecayFactor(decayFactor),
    mTimeShift(timeShift), mDerivativeOrder(derivativeOrder) {
        // nothing
    }
    
    // get start time
    double getStartTime() const {
        return -mHalfDuration * 2.5 + mTimeShift;
    }
    
    // get value
    numerical::Real getValue(double time) {
        static const double invSqrtPi = 1. / sqrt(numerical::dPi);
        const double f = mDecayFactor / mHalfDuration;
        double t = time - mTimeShift;
        double ft = f * t;
        double value = 0.;
        if (mDerivativeOrder == -1) {
            // erf
            value = erf(ft) * 0.5 + 0.5;
        } else if (mDerivativeOrder == 0) {
            // Gaussian
            value = invSqrtPi * f * exp(-ft * ft);
        } else if (mDerivativeOrder == 1) {
            // 1st derivative of Gaussian
            value = invSqrtPi * pow(f, 3) * exp(-ft * ft) * (-2. * t);
        } else if (mDerivativeOrder == 2) {
            // 2nd derivative of Gaussian (Richer wavelet)
            value = invSqrtPi * pow(f, 3) * exp(-ft * ft) * (4. * ft * ft - 2.);
        } else {
            throw std::runtime_error("GaussianSTF::getValue || "
                                     "Unsupported derivative order.");
        }
        return (numerical::Real)value;
    }
    
    // verbose
    std::string verbose() const {
        using namespace bstring;
        std::stringstream ss;
        ss << boxSubTitle(2, "Source-time function");
        ss << boxEquals(4, 18, "class name", "GaussianSTF");
        ss << boxEquals(4, 18, "half duration", mHalfDuration);
        ss << boxEquals(4, 18, "decay factor", mDecayFactor);
        ss << boxEquals(4, 18, "time shift", mTimeShift);
        // order
        std::stringstream ssd;
        ssd << mDerivativeOrder;
        // expression
        std::stringstream sse;
        static const double invSqrtPi = 1. / sqrt(numerical::dPi);
        const double f = mDecayFactor / mHalfDuration;
        std::string t;
        if (std::abs(mTimeShift) < numerical::dEpsilon) {
            t = "t";
        } else if (mTimeShift > 0.) {
            t = "(t - " + toString(mTimeShift) + ")";
        } else {
            t = "(t + " + toString(-mTimeShift) + ")";
        }
        if (mDerivativeOrder == -1) {
            ssd << " (Error function)";
            sse << "(Erf[" << f << "*" << t + "] + 1) / 2";
        } else if (mDerivativeOrder == 0) {
            ssd << " (Gaussian)";
            sse << invSqrtPi * f;
            sse << " * Exp[-" << f * f << " * " << t << " ^ 2]";
        } else if (mDerivativeOrder == 1) {
            ssd << " (1st derivative of Gaussian)";
            sse << -2. * invSqrtPi * pow(f, 3) << " * " << t;
            sse << " * Exp[-" << f * f << " * " << t << " ^ 2]";
        } else if (mDerivativeOrder == 2) {
            ssd << " (2nd derivative of Gaussian; Ricker)";
            sse << 2. * invSqrtPi * pow(f, 3);
            sse << " * (" << 2. * f * f << " * " << t << " ^ 2 - 1)";
            sse << " * Exp[-" << f * f << " * " << t << " ^ 2]";
        } else {
            throw std::runtime_error("GaussianSTF::verbose || "
                                     "Unsupported derivative order.");
        }
        ss << boxEquals(4, 18, "derivative order", ssd.str());
        ss << boxEquals(4, 18, "expression of f(t)", sse.str());
        return ss.str();
    }
    
private:
    // parameters
    const double mHalfDuration;
    const double mDecayFactor;
    const double mTimeShift;
    const int mDerivativeOrder;
};

#endif /* GaussianSTF_hpp */
