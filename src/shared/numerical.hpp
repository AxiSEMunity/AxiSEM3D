//
//  numerical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/10/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  numerical precision and constants

#ifndef numerical_hpp
#define numerical_hpp

#include <complex>

namespace numerical {
    // large int
    typedef std::ptrdiff_t Int;
    
    // solver precision
#ifdef _USE_DOUBLE
    typedef double Real;
#else
    typedef float Real;
#endif
    
    // epsilon
    const float fEpsilon = 1e-5f;
    const double dEpsilon = 1e-12;
    
    // get epsilon
    template <typename T>
    inline T epsilon() = delete;
    
    // float specialization
    template <>
    inline float epsilon<float>() {
        return fEpsilon;
    }
    
    // double specialization
    template <>
    inline double epsilon<double>() {
        return dEpsilon;
    }
    
    // complex
    typedef std::complex<Real> ComplexR;
    typedef std::complex<double> ComplexD;
    
    // constants
    const double dPi = atan(1.) * 4.;
    const double dDegree = dPi / 180.;
    
    // error or undefined number
    const double dErr = -12345e10;
}

#endif /* numerical_hpp */
