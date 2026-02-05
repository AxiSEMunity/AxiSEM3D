//
//  geodesy.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/16/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  geodesy of the model and transformations between these CSs:
//  * geographic: (lat, lon), north pole on top
//  * geocentric: (theta, phi), north pole on top
//  * source-centered: (distance, azimuth), source on top

#ifndef geodesy_hpp
#define geodesy_hpp

#include "eigen_generic.hpp"
#include "vector_tools.hpp"

class ExodusMesh;

namespace eigen {
    typedef Eigen::Matrix<double, 1, 3> DRow3;
    typedef Eigen::Matrix<double, 3, 3> DMat33;
}

namespace geodesy {
    // internal data
    namespace internal {
        // Cartesain
        extern bool iCartesian;
        
        // radius
        extern double iOuterRadius;
        extern double iOuterSolidRadius;
        extern double iInnerRadius;
        
        // ellipticity
        extern double iOuterFlattening;
        extern std::vector<double> iEllipR;
        extern std::vector<double> iEllipF;
        
        // geographic
        // (lat, lon, r) on positive z-axis
        extern eigen::DRow3 iLatLonRadiusAxisZ;
        // Q matrix from source-centered to geographic
        extern eigen::DMat33 iSrc2GeoQ;
    }
    
    // get
    inline bool isCartesian() {
        return internal::iCartesian;
    }
    
    // get
    inline double getOuterRadius() {
        return internal::iOuterRadius;
    }
    
    // get
    inline double getOuterSolidRadius() {
        return internal::iOuterSolidRadius;
    }
    
    // get
    inline double getInnerRadius() {
        return internal::iInnerRadius;
    }
    
    // get
    inline double getOuterFlattening() {
        return internal::iOuterFlattening;
    }
    
    // compute flattenings at given radii
    template <class Vec>
    Vec computeFlattening(const Vec &r) {
        Vec f = r;
        int index0 = -1, index1 = -1;
        double factor0 = 0., factor1 = 0.;
        for (int ir = 0; ir < r.size(); ir++) {
            try {
                vector_tools::linearInterpSorted(internal::iEllipR, r[ir],
                                                 index0, index1,
                                                 factor0, factor1);
                f[ir] = (internal::iEllipF[index0] * factor0 +
                         internal::iEllipF[index1] * factor1);
            } catch (...) {
                // this may occur if using undulated r
                if (r[ir] > internal::iEllipR.back()) {
                    f[ir] = internal::iEllipF.back();
                } else {
                    f[ir] = internal::iEllipF.front();
                }
            }
        }
        return f;
    }
    
    // setup
    void setup(ExodusMesh &exMesh);
    
    // verbose
    std::string verbose(const eigen::DColX &discontinuities);
    
    
    /////////////////////////////// coords ///////////////////////////////
    // (s, z) -> (r, theta)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat sz2rtheta(const Mat &sz, bool sInFirstColumn,
                  int ids = 0, int idz = 1, int idr = 0, int idt = 1) {
        // numeric
        T epsilon = numerical::epsilon<T>();
        
        // allocate
        Mat rt(sz.rows(), sz.cols());
        if (sInFirstColumn) {
            // s in the first column
            rt.col(idr) = (sz.col(ids).array().square() +
                           sz.col(idz).array().square()).sqrt();
            rt.col(idt) = (rt.col(idr).array() < epsilon).select
            ((T)0., (sz.col(idz).cwiseQuotient
                     (rt.col(idr).cwiseMax(epsilon)).array().acos()));
        } else {
            // s in the first row
            rt.row(idr) = (sz.row(ids).array().square() +
                           sz.row(idz).array().square()).sqrt();
            rt.row(idt) = (rt.row(idr).array() < epsilon).select
            ((T)0., (sz.row(idz).cwiseQuotient
                     (rt.row(idr).cwiseMax(epsilon)).array().acos()));
        }
        return rt;
    }
    
    // (x, y) -> (s, phi)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat xy2sphi(const Mat &xy, bool xInFirstColumn, T undefined,
                int idx = 0, int idy = 1, int ids = 0, int idp = 1) {
        // numeric
        T epsilon = numerical::epsilon<T>();
        
        // lambda
        auto lambda_atan2 = [epsilon, undefined](T x, T y) {
            if (std::abs(x) < epsilon && std::abs(y) < epsilon) {
                // undefined result
                return undefined;
            }
            return atan2(y, x);
        };
        
        // allocate
        Mat sp(xy.rows(), xy.cols());
        if (xInFirstColumn) {
            // x in the first column
            sp.col(ids) = (xy.col(idx).array().square() +
                           xy.col(idy).array().square()).sqrt();
            sp.col(idp) = xy.col(idx).binaryExpr(xy.col(idy), lambda_atan2);
        } else {
            // x in the first row
            sp.row(ids) = (xy.row(idx).array().square() +
                           xy.row(idy).array().square()).sqrt();
            sp.row(idp) = xy.row(idx).binaryExpr(xy.row(idy), lambda_atan2);
        }
        return sp;
    }
    
    // (t, p, r) -> (lat, lon, r)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat tpr2llr(const Mat &tpr, bool ellipticity) {
        static const T halfPi = (T)(numerical::dPi / 2.);
        static const T degree = (T)(numerical::dDegree);
        static const T one = (T)(1.);
        
        // allocate
        Mat llr = Mat(tpr.rows(), tpr.cols());
        // latitude
        llr.col(0) = halfPi - tpr.col(0).array();
        if (ellipticity) {
            // ellipticity correction
            typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TColX;
            const TColX &f = computeFlattening(tpr.col(2).eval());
            const TColX &f1sq = (f.array() - one).square();
            llr.col(0) = (llr.col(0).array().tan() / f1sq.array()).atan();
        }
        llr.col(0) /= degree;
        // longitude
        llr.col(1) = tpr.col(1) / degree;
        // radius
        llr.col(2) = tpr.col(2);
        return llr;
    }
    
    // (lat, lon, r) -> (t, p, r)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat llr2tpr(const Mat &llr, bool ellipticity) {
        static const T halfPi = (T)(numerical::dPi / 2.);
        static const T degree = (T)(numerical::dDegree);
        static const T one = (T)(1.);
        
        // allocate
        Mat tpr = Mat(llr.rows(), llr.cols());
        // latitude
        tpr.col(0) = llr.col(0) * degree;
        if (ellipticity) {
            // ellipticity correction
            typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TColX;
            const TColX &f = computeFlattening(llr.col(2).eval());
            const TColX &f1sq = (f.array() - one).square();
            tpr.col(0) = (tpr.col(0).array().tan() * f1sq.array()).atan();
        }
        tpr.col(0) = halfPi - tpr.col(0).array();
        // longitude
        tpr.col(1) = llr.col(1) * degree;
        // radius
        tpr.col(2) = llr.col(2);
        return tpr;
    }
    
    // (s, phi, z) -> (lat, lon, r)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat spz2llr(const Mat &spzSrc, bool ellipticity, bool lon360) {
        // (s, phi, z) -> (x, y, z)
        Mat xyzSrc(spzSrc.rows(), spzSrc.cols());
        xyzSrc.col(0) = spzSrc.col(0).array() * spzSrc.col(1).array().cos();
        xyzSrc.col(1) = spzSrc.col(0).array() * spzSrc.col(1).array().sin();
        xyzSrc.col(2) = spzSrc.col(2);
        
        // (x, y, z) source-centered -> (x, y, z) geographic
        const Mat &xyzGeo = xyzSrc * internal::iSrc2GeoQ;
        
        // (x, y, z) -> (s, p, z)
        T err = std::numeric_limits<T>::lowest();
        Mat spzGeo = xy2sphi(xyzGeo, true, err);
        // use source-centered for undefined
        spzGeo.col(1) = ((spzGeo.col(1).array() < err * .9).
                         select(spzSrc.col(1), spzGeo.col(1)));
        spzGeo.col(2) = xyzGeo.col(2);
        
        // (s, p, z) -> (t, p, r)
        Mat tprGeo = sz2rtheta(spzGeo, true, 0, 2, 2, 0);
        tprGeo.col(1) = spzGeo.col(1);
        
        // (t, p, r) -> (lat, lon, r)
        Mat llr = tpr2llr(tprGeo, ellipticity);
        // re-compute radius with initial input to improve accuracy
        llr.col(2) = (spzSrc.col(0).array().square() +
                      spzSrc.col(2).array().square()).sqrt();
        // longitude range
        if (lon360) {
            // 0 ~ 360
            llr.col(1) = (llr.col(1).array() < 0.).
            select(llr.col(1).array() + (T)360., llr.col(1).array());
        } else {
            // -180 ~ 180
            llr.col(1) = (llr.col(1).array() > 180.).
            select(llr.col(1).array() - (T)360., llr.col(1).array());
        }
        return llr;
    }
    
    // (lat, lon, r) -> (s, phi, z)
    template <typename Mat, typename T = typename Mat::Scalar>
    Mat llr2spz(const Mat &llr, bool ellipticity) {
        // (lat, lon, r) -> (t, p, r)
        Mat tprGeo = llr2tpr(llr, ellipticity);
        
        // (t, p, r) -> (x, y, z)
        Mat xyzGeo(tprGeo.rows(), tprGeo.cols());
        xyzGeo.col(0) = (tprGeo.col(2).array() * tprGeo.col(0).array().sin() *
                         tprGeo.col(1).array().cos());
        xyzGeo.col(1) = (tprGeo.col(2).array() * tprGeo.col(0).array().sin() *
                         tprGeo.col(1).array().sin());
        xyzGeo.col(2) = (tprGeo.col(2).array() * tprGeo.col(0).array().cos());
        
        // (x, y, z) source-centered -> (x, y, z) geographic
        const Mat &xyzSrc = xyzGeo * internal::iSrc2GeoQ.transpose();
        
        // (x, y, z) -> (s, p, z)
        T err = std::numeric_limits<T>::lowest();
        Mat spzSrc = xy2sphi(xyzSrc, true, err);
        // use geographic for undefined
        spzSrc.col(1) = ((spzSrc.col(1).array() < err * .9).
                         select(tprGeo.col(1), spzSrc.col(1)));
        spzSrc.col(2) = xyzSrc.col(2);
        return spzSrc;
    }
    
    // back azimuth
    template <typename Mat, typename T = typename Mat::Scalar,
    typename Col = Eigen::Matrix<T, Eigen::Dynamic, 1>>
    Col backAzimuth(const Mat &llr, bool ellipticity) {
        // event
        const eigen::DRow3 &srctpr =
        llr2tpr(internal::iLatLonRadiusAxisZ, ellipticity);
        double d = sin(srctpr(1));
        double e = -cos(srctpr(1));
        double f = -sin(srctpr(0));
        double c = cos(srctpr(0));
        double a = f * e;
        double b = -f * d;
        // station
        const Mat &tpr = llr2tpr(llr, ellipticity);
        Col dd = tpr.array().col(1).sin();
        Col ee = -tpr.array().col(1).cos();
        Col ff = -tpr.array().col(0).sin();
        Col cc = tpr.array().col(0).cos();
        Col gg = -cc.array() * ee.array();
        Col hh = cc.array() * dd.array();
        // baz
        Mat xy = Mat(llr.rows(), llr.cols());
        xy.col(0) = ((a - gg.array()).square() +
                     (b - hh.array()).square() +
                     (c - ff.array()).square() - (T)2.);
        xy.col(1) = ((a - dd.array()).square() +
                     (b - ee.array()).square() +
                     c * c - (T)2.);
        Mat sp = xy2sphi(xy, true, (T)0.);
        return sp.col(1);
    }
}

#endif /* geodesy_hpp */
