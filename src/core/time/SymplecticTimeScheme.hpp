//
//  SymplecticTimeScheme.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  symplectic time scheme

#ifndef SymplecticTimeScheme_hpp
#define SymplecticTimeScheme_hpp

#include "TimeScheme.hpp"
#include "numerical.hpp"
#include <vector>

class SymplecticTimeScheme: public TimeScheme {
public:
    // constructor
    SymplecticTimeScheme(int verboseInterval, int stabilityInterval):
    TimeScheme(verboseInterval, stabilityInterval) {
        // nothing
    }
    
    // solve
    void solve() const;
    
    
    //////////////////// point ////////////////////
    // create fields on a point
    template <class Point>
    static void createFields(Point &point) {
        auto &f = point.getFields();
        int ndim = (int)f.mStiff.cols();
        int nu_1 = point.getNu_1();
        f.mStiff.resize(nu_1, ndim);
        f.mDispl.resize(nu_1, ndim);
        f.mVeloc.resize(nu_1, ndim);
        f.mStiff.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
    }
    
    // update fields on a point
    template <class Point>
    static void update(Point &point,
                       numerical::Real pi_dt, numerical::Real kappa_dt) {
        auto &f = point.getFields();
        // update dt
        f.mVeloc += pi_dt * f.mStiff;
        f.mDispl += kappa_dt * f.mVeloc;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
    
private:
    // launch fields on points
    template <class Point>
    void launch(const std::vector<std::shared_ptr<Point>> &points) const {
        static const numerical::Real kappa_dt = sKappa[0] * mDt;
        for (const std::shared_ptr<Point> &point: points) {
            auto &f = point->getFields();
            f.mDispl += kappa_dt * f.mVeloc;
        }
    }
    
    // update fields on points
    template <class Point>
    void update(const std::vector<std::shared_ptr<Point>> &points,
                int iSubIter) const {
        const numerical::Real pi_dt = sPi[iSubIter] * mDt;
        const numerical::Real kappa_dt = sKappa[iSubIter] * mDt;
        for (const std::shared_ptr<Point> &point: points) {
            update(*point, pi_dt, kappa_dt);
        }
    }
    
    
    //////////////////////////////////////////
    ////////////////// static ////////////////
    //////////////////////////////////////////
    
private:
    inline static const double sAlpha = 0.1786178958448091;
    inline static const double sBeta = -0.2123418310626054;
    inline static const double sGamma = -0.06626458266981849;
    
    inline static const std::vector<numerical::Real> sKappa = {
        (numerical::Real)sAlpha,
        (numerical::Real)sGamma,
        (numerical::Real)(1. - 2. * (sGamma + sAlpha)),
        (numerical::Real)sGamma,
        (numerical::Real)sAlpha};
    
    inline static const std::vector<numerical::Real> sPi = {
        (numerical::Real)0., // unused front padding
        (numerical::Real)(0.5 - sBeta),
        (numerical::Real)sBeta,
        (numerical::Real)sBeta,
        (numerical::Real)(0.5 - sBeta)};
};

#endif /* SymplecticTimeScheme_hpp */
