//
//  NewmarkTimeScheme.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Newmark time scheme

#ifndef NewmarkTimeScheme_hpp
#define NewmarkTimeScheme_hpp

#include "TimeScheme.hpp"
#include "numerical.hpp"

class NewmarkTimeScheme: public TimeScheme {
public:
    // constructor
    NewmarkTimeScheme(int verboseInterval, int stabilityInterval):
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
        f.mAccel.resize(nu_1, ndim);
        f.mStiff.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
        f.mAccel.setZero();
    }
    
    // update fields on a point
    template <class Point>
    static void update(Point &point,
                       numerical::Real dt, numerical::Real half_dt,
                       numerical::Real half_dt_dt) {
        auto &f = point.getFields();
        // update dt
        f.mVeloc += half_dt * (f.mAccel + f.mStiff);
        f.mAccel = f.mStiff;
        f.mDispl += dt * f.mVeloc + half_dt_dt * f.mAccel;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
    
private:
    // update fields on points
    template <class Point>
    void update(const std::vector<std::shared_ptr<Point>> &points) const {
        static const numerical::Real dt = mDt;
        static const numerical::Real half_dt = .5 * mDt;
        static const numerical::Real half_dt_dt = .5 * mDt * mDt;
        for (const std::shared_ptr<Point> &point: points) {
            update(*point, dt, half_dt, half_dt_dt);
        }
    }
};

#endif /* NewmarkTimeScheme_hpp */
