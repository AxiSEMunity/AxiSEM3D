//
//  point_time.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/11/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  time-scheme-dependent point methods

#ifndef point_time_hpp
#define point_time_hpp

#include "NewmarkTimeScheme.hpp"
#include "SymplecticTimeScheme.hpp"
#include "timer.hpp"
#include "bstring.hpp"

namespace point_time {
    // create fields
    template <class Point>
    void createFields(Point &point, const TimeScheme &timeScheme) {
        const std::string &className = bstring::typeName(timeScheme);
        if (className == "NewmarkTimeScheme") {
            NewmarkTimeScheme::createFields<Point>(point);
        } else if (className == "SymplecticTimeScheme") {
            SymplecticTimeScheme::createFields<Point>(point);
        } else {
            throw std::runtime_error("point_time::createFields || "
                                     "Unknown derived class of TimeScheme: "
                                     + className);
        }
    }
    
    // measure cost of a point
    template <class Point>
    double measure(Point &point, int count, const TimeScheme &timeScheme) {
        const numerical::Real half = (numerical::Real).5;
        // set stiffness to random
        point.randomStiff();
        // class judgement must be outside measurement
        SimpleTimer tm;
        const std::string &className = bstring::typeName(timeScheme);
        if (className == "NewmarkTimeScheme") {
            tm.start();
            for (int irep = 0; irep < count; irep++) {
                point.computeStiffToAccel();
                NewmarkTimeScheme::update(point, half, half, half);
            }
            tm.pause();
        } else if (className == "SymplecticTimeScheme") {
            tm.start();
            for (int irep = 0; irep < count; irep++) {
                point.computeStiffToAccel();
                SymplecticTimeScheme::update(point, half, half);
            }
            tm.pause();
        } else {
            throw std::runtime_error("point_time::measure || "
                                     "Unknown derived class of TimeScheme: "
                                     + className);
        }
        // reset fields to zero
        point.resetToZero();
        // return total, not divided by count
        return tm.elapsedTotal();
    }
}

#endif /* point_time_hpp */
