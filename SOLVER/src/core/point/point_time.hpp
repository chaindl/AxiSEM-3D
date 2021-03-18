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
#include "FluidPointWindow.hpp"
#include "SolidPointWindow.hpp"

namespace point_time {
    // create fields
    template <class SFPointWindow>
    void createFields(SFPointWindow &pw, const TimeScheme &timeScheme) {
        const std::string &className = bstring::typeName(timeScheme);
        if (className == "NewmarkTimeScheme") {
            NewmarkTimeScheme::createFields<SFPointWindow>(pw);
        } else if (className == "SymplecticTimeScheme") {
            SymplecticTimeScheme::createFields<SFPointWindow>(pw);
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
                point.combineWindows();
                point.separateWindows();
                point.computeStiffToAccel();
                for (const std::shared_ptr<PointWindow> &pw: point.getWindows()) {
                    if (pw->getFluidPointWindow()) NewmarkTimeScheme::update(*pw->getFluidPointWindow(), half, half, half);
                    if (pw->getSolidPointWindow()) NewmarkTimeScheme::update(*pw->getSolidPointWindow(), half, half, half);
                }
            }
            tm.pause();
        } else if (className == "SymplecticTimeScheme") {
            tm.start();
            for (int irep = 0; irep < count; irep++) {
                point.combineWindows();
                point.separateWindows();
                point.computeStiffToAccel();
                for (const std::shared_ptr<PointWindow> &pw: point.getWindows()) {
                    if (pw->getFluidPointWindow()) SymplecticTimeScheme::update(*pw->getFluidPointWindow(), half, half);
                    if (pw->getSolidPointWindow()) SymplecticTimeScheme::update(*pw->getSolidPointWindow(), half, half);
                }
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
