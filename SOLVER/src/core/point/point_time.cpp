//
//  point_time.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/11/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  time-scheme-dependent point methods

#include "NewmarkTimeScheme.hpp"
#include "SymplecticTimeScheme.hpp"
#include "timer.hpp"
#include "bstring.hpp"
#include "Point.hpp"

namespace point_time {
    // create fields
    void createFields(PointWindow &pw, const TimeScheme &timeScheme) {
        const std::string &className = bstring::typeName(timeScheme);
        if (className == "NewmarkTimeScheme") {
            if (pw.isFluid()) {
                NewmarkTimeScheme::createFields(pw.getFluidFields(), pw.getNu_1(), pw.getNr());
            } else {
                NewmarkTimeScheme::createFields(pw.getSolidFields(), pw.getNu_1(), pw.getNr());
            }
        } else if (className == "SymplecticTimeScheme") {
            if (pw.isFluid()) {
                SymplecticTimeScheme::createFields(pw.getFluidFields(), pw.getNu_1(), pw.getNr());
            } else {
                SymplecticTimeScheme::createFields(pw.getSolidFields(), pw.getNu_1(), pw.getNr());
            }
        } else {
            throw std::runtime_error("point_time::createFields || "
                                     "Unknown derived class of TimeScheme: "
                                     + className);
        }
    }
    
    // measure cost of a point
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
                    if (pw->isFluid()) {
                        NewmarkTimeScheme::update(pw->getFluidFields(), half, half, half);
                    } else {
                        NewmarkTimeScheme::update(pw->getSolidFields(), half, half, half);
                    }
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
                    if (pw->isFluid()) {
                        SymplecticTimeScheme::update(pw->getFluidFields(), half, half);
                    } else {
                        SymplecticTimeScheme::update(pw->getSolidFields(), half, half);
                    }
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
