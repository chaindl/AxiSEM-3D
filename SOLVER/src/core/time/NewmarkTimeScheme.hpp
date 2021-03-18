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
#include "PointWindow.hpp"

class point;

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
    template <class SFPointWindow>
    static void createFields(SFPointWindow &pw) {
        auto &f = pw.getFields();
        int ndim = (int)f.mStiff.cols();
        int nu_1 = pw.getNu_1();
        int nr = pw.getNr();
        f.mStiff.resize(nu_1, ndim);
        f.mStiffR.resize(nr, ndim);
        f.mDispl.resize(nu_1, ndim);
        f.mVeloc.resize(nu_1, ndim);
        f.mAccel.resize(nu_1, ndim);
        f.mStiff.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
        f.mAccel.setZero();
    }
    
    // update fields on a point
    template <class SFPointWindow>
    static void update(SFPointWindow &pw,
                       numerical::Real dt, numerical::Real half_dt,
                       numerical::Real half_dt_dt) {
        auto &f = pw.getFields();
        // update dt
        f.mVeloc += half_dt * (f.mAccel + f.mStiff);
        f.mAccel = f.mStiff;
        f.mDispl += dt * f.mVeloc + half_dt_dt * f.mAccel;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
    }
    
private:
    // update fields on points
    void update(const std::vector<std::shared_ptr<Point>> &points) const;
};

#endif /* NewmarkTimeScheme_hpp */
