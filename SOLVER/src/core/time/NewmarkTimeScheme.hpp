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
#include <iostream>
class Point;

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
    template <typename Fields>
    static void createFields(Fields &f, int nu_1, int nr, int n_displ) {
        int ndim = (int)f.mStiff.cols();
        f.mStiff.resize(nu_1, ndim);
        f.mStiffR.resize(nr, ndim);
        f.mDispl.resize(n_displ, ndim);
        f.mVeloc.resize(n_displ, ndim);
        f.mAccel.resize(n_displ, ndim);   
        f.mStiff.setZero();
        f.mStiffR.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
        f.mAccel.setZero();
    }
    
    // update fields on a point
    template <typename Fields>
    static void update(Fields &f,
                       numerical::Real dt, numerical::Real half_dt,
                       numerical::Real half_dt_dt) {
        // update dt
        f.mVeloc += half_dt * (f.mAccel + f.mStiffUpdate);
        f.mAccel = f.mStiffUpdate;
        f.mDispl += dt * f.mVeloc + half_dt_dt * f.mAccel;
        
        // zero stiffness for next time step
        f.mStiffUpdate.setZero();
    }
    
private:
    // update fields on points
    void update(const std::vector<std::shared_ptr<Point>> &points) const;
};

#endif /* NewmarkTimeScheme_hpp */
