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

class Point;

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
    template <typename Fields>
    static void createFields(Fields &f, int nu_1, int nr, int n_displ) {
        int ndim = (int)f.mStiff.cols();
        f.mStiff.resize(nu_1, ndim);
        f.mStiffR.resize(nr, ndim);
        f.mDispl.resize(n_displ, ndim);
        f.mVeloc.resize(n_displ, ndim);
        f.mStiff.setZero();
        f.mStiffR.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
    }
    
    // update fields on a point
    template <typename Fields>
    static void update(Fields &f,
                       numerical::Real pi_dt, numerical::Real kappa_dt) {
        // update dt
        f.mVeloc += pi_dt * f.mStiffUpdate;
        f.mDispl += kappa_dt * f.mVeloc;
        
        // zero stiffness for next time step
        f.mStiffUpdate.setZero();
    }
    
private:
    // launch fields on points
    void launch(const std::vector<std::shared_ptr<Point>> &points) const;
    // update fields on points
    void update(const std::vector<std::shared_ptr<Point>> &points, int iSubIter) const;
    
    
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
