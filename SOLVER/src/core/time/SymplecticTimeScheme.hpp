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
#include <memory>

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
        f.mStiff.setZero();
        f.mDispl.setZero();
        f.mVeloc.setZero();
    }
    
    // update fields on a point
    template <class SFPointWindow>
    static void update(SFPointWindow &pw,
                       numerical::Real pi_dt, numerical::Real kappa_dt) {
        auto &f = pw.getFields();
        // update dt
        f.mVeloc += pi_dt * f.mStiff;
        f.mDispl += kappa_dt * f.mVeloc;
        
        // zero stiffness for next time step
        f.mStiff.setZero();
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
