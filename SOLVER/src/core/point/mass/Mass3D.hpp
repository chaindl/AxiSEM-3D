//
//  Mass3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/25/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  3D mass

#ifndef Mass3D_hpp
#define Mass3D_hpp

#include "Mass.hpp"

class Mass3D: public Mass {
public:
    // constructor
    Mass3D(const eigen::DColX &mass):
    mInvMass(mass.cwiseInverse().cast<numerical::Real>()) {
        // nothing
    }
    
    // check compatibility
    void checkCompatibility(int nr, bool solid) const;
    
    // compute accel in-place for fluid
    void computeAccel(eigen::RColX &stiff1) const;
    
    // compute accel in-place for solid
    void computeAccel(eigen::RMatX3 &stiff3) const;
    
    bool is3D() const {return true;};
private:
    // inverse of mass
    const eigen::RColX mInvMass;
};

#endif /* Mass3D_hpp */
