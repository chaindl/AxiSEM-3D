//
//  ClaytonFluid3D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 3D

#ifndef ClaytonFluid3D_hpp
#define ClaytonFluid3D_hpp

#include "ClaytonFluid.hpp"
#include "eigen_point.hpp"

class ClaytonFluid3D: public ClaytonFluid {
public:
    // constructor
    ClaytonFluid3D(const std::shared_ptr<PointWindow> &fpw,
                   const eigen::DColX &rhoVp, const eigen::DColX &area):
    ClaytonFluid(fpw),
    mAreaOverRhoVp(area.cwiseQuotient(rhoVp).cast<numerical::Real>()) {
    };
    
public:
  // check compatibility
    virtual void checkCompatibility() const = 0;
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // area / (rho * vp)
    const eigen::RColX mAreaOverRhoVp;
    
};

class ClaytonFluid3D_C: public ClaytonFluid3D {
public:
    // constructor
    using ClaytonFluid3D::ClaytonFluid3D;
    // apply ABC
    void apply() const;
    
    void checkCompatibility() const;
    
    ////////////////////////////////////////
    //////////////// static ////////////////
    ////////////////////////////////////////
    
    // workspace
    inline static eigen::RColX sVecR = eigen::RColX(0);
    inline static eigen::CColX sVecC = eigen::CColX(0);
};

class ClaytonFluid3D_R: public ClaytonFluid3D {
public:
    // constructor
    using ClaytonFluid3D::ClaytonFluid3D;
    // apply ABC
    void apply() const;
    
    void checkCompatibility() const;    
};

#endif /* ClaytonFluid3D_hpp */
