//
//  SolidFluidCoupling.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid-fluid boundary condition

#ifndef SolidFluidCoupling_hpp
#define SolidFluidCoupling_hpp

// point
#include "eigen_point.hpp"
#include <memory>
class PointWindow;

class SolidFluidCoupling {
public:
    // constructor
    SolidFluidCoupling(const std::shared_ptr<PointWindow> &spw,
                       const std::shared_ptr<PointWindow> &fpw);
    
    // destructor
    virtual ~SolidFluidCoupling() = default;
    
    // get solid point
    const std::shared_ptr<PointWindow> &getSolidPointWindow() const {
        return mSolidPointWindow;
    }
    
    // apply coupling
    void apply() const;
    
    
    ////////////////////////////// virtual //////////////////////////////
protected:
    // check compatibility
    virtual void checkCompatibility(int nr) const;
    
public:
    // solid => fluid
    virtual void coupleSolidToFluid(const eigen::CMatX3 &solidDispl,
                                    eigen::CColX &fluidStiff) const = 0;
    virtual void coupleSolidToFluid(const eigen::RMatX3 &solidDispl,
                                    eigen::RColX &fluidStiff) const = 0;
    // fluid => solid
    virtual void coupleFluidToSolid(const eigen::CColX &fluidStiff,
                                    eigen::CMatX3 &solidStiff) const = 0;
    virtual void coupleFluidToSolid(const eigen::RColX &fluidStiff,
                                    eigen::RMatX3 &solidStiff) const = 0;
protected:
    // coupled solid-fluid pair
    const std::shared_ptr<PointWindow> mSolidPointWindow;
    const std::shared_ptr<PointWindow> mFluidPointWindow;
};

#endif /* SolidFluidCoupling_hpp */
