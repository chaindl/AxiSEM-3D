//
//  ClaytonFluid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points

#ifndef ClaytonFluid_hpp
#define ClaytonFluid_hpp

// point
#include <memory>
#include "FluidPointWindow.hpp"
class Point;

class ClaytonFluid {
public:
    // constructor
    ClaytonFluid(const std::shared_ptr<FluidPointWindow> &fpw):
    mFluidPointWindow(fpw) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<const FluidPointWindow> getPointWindow() const {
        return mFluidPointWindow;
    }
    
    // destructor
    virtual ~ClaytonFluid() = default;
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<FluidPointWindow> mFluidPointWindow;
};

#endif /* ClaytonFluid_hpp */
