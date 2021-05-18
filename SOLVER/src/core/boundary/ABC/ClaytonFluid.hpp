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
#include "PointWindow.hpp"
class Point;

class ClaytonFluid {
public:
    // constructor
    ClaytonFluid(const std::shared_ptr<PointWindow> &fpw):
    mPointWindow(fpw) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<const PointWindow> getPointWindow() const {
        return mPointWindow;
    }
    
    // destructor
    virtual ~ClaytonFluid() = default;
    
    virtual void checkCompatibility() const {};
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<PointWindow> mPointWindow;
};

#endif /* ClaytonFluid_hpp */
