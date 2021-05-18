//
//  ClaytonSolid.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points

#ifndef ClaytonSolid_hpp
#define ClaytonSolid_hpp

// point
#include <memory>
#include "PointWindow.hpp"
class Point;

class ClaytonSolid {
public:
    // constructor
    ClaytonSolid(const std::shared_ptr<PointWindow> &spw):
    mPointWindow(spw) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<const PointWindow> getPointWindow() const {
        return mPointWindow;
    }
    
    // destructor
    virtual ~ClaytonSolid() = default;
    
    virtual void checkCompatibility() const {};
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<PointWindow> mPointWindow;
};

#endif /* ClaytonSolid_hpp */
