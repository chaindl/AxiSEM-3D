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
#include "SolidPointWindow.hpp"
class Point;

class ClaytonSolid {
public:
    // constructor
    ClaytonSolid(const std::shared_ptr<SolidPointWindow> &spw):
    mSolidPointWindow(spw) {
        // nothing
    }
    
    // get point
    const std::shared_ptr<const SolidPointWindow> getPointWindow() const {
        return mSolidPointWindow;
    }
    
    // destructor
    virtual ~ClaytonSolid() = default;
    
    // apply ABC
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<SolidPointWindow> mSolidPointWindow;
};

#endif /* ClaytonSolid_hpp */
