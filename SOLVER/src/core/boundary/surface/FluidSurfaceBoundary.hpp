//
//  FluidSurfaceBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/30/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  stress-free boundary condition on fluid

#ifndef FluidSurfaceBoundary_hpp
#define FluidSurfaceBoundary_hpp

// point
#include <memory>
#include <vector>
class FluidPointWindow;

// domain
class Messaging;

class FluidSurfaceBoundary {
public:
    // add fluid point
    void addPointWindow(const std::shared_ptr<FluidPointWindow> &fpw) {
        mFluidPointWindows.push_back(fpw);
    }
    
    // apply stress-free boundary condition on fluid
    void apply() const;
    
    // count info
    int countInfo(const Messaging &msg) const;
    
private:
    // points on surface
    std::vector<std::shared_ptr<FluidPointWindow>> mFluidPointWindows;
};

#endif /* FluidSurfaceBoundary_hpp */
