//
//  AxialBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  axial boundary condition

#ifndef AxialBoundary_hpp
#define AxialBoundary_hpp

// point
#include <memory>
#include <vector>
class PointWindow;

// domain
class Messaging;

class AxialBoundary {
public:
    // add solid point
    void addPointWindow(const std::shared_ptr<PointWindow> &pw) {
        mPointWindows.push_back(pw);
    }
    
    // apply axial masking
    void apply() const;
    
    // count info
    void countInfo(const Messaging &msg, int &solid, int &fluid) const;
    
private:
    // points on axis
    std::vector<std::shared_ptr<PointWindow>> mPointWindows;
};

#endif /* AxialBoundary_hpp */
