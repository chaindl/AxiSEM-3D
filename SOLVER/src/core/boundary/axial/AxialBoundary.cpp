//
//  AxialBoundary.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  axial boundary condition

#include "AxialBoundary.hpp"

// domain
#include "Messaging.hpp"

// apply axial masking
void AxialBoundary::apply() const {
    for (const std::shared_ptr<PointWindow> &pw: mPointWindows) {
        pw->applyAxialBC();
    }
}

// count info
void AxialBoundary::
countInfo(const Messaging &msg, int &solid, int &fluid) const {
    solid = 0;
    fluid = 0;
    for (const auto &pw: mPointWindows) {
        if (!msg.pointInSmallerRank(pw)) {
            if (pw->isFluid()) {
                fluid++;
            } else {
                solid++;
            }
        }
    }
}
