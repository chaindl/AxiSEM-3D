//
//  AbsorbingBoundary.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/1/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  absorbing boundary condition

#ifndef AbsorbingBoundary_hpp
#define AbsorbingBoundary_hpp

#include <vector>
#include "ClaytonSolid.hpp"
#include "ClaytonFluid.hpp"

#include "Sponge.hpp"
// domain
#include <map>
class Messaging;

class AbsorbingBoundary {
public:
    // add Clayton in solid
    void addClaytonSolid(std::unique_ptr<const ClaytonSolid> &clayton) {
        clayton->checkCompatibility();
        mClaytonSolids.push_back(std::move(clayton));
    }
    
    // add Clayton in fluid
    void addClaytonFluid(std::unique_ptr<const ClaytonFluid> &clayton) {
        clayton->checkCompatibility();
        mClaytonFluids.push_back(std::move(clayton));
    }
    
    // add sponge in solid
    void addSpongeSolid(std::unique_ptr<const Sponge<3>> &sponge) {
        mSpongeSolids.push_back(std::move(sponge));
    }
    
    // add sponge in fluid
    void addSpongeFluid(std::unique_ptr<const Sponge<1>> &sponge) {
        mSpongeFluids.push_back(std::move(sponge));
    }
    
    // apply Clayton ABC
    void applyClayton() const {
        for (const auto &clayton: mClaytonSolids) {
            clayton->apply();
        }
        for (const auto &clayton: mClaytonFluids) {
            clayton->apply();
        }
    }
    
    // apply sponge ABC
    void applySponge() const {
        for (const auto &sponge: mSpongeSolids) {
            sponge->apply();
        }
        for (const auto &sponge: mSpongeFluids) {
            sponge->apply();
        }
    }
    
    // count info
    std::map<std::string, int> countInfo(const Messaging &msg) const;
    
private:
    // Claytons
    std::vector<std::unique_ptr<const ClaytonSolid>> mClaytonSolids;
    std::vector<std::unique_ptr<const ClaytonFluid>> mClaytonFluids;
    
    // sponges
    std::vector<std::unique_ptr<const Sponge<3>>> mSpongeSolids;
    std::vector<std::unique_ptr<const Sponge<1>>> mSpongeFluids;
};

#endif /* AbsorbingBoundary_hpp */
