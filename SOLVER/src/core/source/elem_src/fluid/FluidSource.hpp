//
//  FluidSource.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/4/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  source on fluid element

#ifndef FluidSource_hpp
#define FluidSource_hpp

#include "ElementSource.hpp"
class Element;

class FluidSource: public ElementSource {
public:
    // constructor
    FluidSource(std::unique_ptr<STF> &stf,
                const std::shared_ptr<const Element> &element, int m):
    ElementSource(stf), mElement(element), mM(m) {
        // nothing
    }
    
    // destructor
    virtual ~FluidSource() = default;
    
protected:
    // element pointer
    const std::shared_ptr<const Element> mElement;
    const int mM;
};

#endif /* FluidSource_hpp */
