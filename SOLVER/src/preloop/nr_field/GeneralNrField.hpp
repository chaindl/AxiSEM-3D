//
//  GeneralNrField.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#ifndef GeneralNrField_hpp
#define GeneralNrField_hpp

#include "eigen_generic.hpp"
#include <memory>

namespace eigen {
    // coords
    typedef Eigen::Matrix<double, 2, 1> DCol2;
}

class GeneralNrField {
public:
    // build from inparam
    static std::shared_ptr<const GeneralNrField>
    buildInparam(const double distTolerance);
  
    // destructor
    virtual ~GeneralNrField() = default;
    
    // get nr by (s, z)
    virtual int getNrAtPoint(const eigen::DCol2 &sz) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
};

#endif /* GeneralNrField_hpp */
