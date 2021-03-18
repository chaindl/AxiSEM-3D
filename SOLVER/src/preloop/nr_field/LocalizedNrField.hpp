//
//  LocalizedNrField.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#ifndef LocalizedNrField_hpp
#define LocalizedNrField_hpp

#include "eigen_generic.hpp"
#include <memory>

// need distance tolerance for LocalizedNrFieldPointwise
class GeneralNrField;
class LocalizedNr;

namespace eigen {
    // coords
    typedef Eigen::Matrix<double, 2, 1> DCol2;
    typedef Eigen::Matrix<double, Dynamic, 3> DMatX3;
    typedef Eigen::Matrix<int, 1, 4> IRow4;
}

class LocalizedNrField {
public:
    // build from inparam
    static void buildInparam(std::vector<std::unique_ptr<LocalizedNrField>> &localFields, 
        const std::shared_ptr<const GeneralNrField> &gf, const double distTol);
  
    // destructor
    virtual ~LocalizedNrField() = default;
    
    // get nr by (s, z)
    virtual LocalizedNr getWindowsAtPoint(const eigen::DCol2 &sz) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
};

#endif /* LocalizedNrField_hpp */
