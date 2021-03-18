//
//  NrFieldConstant.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  constant Nr(s,z)

#ifndef GeneralNrFieldConstant_hpp
#define GeneralNrFieldConstant_hpp

#include "GeneralNrField.hpp"

class GeneralNrFieldConstant: public GeneralNrField {
public:
    // constructor
    GeneralNrFieldConstant(int nr): mNr(nr) {
        // nothing
    }
    
    // get nr by (s, z)
    int getNrAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // the constant nr value
    const int mNr;
};

#endif /* NrFieldConstant_hpp */
