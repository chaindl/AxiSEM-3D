//
//  LocalizedNrField.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#ifndef LocalizedNrFieldPointwise_hpp
#define LocalizedNrFieldPointwise_hpp

#include "LocalizedNrField.hpp"

class LocalizedNr;

class LocalizedNrFieldPointwise: public LocalizedNrField {
public:
    LocalizedNrFieldPointwise(const std::string &fname, double distTol);
  
    // get nr by (s, z)
    LocalizedNr getWindowsAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
};

#endif /* LocalizedNrFieldPointwise_hpp */
