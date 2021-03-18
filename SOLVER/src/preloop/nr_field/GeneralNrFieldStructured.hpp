//
//  GeneralNrFieldStructured.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  structured Nr(s,z)

#ifndef GeneralNrFieldStructured_hpp
#define GeneralNrFieldStructured_hpp

#include "GeneralNrField.hpp"
#include "StructuredGrid.hpp"

class GeneralNrFieldStructured: public GeneralNrField {
public:
    // constructor
    GeneralNrFieldStructured(const std::string &fname, int valOutOfRange);
    
    // get nr by (s, z)
    int getNrAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // file name
    const std::string mFilename;
    
    // factor
    const int mValueOutOfRange;
    
    // grid
    std::unique_ptr<const StructuredGrid<2, int>> mGrid = nullptr;
};

#endif /* NrFieldStructured_hpp */
