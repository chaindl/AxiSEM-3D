//
//  GeneralNrFieldAnalytical.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  analytical Nr(s,z)

#ifndef GeneralNrFieldAnalytical_hpp
#define GeneralNrFieldAnalytical_hpp

#include "GeneralNrField.hpp"
#include <vector>

class GeneralNrFieldAnalytical: public GeneralNrField {
public:
    // constructor
    GeneralNrFieldAnalytical();
    
    // get nr by (s, z)
    int getNrAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // code ID
    static std::string sCodeID;
    
    // TODO: add your data here
    // below are data for
    // sCodeID = "depth-dependent (AxiSEM3D default)"
    std::vector<double> mControlDepths;
    std::vector<double> mControlNrs;
};

#endif /* GeneralNrFieldAnalytical_hpp */
