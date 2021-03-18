//
//  LocalizedNrField.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#ifndef LocalizedNrFieldRefinement_hpp
#define LocalizedNrFieldRefinement_hpp

#include "LocalizedNrField.hpp"

class LocalizedNr;
class GeneralNrField;

class LocalizedNrFieldRefinement: public LocalizedNrField {
public:
    LocalizedNrFieldRefinement(const std::shared_ptr<const GeneralNrField> &gnr, 
    eigen::DMatX3 &coords, double &value, bool relative, bool geographic, bool ellipticity);
  
    // get nr by (s, z)
    LocalizedNr getWindowsAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    const std::shared_ptr<const GeneralNrField> mGeneralNrField;
    const eigen::DMatX3 mCoords;
    const double mValue; 
    bool mIsRelative, mGeographic, mEllipticity;
};

#endif /* LocalizedNrFieldRefinement_hpp */
