//
//  NrFieldPointwise.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  pointwise Nr(s,z)

#ifndef GeneralNrFieldPointwise_hpp
#define GeneralNrFieldPointwise_hpp

#include "GeneralNrField.hpp"
#include "RTreeND.hpp"

class GeneralNrFieldPointwise: public GeneralNrField {
public:
    // constructor
    GeneralNrFieldPointwise(const std::string &fname, double factor,
                     double distTolExact);
    
    // get nr by (s, z)
    int getNrAtPoint(const eigen::DCol2 &sz) const;
    
    // verbose
    std::string verbose() const;
    
private:
    // file name
    const std::string mFilename;
    
    // factor
    const double mFactor;
    
    // dist tolerance for exact match
    const double mDistTolExact;
    
    // rtree
    std::unique_ptr<const RTreeND<2, 1, int>> mRTree = nullptr;
    
    // scanning
    long mSumNrStart = -1;
};

#endif /* GeneralNrFieldPointwise_hpp */
