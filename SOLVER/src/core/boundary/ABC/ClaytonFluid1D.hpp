//
//  ClaytonFluid1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for fluid points in 1D

#ifndef ClaytonFluid1D_hpp
#define ClaytonFluid1D_hpp

#include "ClaytonFluid.hpp"
#include "numerical.hpp"
#include "PointWindow.hpp"

template <typename T>
struct InterfaceRC_ClaytonFluid {
};

template <>
struct InterfaceRC_ClaytonFluid<numerical::Real> {
    inline static auto &getVeloc(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsR().mVeloc;};
    inline static auto &getStiff(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsR().mStiffR;};
};

template <>
struct InterfaceRC_ClaytonFluid<numerical::ComplexR> {
    inline static auto &getVeloc(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsC().mVeloc;};
    inline static auto &getStiff(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsC().mStiff;};
};


template<typename T>
class ClaytonFluid1D: public ClaytonFluid {
public:
    // constructor
    ClaytonFluid1D(const std::shared_ptr<PointWindow> &fpw,
                   double rhoVp, double area):
    ClaytonFluid(fpw), mAreaOverRhoVp(area / rhoVp) {
        // nothing
    }
    
    // apply ABC
    void apply() const;
    
private:
    // area / (rho * vp)
    const numerical::Real mAreaOverRhoVp;
};

template<typename T>
void ClaytonFluid1D<T>::apply() const {
    // get fields
    Eigen::Matrix<T, Eigen::Dynamic, 1> &veloc = InterfaceRC_ClaytonFluid<T>::getVeloc(mPointWindow);
    Eigen::Matrix<T, Eigen::Dynamic, 1> &stiff = InterfaceRC_ClaytonFluid<T>::getStiff(mPointWindow);
    
    // apply
    stiff -= veloc * mAreaOverRhoVp;
}

#endif /* ClaytonFluid1D_hpp */
