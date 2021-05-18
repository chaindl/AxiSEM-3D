//
//  ClaytonSolid1D.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 7/29/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Clayton-Enquist ABC for solid points in 1D
//  theta: angle between surface normal and z-axis

#ifndef ClaytonSolid1D_hpp
#define ClaytonSolid1D_hpp

#include "ClaytonSolid.hpp"
#include "numerical.hpp"
#include "PointWindow.hpp"

template <typename T>
struct InterfaceRC_ClaytonSolid {
};

template <>
struct InterfaceRC_ClaytonSolid<numerical::Real> {
    inline static auto &getVeloc(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsR().mVeloc;};
    inline static auto &getStiff(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsR().mStiffR;};
};

template <>
struct InterfaceRC_ClaytonSolid<numerical::ComplexR> {
    inline static auto &getVeloc(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsC().mVeloc;};
    inline static auto &getStiff(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsC().mStiff;};
};

template<typename T>
class ClaytonSolid1D: public ClaytonSolid {
public:
    // constructor
    ClaytonSolid1D(const std::shared_ptr<PointWindow> &spw,
                   double rhoVp, double rhoVs, double area, double theta):
    ClaytonSolid(spw),
    mRSA_CosT2_p_RPA_SinT2(rhoVs * area * cos(theta) * cos(theta) +
                           rhoVp * area * sin(theta) * sin(theta)),
    mRSA_SinT2_p_RPA_CosT2(rhoVs * area * sin(theta) * sin(theta) +
                           rhoVp * area * cos(theta) * cos(theta)),
    mRPA_m_RSA_x_CosT_SinT((rhoVp * area - rhoVs * area) *
                           cos(theta) * sin(theta)),
    mRSA(rhoVs * area) {
        // nothing
    }
    
    // apply ABC
    void apply() const;
    
private:
    // RSA = rho * vs * area
    // RPA = rho * vp * area
    // RSA Cos[t]^2 + RPA Sin[t]^2
    const numerical::Real mRSA_CosT2_p_RPA_SinT2;
    // RSA Sin[t]^2 + RPA Cos[t]^2
    const numerical::Real mRSA_SinT2_p_RPA_CosT2;
    // (RPA - RSA) Cos[t] Sin[t]
    const numerical::Real mRPA_m_RSA_x_CosT_SinT;
    // RSA (for transverse component)
    const numerical::Real mRSA;
};

// apply ABC
template<typename T>
void ClaytonSolid1D<T>::apply() const {
    // get fields
    const Eigen::Matrix<T, Eigen::Dynamic, 3> &veloc = InterfaceRC_ClaytonSolid<T>::getVeloc(mPointWindow);
    Eigen::Matrix<T, Eigen::Dynamic, 3> &stiff = InterfaceRC_ClaytonSolid<T>::getStiff(mPointWindow);
    
    // s, z
    stiff.col(0) -= (mRSA_CosT2_p_RPA_SinT2 * veloc.col(0) +
                     mRPA_m_RSA_x_CosT_SinT * veloc.col(2));
    stiff.col(2) -= (mRPA_m_RSA_x_CosT_SinT * veloc.col(0) +
                     mRSA_SinT2_p_RPA_CosT2 * veloc.col(2));
    
    // phi
    stiff.col(1) -= mRSA * veloc.col(1);
}


#endif /* ClaytonSolid1D_hpp */
