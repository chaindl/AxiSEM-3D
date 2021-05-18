#include "SFRCPointWindow.hpp"

template <>
const Fields<3, numerical::Real> &SolidRCPointWindow<3, numerical::Real>::getSolidFieldsR() const {
    return this->mFields;
};

template <>
Fields<3, numerical::Real> &SolidRCPointWindow<3, numerical::Real>::getSolidFieldsR() {
    return this->mFields;
};

template <>
const Fields<1, numerical::Real> &FluidRCPointWindow<1, numerical::Real>::getFluidFieldsR() const {
    return this->mFields;
};

template <>
Fields<1, numerical::Real> &FluidRCPointWindow<1, numerical::Real>::getFluidFieldsR() {
    return this->mFields;
};

template <>
const Fields<3, numerical::ComplexR> &SolidRCPointWindow<3, numerical::ComplexR>::getSolidFieldsC() const {
    return this->mFields;
};

template <>
Fields<3, numerical::ComplexR> &SolidRCPointWindow<3, numerical::ComplexR>::getSolidFieldsC() {
    return this->mFields;
};

template <>
const Fields<1, numerical::ComplexR> &FluidRCPointWindow<1, numerical::ComplexR>::getFluidFieldsC() const {
    return this->mFields;
};

template <>
Fields<1, numerical::ComplexR> &FluidRCPointWindow<1, numerical::ComplexR>::getFluidFieldsC() {
    return this->mFields;
};
