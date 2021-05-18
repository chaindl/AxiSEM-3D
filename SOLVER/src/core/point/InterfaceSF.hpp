#ifndef InterfaceSF_hpp
#define InterfaceSF_hpp

#include "eigen_generic.hpp"
#include "fft.hpp"

template <int dims>
struct InterfaceSF {
};

template <>
struct InterfaceSF<3> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 3> RMatSF;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> CMatSF;
    inline static SolverFFTW<numerical::Real, 3> &fftSF = fft::gFFT_3;
    inline static auto &getFieldsR(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsR();};
    inline static auto &getFieldsC(const std::shared_ptr<PointWindow> pw) {return pw->getSolidFieldsC();};
    inline static auto &getStiffForWindowSum(const std::shared_ptr<PointWindow> pw) {return pw->getSolidStiffForWindowSum();};
};

template <>
struct InterfaceSF<1> {
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> RMatSF;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> CMatSF;
    inline static SolverFFTW<numerical::Real, 1> &fftSF = fft::gFFT_1;
    inline static auto &getFieldsR(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsR();};
    inline static auto &getFieldsC(const std::shared_ptr<PointWindow> pw) {return pw->getFluidFieldsC();};
    inline static auto &getStiffForWindowSum(const std::shared_ptr<PointWindow> pw) {return pw->getFluidStiffForWindowSum();};
};

#endif
