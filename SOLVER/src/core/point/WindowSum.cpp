#include "WindowSum.hpp"

template <>
void SFWindowSum<3>::feedComm(eigen::CColX &buffer, int &row) const {
    if (mStiffR.rows() > 0) {
        buffer.block(row, 0, mStiffR.size(), 1).real() =
        Eigen::Map<const eigen::RColX>(mStiffR.data(),
                                       mStiffR.size());
        row += mStiffR.size();
    } else if (mStiff.rows() > 0) {
         buffer.block(row, 0, mStiff.size(), 1) =
         Eigen::Map<const eigen::CColX>(mStiff.data(),
                                        mStiff.size());
         row += mStiff.size();
    } else {
        throw std::runtime_error("SolidWindowSum::feedComm || Buffer allocation failed.");
    }
};

template <>
void SFWindowSum<1>::feedComm(eigen::CColX &buffer, int &row) const {
    if (mStiffR.rows() > 0) {
        buffer.block(row, 0, mStiffR.rows(), 1).real() = mStiffR;
        row += mStiffR.rows();
    } else if (mStiff.rows() > 0) {
        buffer.block(row, 0, mStiff.rows(), 1) = mStiff;
        row += mStiff.rows();
    } else {
        throw std::runtime_error("FluidWindowSum::feedComm || Buffer allocation failed.");
    }
};

template <>
void SFWindowSum<3>::extractComm(const eigen::CColX &buffer, int &row) {
  if (mStiffR.rows() > 0) {
        mStiffR +=
        Eigen::Map<const CMat>(&buffer(row), mStiffR.rows(), 3).real();
        row += mStiffR.size();
    } else if (mStiff.rows() > 0) {
        mStiff +=
        Eigen::Map<const CMat>(&buffer(row), mStiff.rows(), 3);
        row += mStiff.size();
    } else {
        throw std::runtime_error("SolidWindowSum::extractComm || Buffer allocation failed.");
    }
};

template <>
void SFWindowSum<1>::extractComm(const eigen::CColX &buffer, int &row) {
    if (mStiffR.rows() > 0) {
        mStiffR += buffer.block(row, 0, mStiffR.rows(), 1).real();
        row += mStiffR.rows();
    } else if (mStiff.rows() > 0) {
        mStiff += buffer.block(row, 0, mStiff.rows(), 1);
        row += mStiff.rows();
    } else {
        throw std::runtime_error("FluidWindowSum::extractComm || Buffer allocation failed.");
    }
};
