// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 20010-2011 Hauke Heibel <hauke.heibel@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Spline_hpp
#define Spline_hpp

#include "eigen_generic.hpp"
#include <iostream>
template <int _Degree, typename _Scalar>
class Spline {
public:
    /** \brief The data type used to store non-zero basis functions. */
    typedef Eigen::Array<_Scalar, 1, _Degree + 1> BasisVectorType;
    typedef Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> RColX;
    typedef Eigen::Matrix<_Scalar, 1, Eigen::Dynamic> RRowX;
    typedef Eigen::Array<_Scalar, Eigen::Dynamic, Eigen::Dynamic> RArrXX;

    Spline() = default;
    
    void expandWorkspace(int dims, int nr = 0) const {
        if (s_ctrls.rows() < dims) s_ctrls.resize(dims, Eigen::NoChange);
        if (s_ctrls.cols() < nr) s_ctrls.resize(Eigen::NoChange, nr);
    }

    void setKnots(const RColX& knots) {m_knots = knots;};
    
    template <typename Mat>
    void setControls(const Mat& ctrls) const {s_ctrls.block(0, 0, ctrls.rows(), ctrls.cols()) = ctrls;};

    const RRowX& knots() const { return m_knots; }
    
    template <typename Mat>
    void getValue(const _Scalar u, Mat &out, const int row, const int dims, const bool accelerate) const;

    int degree() const {return _Degree;};

    int span(const _Scalar u) const;
    
    BasisVectorType BasisFunctions(const _Scalar u) const;

private:
    RRowX m_knots; /*!< Knot vector. */
    inline static RArrXX s_ctrls = RArrXX::Zero(0, 0); /*!< Control points. */
};


template <int _Degree, typename _Scalar>
typename Spline<_Degree, _Scalar>::BasisVectorType
    Spline<_Degree, _Scalar>::BasisFunctions(_Scalar u) const
{
    const int p = degree();
    const int i = span(u);

    const RRowX& U = m_knots;

    BasisVectorType left(p + 1); left(0) = 0.;
    BasisVectorType right(p + 1); right(0) = 0.;

    Eigen::VectorBlock<BasisVectorType, _Degree>(left, 1, p) = u - Eigen::VectorBlock<const RRowX, _Degree>(U, i + 1 - p, p).reverse().array();
    Eigen::VectorBlock<BasisVectorType, _Degree>(right, 1, p) = Eigen::VectorBlock<const RRowX, _Degree>(U, i + 1, p).array() - u;

    BasisVectorType N(1, p + 1);
    N(0) = 1.;
    for (int j = 1; j <= p; ++j)
    {
        _Scalar saved = 0.;
        for (int r = 0; r < j; r++)
        {
            const _Scalar tmp = N(r) / (right(r + 1) + left(j - r));
            N[r] = saved + right(r + 1) * tmp;
            saved = left(j - r) * tmp;
        }
        N(j) = saved;
    }
    return N;
}

template <int _Degree, typename _Scalar>
int Spline<_Degree, _Scalar>::span(const _Scalar u) const
{
      // Piegl & Tiller, "The NURBS Book", A2.1 (p. 68)
    if (u <= m_knots(0)) return degree();
    const _Scalar* pos = std::upper_bound(m_knots.data() + degree() - 1, m_knots.data() + m_knots.size() - degree() - 1, u);
    return static_cast<int>( std::distance(m_knots.data(), pos) - 1 );
}

template <int _Degree, typename _Scalar>
template <typename Mat>
void Spline<_Degree, _Scalar>::getValue(const _Scalar u, Mat &out, const int row, const int dims, const bool accelerate) const
{
  
  static int span_accel;
  static BasisVectorType basis_funcs_accel;
  
  const int p = degree();
  const int span = this->span(u);
  if (!accelerate || span != span_accel) {
      basis_funcs_accel = BasisFunctions(u);
      span_accel = span;
  }
  
  out.block(row, 0, 1, dims).transpose() = (s_ctrls.block(0, span - p, dims, p + 1).rowwise() * basis_funcs_accel).rowwise().sum();
}

#endif // SPLINE_H
