// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 20010-2011 Hauke Heibel <hauke.heibel@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SplineFitting_hpp
#define SplineFitting_hpp

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "Spline.hpp"
#include <iostream>
template <int _Degree, typename _Scalar>
class SplineFitting {
public:
    typedef Eigen::Matrix<_Scalar, 1, Eigen::Dynamic> RRowX;
    typedef Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> RMatXX;
    
    SplineFitting(const RRowX& knot_parameters, const int dim):
        mQr(knot_parameters.cols(), knot_parameters.cols()) {
        
        RRowX knots;
        int nr = knot_parameters.cols();
        KnotAveraging(knot_parameters, mSpline.degree(), knots);
        
        mSpline.setKnots(knots);
        mSpline.expandWorkspace(dim, nr);
        
        RMatXX A = RMatXX::Zero(nr, nr);
        for (int i = 1; i < nr - 1; ++i)
        {
          const int span = mSpline.span(knot_parameters(i));
          
          // The segment call should somehow be told the spline order at compile time.
          A.row(i).segment(span - mSpline.degree(), mSpline.degree() + 1) = mSpline.BasisFunctions(knot_parameters(i));
        }

        A(0, 0) = 1.0;
        A(nr - 1, nr - 1) = 1.0;

        mQr.compute(A);
    }
    
    void expandWorkspace(int dim) {
        mSpline.expandWorkspace(dim);
    }
    
    void interpolatePreloop(const RMatXX& pts_in) const {
        mSpline.setControls(mQr.solve(pts_in).transpose());
    }
    
    void interpolate(const Eigen::Ref<RMatXX>& pts_in) const {
        Eigen::internal::set_is_malloc_allowed(true);
        mSpline.setControls(mQr.solve(pts_in).transpose());
        Eigen::internal::set_is_malloc_allowed(false);
    }
    
    template <typename Mat>
    void getValue(const _Scalar u, Mat &out, const int row, 
                  const int dims, const bool accelerate) const {        
        mSpline.getValue(u, out, row, dims, accelerate);
    }
    
private:
    void KnotAveraging(const RRowX& parameters, int degree, RRowX& knots)
    {
        knots.resize(parameters.size() + degree + 1);      
        for (int j = 1; j < parameters.size() - degree; ++j)
          knots(j + degree) = parameters.segment(j, degree).mean();

        knots.segment(0, degree + 1) = RRowX::Zero(degree + 1);
        knots.segment(knots.size() - degree - 1, degree + 1) = RRowX::Constant(degree + 1, 1);
    }
    
    Eigen::HouseholderQR<RMatXX> mQr;
    Spline<_Degree, _Scalar> mSpline;
};

#endif // SplineFitting_hpp
