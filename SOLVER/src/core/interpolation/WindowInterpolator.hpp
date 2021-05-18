
#ifndef WindowInterpolator_hpp
#define WindowInterpolator_hpp

#include "SplineFitting.hpp"
#include "numerical.hpp"
#include <iostream>

namespace interpolation {
    const int SplineOrder = 8;
    const int NrExtend = 4;
    const static eigen::RRowX BufferQuery = eigen::RRowX::LinSpaced(SplineOrder, 
            ((double)SplineOrder / 2) / (2 * SplineOrder - 1), (3 * (double)SplineOrder / 2) / (2 * SplineOrder - 1));
};

template <typename _Scalar>
class WindowInterpolator {
public:
    typedef SplineFitting<interpolation::SplineOrder - 1, _Scalar> SplineFittingType;
    typedef Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> RColX;
    typedef Eigen::Matrix<_Scalar,  Eigen::Dynamic, Eigen::Dynamic> RMatXX;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IColX;
    
    WindowInterpolator() = default;
    
    int addSplineFitting(const RColX &knots, int dim) {
        auto mapCheck = mFittingsCollection.insert(std::make_pair(knots, mSplineFittings.size()));
        int tag;
        if (mapCheck.second) {
            mSplineFittings.push_back(std::make_unique<SplineFittingType>(knots.transpose(), dim));
            tag = mSplineFittings.size() - 1;
        } else {
            tag = mapCheck.first->second;
        }

        expandWorkspace(tag, dim, knots.rows());
        return tag;
    }
    
    void expandWorkspace(int tag, int dim, int nr = 0) {
        mSplineFittings[tag]->expandWorkspace(dim);
        if (sInputStorage.rows() < nr) sInputStorage.resize(nr, Eigen::NoChange);
        if (sInputStorage.cols() < dim) sInputStorage.resize(Eigen::NoChange, dim);
    }
    
    void finalize() {mFittingsCollection.clear();};
    
    template <typename Mat1, typename Mat2>
    void interpolatePreloop(int tag, const Mat1 &in, Mat2 &out, const RColX &query_positions, const int row, const int dim) const {
        mSplineFittings[tag]->interpolatePreloop(RMatXX(in));
        for (int i = 0; i < query_positions.cols(); i++) {
            mSplineFittings[tag]->getValue(query_positions(i), out, row + i, dim, (i > 0));
        }
    }
    
    template <typename Mat>
    void newDataWithIndexing(const int tag, const Mat &in, const Eigen::Ref<const IColX>& indices) const {
        for (int i = 0; i < indices.rows(); i++) {
            sInputStorage.block(i, 0, 1, in.cols()) = in.row(indices(i));
        }
        mSplineFittings[tag]->interpolate(sInputStorage.block(0, 0, indices.rows(), in.cols()));
    }
    
    template <typename Mat>
    void newDataAsBlock(const int tag, const Mat &in, int row1, int col1, int nrows, int ncols) const {
        sInputStorage.block(0, 0, nrows, ncols) = in.block(row1, col1, nrows, ncols);
        mSplineFittings[tag]->interpolate(sInputStorage.block(0, 0, nrows, ncols));
    }
    
    template <typename Mat>
    void newDataWithScaling(const int tag, const Mat &in, const Eigen::Ref<const RColX>& winFunc) const {
        for (int i = 0; i < in.cols(); i++) {
            sInputStorage.block(0, i, in.rows(), 1) = in.array().col(i) * winFunc.array();
        }
        mSplineFittings[tag]->interpolate(sInputStorage.block(0, 0, in.rows(), in.cols()));
    }
    
    template <typename Mat>
    void interpolate(const int tag, Mat &out, const Eigen::Ref<const RColX>& query_positions, 
                     const Eigen::Ref<const IColX>& rows, const int dim) const {
        for (int i = 0; i < query_positions.cols(); i++) {
            mSplineFittings[tag]->getValue(query_positions(i), out, rows(i), dim, (i > 0));
        }
    }
    
    template <typename Mat>
    void interpolate(const int tag, Mat &out, const Eigen::Ref<const RColX> &query_positions, 
                     const int row, const int dim) const {
        for (int i = 0; i < query_positions.cols(); i++) {
            mSplineFittings[tag]->getValue(query_positions(i), out, row + i, dim, (i > 0));
        }
    }
    
    static RColX relPhi(const RColX &phi, double phi1, double phi2) {
        if (phi2 - phi1 < numerical::dEpsilon) phi2 += 2 * numerical::dPi;
          
        RColX relPhi = phi;
        for (int alpha = 1; alpha < relPhi.rows(); alpha++) {
            if (relPhi(alpha) < relPhi(alpha - 1)) {
                relPhi.bottomRows(relPhi.rows() - alpha).array() += 2 * numerical::dPi;
                break;
            }
        }
        
        if ((phi2 < relPhi.array()).any() || (relPhi.array() < phi1).any()) {
            throw std::runtime_error("WindowInterpolator::relPhi || phi is out of bounds.");
        }
        
        return (relPhi.array() - phi1) * (1. / (phi2 - phi1));
    }
    
private:
    struct cmp {
        bool operator()(const RColX &a, const RColX &b) const {
            if (a.rows() == b.rows()) {
                return ((a - b).norm() > numerical::epsilon<_Scalar>()) && ((a - b).mean() <= 0);
            } else {
                return a.rows() < b.rows();
            }
        }
    };
  
    std::vector<std::unique_ptr<SplineFittingType>> mSplineFittings;
    std::map<RColX, int, cmp> mFittingsCollection;
    inline static RMatXX sInputStorage = RMatXX::Zero(0,0);
};

#endif
