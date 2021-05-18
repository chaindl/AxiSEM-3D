#ifndef WindowSum_hpp
#define WindowSum_hpp

#include "eigen_point.hpp"
#include "PointWindow.hpp"
#include "WindowInterpolator.hpp"
#include "InterfaceSF.hpp"
#include <iostream>

using interpolation::NrExtend;

class WindowSum {
public:
    // constructor
    WindowSum(const std::shared_ptr<WindowInterpolator<numerical::Real>> &interpolator, int interpTagWhole, int nr):
    mInterpolator(interpolator), mInterpTagWhole(interpTagWhole), mNr(nr) {
        mPosIdxWhole.resize(mNr + 2 * NrExtend);
        mPosIdxWhole << eigen::IColX::LinSpaced(NrExtend, mNr - NrExtend, mNr - 1),
                        eigen::IColX::LinSpaced(mNr, 0, mNr - 1),
                        eigen::IColX::LinSpaced(NrExtend, 0, NrExtend - 1);
    }
    
    // destructor
    virtual ~WindowSum() = default;
    
    virtual bool onlyOneWindow() const = 0;
    virtual bool needsFT() const = 0;
    virtual void allocateComm() = 0;
    virtual int sizeComm() const = 0;
    
    virtual void setInFourier() = 0;
    
    // feed to mpi buffer
    virtual void feedComm(eigen::CColX &buffer, int &row) const = 0;

    // extract from mpi buffer
    virtual void extractComm(const eigen::CColX &buffer, int &row) = 0;
    
    virtual void overlapAndAddStiff() = 0;
        
    virtual void scatterStiffToWindows() = 0;

    virtual int getTag() const = 0;
        
    void addWindow(std::shared_ptr<PointWindow> &win, int interpTag, 
                   eigen::IColX posIdx, eigen::RColX relPhi) {
        checkCompatibility(win->isFluid(), posIdx.size() > 0);
        mWindows.push_back(win);
        mInterpTags.push_back(interpTag);
        mPosIdx.push_back(posIdx);
        mRelPhis.push_back(relPhi);
    }
    
    virtual void checkCompatibility(bool fluid, bool needRealStorage) = 0;
    
protected:
    const int mNr;
    const std::shared_ptr<WindowInterpolator<numerical::Real>> mInterpolator;
    std::vector<std::shared_ptr<PointWindow>> mWindows;
    const int mInterpTagWhole;
    eigen::IColX mPosIdxWhole;
    std::vector<int> mInterpTags;
    std::vector<eigen::IColX> mPosIdx;
    std::vector<eigen::RColX> mRelPhis;
};

template <int dims>
class SFWindowSum: public WindowSum {
public:
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, dims> RMat;
    typedef Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, dims> CMat;
    
    using WindowSum::WindowSum;
    
    int getTag() const {return mWindows[0]->getMeshTag();};
    
    void setInFourier() {
        for (std::shared_ptr<PointWindow> &win: mWindows) {
            win->setInFourier(!needsFT());
        }
    }
    
    bool onlyOneWindow() const {return (mWindows.size() == 1);};
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    bool needsFT() const {
        return (mStiffR.size() > 0 || mWindows[0]->is3D());
    }
    
    void allocateComm() {
        if (needsFT()) {
            mStiffR = RMat::Zero(mNr, dims);
        } else {
            mStiff = CMat::Zero(mWindows[0]->getNu_1(), dims);
        }
    }
    
    int sizeComm() const {
        if (mStiffR.size() > 0) {
            return (int)mStiffR.size();
        } else {
            return (int)mStiff.size();
        }
    }
    
    // feed to mpi buffer
    void feedComm(eigen::CColX &buffer, int &row) const;
    
    // extract from mpi buffer
    void extractComm(const eigen::CColX &buffer, int &row);
    
    void overlapAndAddStiff() {
        if (mWindows.size() > 1) {
            mStiffR.setZero();
            for (int m = 0; m < mWindows.size(); m++) {
                if (mInterpTags[m] >= 0) {
                    mInterpolator->newDataWithScaling(mInterpTags[m], InterfaceSF<dims>::getFieldsR(mWindows[m]).mStiffR,
                                   mWindows[m]->getFracForWindowSum());
                    mInterpolator->interpolate(mInterpTags[m], mStiffR, mRelPhis[m], mPosIdx[m], dims);
                } else {
                    for (int i = 0; i < mPosIdx[m].rows(); i++) {
                        mStiffR.row(mPosIdx[m](i)) += InterfaceSF<dims>::getStiffForWindowSum(mWindows[m]).row(i) * 
                                                      mWindows[m]->getFracForWindowSum()(i);
                    }
                }
            }
        } else if (mStiffR.rows() > 0) {
            mStiff = InterfaceSF<dims>::getStiffForWindowSum(mWindows[0]);
        } else if (mStiff.rows() > 0) {
            mStiff = InterfaceSF<dims>::getFieldsC(mWindows[0]).mStiff;
        }
    }
    
    void scatterStiffToWindows() {
        if (mWindows.size() > 1) {
            if (mInterpTagWhole >= 0) mInterpolator->newDataWithIndexing(mInterpTagWhole, mStiffR, mPosIdxWhole);
            for (int m = 0; m < mWindows.size(); m++) {
                if (mInterpTags[m] >= 0) {
                    mInterpolator->interpolate(mInterpTagWhole, InterfaceSF<dims>::getStiffForWindowSum(mWindows[m]),
                            mWindows[m]->getPhiForWindowSum(), 0, dims);
                } else {
                    mWindows[m]->collectStiffFromWindowSum(mStiffR, mPosIdx[m]);
                }
            }
        } else if (mStiffR.rows() > 0) {
            mWindows[0]->collectStiffFromWindowSum(mStiffR);
        } else if (mStiff.size() > 0) {
            mWindows[0]->collectStiffFromWindowSum(mStiff);
        }
    }
    
    void checkCompatibility(bool fluid, bool needRealStorage) {
        if (fluid && dims!=1) {
            throw std::runtime_error("SFWindowSum::addWindow || adding fluid window to solid window sum.");
        } 
        if (!fluid && dims!=3) {
            throw std::runtime_error("SFWindowSum::addWindow || adding solid window to fluid window sum.");
        }
        
        if (needRealStorage) mStiffR.resize(mNr, dims);
    }
    
protected:
    RMat mStiffR = RMat::Zero(0, dims);
    CMat mStiff = CMat::Zero(0, dims);
};

#endif
