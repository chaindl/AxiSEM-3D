#ifndef WindowSum_hpp
#define WindowSum_hpp

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include "numerical.hpp"
#include "eigen_point.hpp"
#include "SolidPointWindow.hpp"
#include "FluidPointWindow.hpp"
#include <iostream>
class WindowSum {
public:
    // constructor
    WindowSum(const eigen::RColX phi):
    mPhi(phi) {
        // nothing
    }
    
    // destructor
    virtual ~WindowSum() = default;
    
    const eigen::RColX getPhi() const {return mPhi;};
    
    virtual bool onlyOneWindow() const = 0;
    virtual bool needsFT() const = 0;
    virtual void allocateComm() = 0;
    virtual int sizeComm() const = 0;
    
    virtual void setInFourier() = 0;
    
    // feed to mpi buffer
    virtual void feedComm(eigen::RColX &buffer, int &row) const = 0;

    // extract from mpi buffer
    virtual void extractComm(const eigen::RColX &buffer, int &row) = 0;
    
    virtual void overlapAndAddStiff() = 0;
        
    virtual void scatterStiffToWindows() = 0;
    
    eigen::RColX interpolateWinToWhole(const eigen::RColX &phi_whole, const eigen::RColX &fun, const eigen::RColX &phi_win) const {
        double dphi = phi_win(1) - phi_win(0);
        if (dphi < 0) dphi += 2 * numerical::dPi;
      
        boost::math::interpolators::cardinal_cubic_b_spline<numerical::Real> spline(fun.begin(), fun.end(), phi_win(0), dphi);
        
        eigen::RColX fun_whole = eigen::RColX::Zero(phi_whole.rows());
        for (int i = 0; i < phi_whole.rows(); i++) {
            double phi_with_wrapping = phi_whole(i);
            if (phi_win(phi_win.size() - 1) < phi_win(0) && phi_with_wrapping < phi_win(0)) {
                phi_with_wrapping += 2 * numerical::dPi;
            }
            if (phi_with_wrapping <= phi_win(phi_win.size() - 1) + 2 * numerical::dPi) {
                fun_whole(i) = spline(phi_with_wrapping);  
            }
        }
        return fun_whole;
    };

    eigen::RColX interpolateWholeToWin(const eigen::RColX &phi_win, const eigen::RColX &fun, const eigen::RColX &phi_whole) const {
        double dphi = phi_whole(1) - phi_whole(0);
        if (dphi < 0) dphi += 2 * numerical::dPi;
        
        // periodic bc
        double dfun_dphi_left = 0.5 * (fun(fun.size() - 1) - fun(1)) / dphi;
        double dfun_dphi_right = 0.5 * (fun(fun.size() - 2) - fun(0)) / dphi;
        boost::math::interpolators::cardinal_cubic_b_spline<numerical::Real> spline(fun.begin(), fun.end(), phi_whole(0), dphi, dfun_dphi_left, dfun_dphi_right);
        
        eigen::RColX fun_win = eigen::RColX::Zero(phi_win.rows());
        for (int i = 0; i < phi_win.rows(); i++) {
            fun_win(i) = spline(phi_win(i));  
        }
        return fun_win;
    };
    
    static void addSolidWindow(const std::shared_ptr<SolidPointWindow> &win, const int aligned) {
        throw std::runtime_error("WindowSum::addSolidWindow || attempting to add solid point window to fluid window sum.");
    };
    
    static void addFluidWindow(const std::shared_ptr<FluidPointWindow> &win, const int aligned) {
        throw std::runtime_error("WindowSum::addFluidWindow || attempting to add fluid point window to solid window sum.");
    };
    
protected:
    const eigen::RColX mPhi;
};

class SolidWindowSum: public WindowSum {
public:
    SolidWindowSum(const eigen::RColX phi, bool onlyOneWindow):
    WindowSum(phi) {
        if (!onlyOneWindow) {
            mStiffR3 = eigen::RMatX3::Zero(phi.rows(), 3);
        }
    }
    
    void setInFourier() {
        for (auto &win: mWindows) {
            win->setInFourier(!needsFT());
        }
    }
    
    bool onlyOneWindow() const {return (mWindows.size() == 1);};
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    bool needsFT() const {
        return (mWindows.size() > 1 || mWindows[0]->is3D());
    }
    
    void allocateComm() {
        //if (needsFT()) {
            mStiffR3 = eigen::RMatX3::Zero(mPhi.rows(), 3);
            for (auto &win: mWindows) {
                win->setInFourier(false);
            }
        //} else {
        //    mStiff3 = eigen::CMatX3::Zero(mWindows[0]->getNu_1(), 3);
        //}
    }
    
    int sizeComm() const {
        if (mStiffR3.size() > 0) {
            return (int)mStiffR3.size();
        } else {
            return (int)mStiff3.size();
        }
    }
    
    // feed to mpi buffer
    void feedComm(eigen::RColX &buffer, int &row) const {
        if (mStiffR3.rows() > 0) {
            buffer.block(row, 0, mStiffR3.size(), 1) =
            Eigen::Map<const eigen::RColX>(mStiffR3.data(),
                                           mStiffR3.size());
            row += mStiffR3.size();
        //} else if (mStiff3.rows() > 0) {
            // buffer.block(row, 0, mStiff3.size(), 1) =
            // Eigen::Map<const eigen::CColX>(mStiff3.data(),
            //                                mStiff3.size());
            // row += mStiff3.size();
        } else {
            throw std::runtime_error("SolidWindowSum::feedComm || Buffer allocation failed.");
        }
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::RColX &buffer, int &row) {
        if (mStiffR3.rows() > 0) {
            mStiffR3 +=
            Eigen::Map<const eigen::RMatX3>(&buffer(row), mStiffR3.rows(), 3);
            row += mStiffR3.size();
        //} else if (mStiff3.rows() > 0) {
            // mStiff3 +=
            // Eigen::Map<const eigen::CMatX3>(&buffer(row), mStiff3.rows(), 3);
            // row += mStiff3.size();
        } else {
            throw std::runtime_error("SolidWindowSum::extractComm || Buffer allocation failed.");
        }
    }
    
    void overlapAndAddStiff() {
        if (mWindows.size() > 1) {
            mStiffR3.setZero();
            for (int m = 0; m < mWindows.size(); m++) {
                if (mAlignment[m] < 0) {
                    for (int idim = 0; idim < 3; idim++) {
                        mStiffR3.col(idim) += interpolateWinToWhole(mPhi, 
                                              mWindows[m]->getStiffForWindowSum(idim), 
                                              mWindows[m]->getPhiForWindowSum());
                    }
                } else {
                    for (int idim = 0; idim < 3; idim++) {
                        mStiffR3.block(mAlignment[m], idim, mWindows[m]->getNr(), 1) += mWindows[m]->getStiffForWindowSum(idim);
                    }
                }
            }
        } else if (mStiffR3.rows() > 0) {
            mStiffR3 = mWindows[0]->getStiffForCommR();
        } else if (mStiff3.rows() > 0) {
            mStiff3 = mWindows[0]->getStiffForCommC();
        }
    }
    
    void scatterStiffToWindows() {
        if (mWindows.size() > 1) {
            for (int m = 0; m < mWindows.size(); m++) {
                if (mAlignment[m] < 0) {
                    for (int idim = 0; idim < 3; idim++) {
                        mWindows[m]->collectStiffFromWindowSum(interpolateWholeToWin(mWindows[m]->getPhiForWindowSum(), 
                                              mStiffR3.col(idim), mPhi), idim);
                    }
                } else {
                    for (int idim = 0; idim < 3; idim++) {
                        mWindows[m]->collectStiffFromWindowSum(mStiffR3.block(mAlignment[m], idim, mWindows[m]->getNr(), 1), idim);
                    }
                }
            }
        } else if (mStiffR3.rows() > 0) {
            mWindows[0]->collectStiffFromMessaging(mStiffR3);
        } else if (mStiff3.size() > 0) {
            mWindows[0]->collectStiffFromMessaging(mStiff3);
        }
    }
    
    void addSolidWindow(const std::shared_ptr<SolidPointWindow> &win, const int aligned) {
        mWindows.push_back(win);
        mAlignment.push_back(aligned);
    }
    
private:
    eigen::RMatX3 mStiffR3 = eigen::RMatX3::Zero(0, 3);
    eigen::CMatX3 mStiff3 = eigen::CMatX3::Zero(0, 3);
    std::vector<std::shared_ptr<SolidPointWindow>> mWindows;
    std::vector<int> mAlignment;
    mutable int mSetCount = 0;
};

class FluidWindowSum: public WindowSum {
public:
    FluidWindowSum(const eigen::RColX phi, bool onlyOneWindow):
    WindowSum(phi) {
        if (!onlyOneWindow) {
            mStiffR1 = eigen::RColX::Zero(phi.rows(), 1);
        }
    }
    
    void setInFourier() {
        for (auto &win: mWindows) {
            win->setInFourier(!needsFT());
        }
    }
    
    bool onlyOneWindow() const {return (mWindows.size() == 1);};
    
    /////////////////////////// mpi ///////////////////////////
    bool needsFT() const {
        return (mWindows.size() > 1 || mWindows[0]->is3D());
    }
    
    void allocateComm() {
        //if (needsFT()) {
            mStiffR1 = eigen::RColX::Zero(mPhi.rows(), 1);
            for (auto &win: mWindows) {
                win->setInFourier(false);
            }
        //} else {
        //    mStiff1 = eigen::CColX::Zero(mWindows[0]->getNu_1(), 1);
        //}
    }
    
    // size for mpi communication
    int sizeComm() const {
        if (mStiffR1.rows() > 0) {
            return (int)mStiffR1.size();
        } else {
            return (int)mStiff1.size();
        }
    }
    
    // feed to mpi buffer
    void feedComm(eigen::RColX &buffer, int &row) const {
        if (mStiffR1.rows() > 0) {
            buffer.block(row, 0, mStiffR1.rows(), 1) = mStiffR1;
            row += mStiffR1.rows();
        //} else if (mStiff1.rows() > 0) {
            //buffer.block(row, 0, mStiff1.rows(), 1) = mStiff1;
            //row += mStiff1.rows();
        } else {
            throw std::runtime_error("FluidWindowSum::feedComm || Buffer allocation failed.");
        }
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::RColX &buffer, int &row) {
        if (mStiffR1.rows() > 0) {
            mStiffR1 += buffer.block(row, 0, mStiffR1.rows(), 1);
            row += mStiffR1.rows();
        //} else if (mStiff1.rows() > 0) {
            //mStiff1 += buffer.block(row, 0, mStiff1.rows(), 1);
            //row += mStiff1.rows();
        } else {
            throw std::runtime_error("FluidWindowSum::extractComm || Buffer allocation failed.");
        }
    }
    
    void overlapAndAddStiff() {
        if (mWindows.size() > 1) {
            mStiffR1.setZero();
            for (int m = 0; m < mWindows.size(); m++) {
                if (mAlignment[m] < 0) {
                    mStiffR1 += interpolateWinToWhole(mPhi, 
                                              mWindows[m]->getStiffForWindowSum(), 
                                              mWindows[m]->getPhiForWindowSum());
                } else {
                    mStiffR1.block(mAlignment[m], 0, mWindows[m]->getNr(), 1) += mWindows[m]->getStiffForWindowSum();
                }
            }
        } else if (mStiffR1.size() > 0) {
            mStiffR1 = mWindows[0]->getStiffForCommR();
        } else if (mStiff1.size() > 0) {
            mStiff1 = mWindows[0]->getStiffForCommC();
        }
    }
        
    void scatterStiffToWindows() {
        if (mWindows.size() > 1) {
            for (int m = 0; m < mWindows.size(); m++) {
                if (mAlignment[m] < 0) {
                    mWindows[m]->collectStiffFromWindowSum(interpolateWholeToWin(mWindows[m]->getPhiForWindowSum(), 
                                              mStiffR1, mPhi));
                } else {
                    mWindows[m]->collectStiffFromWindowSum(mStiffR1.block(mAlignment[m], 0, mWindows[m]->getNr(), 1));
                }
            }
        } else if (mStiffR1.rows() > 0) {
            mWindows[0]->collectStiffFromMessaging(mStiffR1);
        } else if (mStiff1.size() > 0) {
            mWindows[0]->collectStiffFromMessaging(mStiff1);
        }
    }
    
    void addFluidWindow(const std::shared_ptr<FluidPointWindow> &win, const int aligned) {
        mWindows.push_back(win);
        mAlignment.push_back(aligned);
    }
    
private:
    eigen::RColX mStiffR1 = eigen::RColX::Zero(0, 1);
    eigen::CColX mStiff1= eigen::CColX::Zero(0, 1);
    std::vector<std::shared_ptr<FluidPointWindow>> mWindows;
    std::vector<int> mAlignment;
};

#endif
