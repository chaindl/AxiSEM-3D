#ifndef WindowSum_hpp
#define WindowSum_hpp

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include "numerical.hpp"
#include "eigen_point.hpp"
#include "SolidPointWindow.hpp"
#include "FluidPointWindow.hpp"

class WindowSum {
public:
    // constructor
    WindowSum(const eigen::RColX phi, bool onlyOneWindow):
    mPhi(phi), mOnlyOneWindow(onlyOneWindow) {
        // nothing
    }
    
    // destructor
    virtual ~WindowSum() = default;
    
    const eigen::RColX getPhi() const {return mPhi;};
    
    virtual int sizeComm() const = 0;
    
    virtual void setInFourier() = 0;
    
    // feed to mpi buffer
    virtual void feedComm(eigen::RColX &buffer, int &row) const = 0;
    
    // extract from mpi buffer
    virtual void extractComm(const eigen::RColX &buffer, int &row) const = 0;
    
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
    
    bool onlyOneWindow() {return mOnlyOneWindow;};
    
protected:
    const eigen::RColX mPhi;
    const bool mOnlyOneWindow; // this is redundant for now but will be needed for discontinuous fluids etc
};

class SolidWindowSum: public WindowSum {
public:
    SolidWindowSum(const eigen::RColX phi, bool onlyOneWindow):
    WindowSum(phi, onlyOneWindow) {
                
    }
    
    void setInFourier() {
        bool needFourier = (mWindows.size() > 1);
        for (auto &win: mWindows) {
            win->setInFourier(needFourier);
        }
    }
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    int sizeComm() const {
        return (int)sStiffR3.size();
    }
    
    // feed to mpi buffer
    void feedComm(eigen::RColX &buffer, int &row) const {
        buffer.block(row, 0, sStiffR3.size(), 1) =
        Eigen::Map<const eigen::RColX>(sStiffR3.data(),
                                       sStiffR3.size());
        row += sStiffR3.size();
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::RColX &buffer, int &row) const {
        sStiffR3 +=
        Eigen::Map<const eigen::RMatX3>(&buffer(row), sStiffR3.rows(), 3);
        row += sStiffR3.size();
    }
    
    void overlapAndAddStiff() {
        if (mWindows.size() == 1) return;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mAlignment[m] < 0) {
                for (int idim = 0; idim < 3; idim++) {
                    sStiffR3.col(idim) += interpolateWinToWhole(mPhi, 
                                          mWindows[m]->getStiffForWindowSum(idim), 
                                          mWindows[m]->getPhiForWindowSum());
                }
            } else {
                for (int idim = 0; idim < 3; idim++) {
                    sStiffR3.block(mAlignment[m], idim, mWindows[m]->getNr(), 1) += mWindows[m]->getStiffForWindowSum(idim);
                }
            }
        }
    }
    
    void scatterStiffToWindows() {
        if (mWindows.size() == 1) return;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mAlignment[m] < 0) {
                for (int idim = 0; idim < 3; idim++) {
                    mWindows[m]->collectStiffFromWindowSum(interpolateWholeToWin(mWindows[m]->getPhiForWindowSum(), 
                                          sStiffR3.col(idim), mPhi), idim);
                }
            } else {
                for (int idim = 0; idim < 3; idim++) {
                    mWindows[m]->collectStiffFromWindowSum(sStiffR3.block(mAlignment[m], idim, mWindows[m]->getNr(), 1), idim);
                }
            }
        }
        sStiffR3.setZero();
    }
    
    void addSolidWindow(const std::shared_ptr<SolidPointWindow> &win, const int aligned) {
        mWindows.push_back(win);
        mAlignment.push_back(aligned);
    }
    
private:
    inline static eigen::RMatX3 sStiffR3;
    std::vector<std::shared_ptr<SolidPointWindow>> mWindows;
    std::vector<int> mAlignment;
};

class FluidWindowSum: public WindowSum {
public:
    FluidWindowSum(const eigen::RColX phi, bool onlyOneWindow):
    WindowSum(phi, onlyOneWindow) {
                
    }
    
    void setInFourier() {
        bool needFourier = (mWindows.size() > 1);
        for (auto &win: mWindows) {
            win->setInFourier(needFourier);
        }
    }
    
    /////////////////////////// mpi ///////////////////////////
    // size for mpi communication
    int sizeComm() const {
        return (int)sStiffR1.size();
    }
    
    // feed to mpi buffer
    void feedComm(eigen::RColX &buffer, int &row) const {
        buffer.block(row, 0, sStiffR1.rows(), 1) = sStiffR1;
        row += sStiffR1.rows();
    }
    
    // extract from mpi buffer
    void extractComm(const eigen::RColX &buffer, int &row) const {
        sStiffR1 += buffer.block(row, 0, sStiffR1.rows(), 1);
        row += sStiffR1.rows();
    }
    
    void overlapAndAddStiff() {
        if (mWindows.size() == 1) return;
        for (int m = 0; m < mWindows.size(); m++) {
            if (mAlignment[m] < 0) {
                sStiffR1 += interpolateWinToWhole(mPhi, 
                                          mWindows[m]->getStiffForWindowSum(), 
                                          mWindows[m]->getPhiForWindowSum());
            } else {
                sStiffR1.block(mAlignment[m], 0, mWindows[m]->getNr(), 1) += mWindows[m]->getStiffForWindowSum();
            }
        }
    }
        
    void scatterStiffToWindows() {
        for (int m = 0; m < mWindows.size(); m++) {
            if (mAlignment[m] < 0) {
                mWindows[m]->collectStiffFromWindowSum(interpolateWholeToWin(mWindows[m]->getPhiForWindowSum(), 
                                          sStiffR1, mPhi));
            } else {
                mWindows[m]->collectStiffFromWindowSum(sStiffR1.block(mAlignment[m], 0, mWindows[m]->getNr(), 1));
            }
        }
        sStiffR1.setZero();
    }
    
    void addFluidWindow(const std::shared_ptr<FluidPointWindow> &win, const int aligned) {
        mWindows.push_back(win);
        mAlignment.push_back(aligned);
    }
    
private:
    inline static eigen::RColX sStiffR1;
    std::vector<std::shared_ptr<FluidPointWindow>> mWindows;
    std::vector<int> mAlignment;
};

#endif
