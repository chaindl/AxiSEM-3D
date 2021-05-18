#ifndef window_tools_hpp
#define window_tools_hpp

#include "numerical.hpp"
#include "eigen_sem.hpp"

namespace window_tools {
    eigen::DColX getWindowSlope(const eigen::DColX &phi);
    double getWindowSlope(const double &phi);
    
    template <typename EigenVec>
    void wrapPhi(EigenVec &phi) {
        int i_wrap = -1;
        for (int i = 0; i < phi.size(); i++) {
            if (phi(i) >= 2 * numerical::dPi) {
                i_wrap = i;
                break;
            }
        }
        if (i_wrap > 0) phi.tail(phi.size() - i_wrap).array() -= 2 * numerical::dPi;
    }
    
    template <typename EigenVec>
    void unwrapPhi(EigenVec &phi) {
        int i_wrap = -1;
        for (int i = 1; i < phi.size(); i++) {
            if (phi(i) <= phi(0)) {
                i_wrap = i;
                break;
            }
        }
        if (i_wrap > 0) phi.tail(phi.size() - i_wrap).array() += 2 * numerical::dPi;
    }
    
    eigen::DRow4 unwrapPhi(eigen::DRow4 phi, double tol);
    double setBounds2Pi(double phi);
    
    inline int nextWin(int m, int M) {return (m + 1) %  M;}
    inline int prevWin(int m, int M) {return (m + M - 1) % M;}
}

#endif
