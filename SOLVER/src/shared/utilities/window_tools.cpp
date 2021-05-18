#include "window_tools.hpp"

namespace window_tools {
    eigen::DColX getWindowSlope(const eigen::DColX &phi) {
        return 0.5 * (1. - (numerical::dPi * phi).array().cos());
    }
    
    double getWindowSlope(const double &phi) {
        return 0.5 * (1. - std::cos(numerical::dPi * phi));
    }
    
    eigen::DRow4 unwrapPhi(eigen::DRow4 phi, double tol) { // more careful treatment of shape function
        if ((phi.array() - phi(0)).abs().maxCoeff() < tol) {
            phi(1) = phi(0);
            phi.tail(2).array() = phi(0) + 2 * numerical::dPi;
            return phi;
        }
        
        for (int i = 1; i <= 3; i++) {
            if (phi(i) < phi(i - 1)) {
                if (phi(i - 1) < phi(i) + tol) {
                    phi(i) = phi(i - 1);
                } else {
                    phi.tail(4 - i).array() += 2 * numerical::dPi;
                    if (std::abs(phi(3) - phi(2)) < tol) phi(2) = phi(3);
                }
            }
        }
        return phi;
    }
    
    double setBounds2Pi(double phi) {
        while (phi < 0) {
            phi += 2 * numerical::dPi;
        }
        while (phi >= 2 * numerical::dPi) {
            phi -= 2 * numerical::dPi;
        }
        return phi;
    }
};
