//
//  Sponge.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/29/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  sponge ABC

#ifndef Sponge_hpp
#define Sponge_hpp

#include "InterfaceSF.hpp"
#include "PointWindow.hpp"

///////////////////// Sponge class /////////////////////
template <int dims>
class Sponge {
public:
    // 1D constructor
    Sponge(const std::shared_ptr<PointWindow> &pw, double gamma):
    mPointWindow(pw), mGamma((eigen::RColX(1) << gamma).finished()) {
        // nothing
    }
    
    // 3D constructor
    Sponge(const std::shared_ptr<PointWindow> &pw, const eigen::DColX &gamma):
    mPointWindow(pw), mGamma(gamma.cast<numerical::Real>()) {
        // check size
        int nr = (int)gamma.rows();
        if (nr != pw->getNr()) {
            throw std::runtime_error("Sponge::Sponge || Incompatible sizes.");
        }
    }
    
    // get point
    const std::shared_ptr<const PointWindow> getPointWindow() const {
        return mPointWindow;
    }
    
    // apply ABC
    // must be called after point->computeStiffToAccel(),
    // so here "stiff" has been converted to acceleration
    virtual void apply() const = 0;
    
protected:
    // point
    const std::shared_ptr<PointWindow> mPointWindow;
    
    // gamma
    const eigen::RColX mGamma;
};

template <int dims>
class Sponge_R: public Sponge<dims> {
public:
    // 1D constructor
    Sponge_R(const std::shared_ptr<PointWindow> &pw, double gamma):
    Sponge<dims>(pw, gamma) {
        // nothing
    }
    
    // 3D constructor
    Sponge_R(const std::shared_ptr<PointWindow> &pw, const eigen::DColX &gamma):
    Sponge<dims>(pw, gamma) {
        // nothing
    }
    
    void apply() const {
        static const numerical::Real two = 2.;
      
        auto &stiff = InterfaceSF<dims>::getFieldsR(Sponge<dims>::mPointWindow).mStiffR;
        const auto &veloc = InterfaceSF<dims>::getFieldsR(Sponge<dims>::mPointWindow).mVeloc;
        const auto &displ = InterfaceSF<dims>::getFieldsR(Sponge<dims>::mPointWindow).mDispl;
      
        // update acceleration
        if (Sponge<dims>::mGamma.rows() == 1) {
            // 1D
            numerical::Real gamma = Sponge<dims>::mGamma(0);
            stiff -= (two * gamma) * veloc + (gamma * gamma) * displ;
        } else {
            // 3D
            stiff -= ((two * Sponge<dims>::mGamma).asDiagonal() * veloc +
                     Sponge<dims>::mGamma.array().square().matrix().asDiagonal() * displ);
        }    
    }
};

template <int dims>
class Sponge_C: public Sponge<dims> {
public:
    // 1D constructor
    Sponge_C(const std::shared_ptr<PointWindow> &pw, double gamma):
    Sponge<dims>(pw, gamma) {
        // nothing
    }
    
    // 3D constructor
    Sponge_C(const std::shared_ptr<PointWindow> &pw, const eigen::DColX &gamma):
    Sponge<dims>(pw, gamma) {
        
        int nr = (int)Sponge<dims>::mGamma.rows();  
        InterfaceSF<dims>::fftSF.addNR(nr);
        
        // workspace
        if (sStiff.rows() < nr) {
            sStiff.resize(nr, sStiff.cols());
            sDispl.resize(nr, sDispl.cols());
            sVeloc.resize(nr, sVeloc.cols());
        }
    }
    
    void apply() const {
        static const numerical::Real two = 2.;
        
        // get fields from point
        auto &stiff = InterfaceSF<dims>::getFieldsC(Sponge<dims>::mPointWindow).mStiff;
        const auto &veloc = InterfaceSF<dims>::getFieldsC(Sponge<dims>::mPointWindow).mVeloc;
        const auto &displ = InterfaceSF<dims>::getFieldsC(Sponge<dims>::mPointWindow).mDispl;
        
        // update acceleration
        if (Sponge<dims>::mGamma.rows() == 1) {
            // 1D
            numerical::Real gamma = Sponge<dims>::mGamma(0);
            stiff -= (two * gamma) * veloc + (gamma * gamma) * displ;
        } else {
            // 3D
            int nr = (int)Sponge<dims>::mGamma.rows();
            InterfaceSF<dims>::fftSF.computeC2R(stiff, sStiff, nr);
            InterfaceSF<dims>::fftSF.computeC2R(displ, sDispl, nr);
            InterfaceSF<dims>::fftSF.computeC2R(veloc, sVeloc, nr);
            sStiff -= ((two * Sponge<dims>::mGamma).asDiagonal() * sVeloc +
                       Sponge<dims>::mGamma.array().square().matrix().asDiagonal() * sDispl);
            InterfaceSF<dims>::fftSF.computeR2C(sStiff, stiff, nr);
        }
    }
    
private:
    //////////////////// static ////////////////////
    // 3D workspace
    inline static typename
    InterfaceSF<dims>::RMatSF sStiff, sVeloc, sDispl;
};

#endif /* Sponge_hpp */
