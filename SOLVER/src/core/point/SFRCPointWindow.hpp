#ifndef SFRCPointWindow_hpp
#define SFRCPointWindow_hpp

#include "PointWindow.hpp"
#include "Scanning1D.hpp"
#include "InterfaceSF.hpp"
#include "point_time.hpp"
#include <iostream>
class TimeScheme;

struct FluidStore {
    // pressure source (to be added to accel)
    eigen::CColX mPressureSource = eigen::CColX(0, 1);
    // acceleration storage for pressure output
    // NOTE: duplicated in Newmark
    eigen::CColX mPressureStore = eigen::CColX(0, 1);
    // stiffness storage for delta output
    eigen::CColX mDeltaStore = eigen::CColX(0, 1);
};

struct InterfaceRC {
    inline static void setTo(numerical::Real &out, const numerical::Real &displ) {out = displ;};
    inline static void setTo(numerical::ComplexR &out, const numerical::Real &displ) {
        throw std::runtime_error("SFRCPointWindow::scatterDisplToElementWindow || " 
            "requesting complex displ to real-type point window.");
    };
    
    inline static void setTo(numerical::Real &displ, const numerical::ComplexR &input) {
        throw std::runtime_error("SFRCPointWindow::scatterDisplToElementWindow || " 
            "requesting real displ to complex-type point window.");
    };
    inline static void setTo(numerical::ComplexR &displ, const numerical::ComplexR &input) {displ = input;};
    
    inline static void maskImag(Eigen::Matrix<numerical::Real, Eigen::Dynamic, 3> &val) {};
    inline static void maskImag(Eigen::Matrix<numerical::Real, Eigen::Dynamic, 1> &val) {};
    inline static void maskImag(Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 3> &val) {val.row(0).imag().setZero();};
    inline static void maskImag(Eigen::Matrix<numerical::ComplexR, Eigen::Dynamic, 1> &val) {val.row(0).imag().setZero();};
};

template <int dims, typename DisplType>
class SFRCPointWindow: public PointWindow {
public:
    typedef typename Fields<dims, numerical::Real>::RMat RMat;
    typedef typename Fields<dims, numerical::ComplexR>::CMat CMat;
    typedef typename std::vector<std::array<eigen::CMatPP_RM, dims>> vec_arSF_CMatPP_RM;
    typedef typename Eigen::Matrix<numerical::Real, Eigen::Dynamic, dims * spectral::nPEM> RMatXNSF;
  
    // constructor
    using PointWindow::PointWindow;
    
    void checkCompatibility(const TimeScheme &timeScheme) {
        // fields
        point_time::createFields(*this, timeScheme);
        // mass
        mMass->checkCompatibility(mNr, false);
        InterfaceSF<dims>::fftSF.addNR(mNr);
    };
    
    bool isFluid() const {return (dims == 1);};
    bool storesFieldsInFourier() const {return !(std::is_floating_point<DisplType>::value);};
    
    /////////////////////////// measure ///////////////////////////
    // random displ
    void randomDispl() {
        mFields.mDispl.setRandom();
        InterfaceRC::maskImag(mFields.mDispl);
    }
    
    // random stiff
    void randomStiff() {
        mFields.mStiffR.setRandom();
        mFields.mStiff.setRandom();
        mFields.mStiff.row(0).imag().setZero();
    }
    
    // reset to zero
    void resetToZero() {
        mFields.mStiff.setZero();
        mFields.mStiffR.setZero();
        mFields.mDispl.setZero();
        mFields.mVeloc.setZero();
        mFields.mAccel.setZero();
    }
    
    void resetStiffToZero() {
        mFields.mStiff.setZero();
        mFields.mStiffR.setZero();
    }
    
    /////////////////////////// time loop ///////////////////////////
    // check stability
    bool stable() const {
        return mFields.mDispl.allFinite();
    }
    
    // stiff to accel
    void transformToPhysical() {
        if (!mInFourier) {
            InterfaceSF<dims>::fftSF.computeC2R(mFields.mStiff, mFields.mStiffR, mNr);
        }
    };
    
    void transformToFourier() {
        if (!mInFourier) InterfaceSF<dims>::fftSF.computeR2C(mFields.mStiffR, mFields.mStiff, mNr);
    };
    
    void maskNyquist() {
        if (mNr % 2 == 0) mFields.mStiff.bottomRows(1).imag().setZero();
    };
    
    void computeStiffToAccel() {
        // stiff to accel in-place
        if (mInFourier) {
            mMass->computeAccel(mFields.mStiff);
        } else {
            mMass->computeAccel(mFields.mStiffR);
        }
    }
    
    void applyAxialBC() {
        transformToFourier();
        applyAxialBC_SF();
        transformToPhysical();
    }
    
    virtual void storeDelta() {};
    virtual void applyPressureSource() {};
    virtual void applyAxialBC_SF() = 0;
    
    /////////////////////////// window sum  ///////////////////////////
    void scatterDisplToElementWindow(vec_arSF_CMatPP_RM &displ, int nu_1_element, int ipol, int jpol) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            for (int idim = 0; idim < dims; idim++) {
                InterfaceRC::setTo(displ[alpha][idim](ipol, jpol), mFields.mDispl(alpha, idim));
            }
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = mNu_1; alpha < nu_1_element; alpha++) {
            for (int idim = 0; idim < dims; idim++) {
                InterfaceRC::setTo(displ[alpha][idim](ipol, jpol), czero);
            }
        }
    };
    
    void scatterDisplToElementWindow(RMatXNSF &displ, int ipnt) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNr; alpha++) {
            for (int idim = 0; idim < dims; idim++) {
                InterfaceRC::setTo(displ(alpha, ipnt + idim * spectral::nPEM), mFields.mDispl(alpha, idim));
            }
        }
    }
    
    void gatherStiffFromElementWindow(const vec_arSF_CMatPP_RM &stiff, int ipol, int jpol) {
        if (!storesFieldsInFourier()) 
            throw std::runtime_error("SFRCPointWindow::gatherStiffFromElementWindow || " 
            "adding complex stiff to real-type point window.");
        
        // add lower orders only
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            for (int idim = 0; idim < dims; idim++) {
                mFields.mStiff(alpha, idim) -= stiff[alpha][idim](ipol, jpol);
            }
        }
    }
    
    void gatherStiffFromElementWindow(const RMatXNSF &stiff, int ipnt) {
        if (storesFieldsInFourier()) 
            throw std::runtime_error("SFRCPointWindow::gatherStiffFromElementWindow || "
            "adding real stiff to complex-type point window.");
        
        for (int idim = 0; idim < dims; idim++) {
            mFields.mStiffR.col(idim) -= stiff.block(0, idim * spectral::nPEM + ipnt, mNr, 1);
        }
    }
                    
    void collectStiffFromWindowSum(const RMat &stiff) {mFields.mStiffR = stiff;};
    
    void collectStiffFromWindowSum(const RMat &stiff, const eigen::IColX &indices) {
        for (int i = 0; i < indices.rows(); i++) {
            mFields.mStiffR.row(i) = stiff.row(indices(i));
        }
    };
    
    void collectStiffFromWindowSum(const CMat &stiff) {
        if (!storesFieldsInFourier()) 
            throw std::runtime_error("SFRCPointWindow::collectStiffFromWindowSum || " 
            "adding complex stiff to real-type point window.");
            
        mFields.mStiff = stiff;
    };
    
    void enableScanning() {
        if (!storesFieldsInFourier()) {
            throw std::runtime_error("SFRCPointWindow::enableScanning || Wavefield scanning not implemented for mutliple windows.");
        }
        if (!mScanning[0]) {
           for (int i = 0; i < dims; i++) mScanning[i] = std::make_unique<Scanning1D>();
           std::make_unique<Scanning1D>();
        }
    }
    
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) {
        if (mScanning[0]) {
            for (int i = 0; i < dims; i++) mScanning[i]->doScanning(relTolFourierH2, relTolH2, absTolH2,
                               maxNumPeaks, mFields.mDispl.col(0));
        }
    }
    
    int reportScanningNr() const {
        int maxNr = -1;
        if (mScanning[0]) {
            for (int i = 0; i < dims; i++) {
                maxNr = std::max(maxNr, mScanning[i]->reportScanningNr(mNr));
            }
        }
        return maxNr;
    }
    
    virtual void disableScanning() {};
    
    virtual eigen::RMatX3 &getSolidStiffForWindowSum() {
        throw std::runtime_error("PointWindow::getSolidStiffFromWindowSum || "
                                 "requesting solid stiff to fluid window.");
    }
    
    virtual eigen::RColX &getFluidStiffForWindowSum() {
        throw std::runtime_error("PointWindow::getFluidStiffFromWindowSum || "
                                 "requesting fluid stiff to solid window.");
    }
    
    virtual Fields<1, numerical::Real> &getFluidFieldsR() {
        throw std::runtime_error("PointWindow::getFluidFieldsR || requesting fluid fields from a solid window.");
    }
    
    virtual Fields<3, numerical::Real> &getSolidFieldsR() {
        throw std::runtime_error("PointWindow::getSolidFieldsR || requesting solid fields from a fluid window.");
    }
    
    virtual const Fields<1, numerical::Real> &getFluidFieldsR() const {
        throw std::runtime_error("PointWindow::getFluidFieldsR || requesting fluid fields from a solid window.");
    }
    
    virtual const Fields<3, numerical::Real> &getSolidFieldsR() const {
        throw std::runtime_error("PointWindow::getSolidFieldsR || requesting solid fields from a fluid window.");
    }
    
    virtual Fields<1, numerical::ComplexR> &getFluidFieldsC() {
        throw std::runtime_error("PointWindow::getFluidFieldsC || requesting fluid fields from a solid window.");
    }
    
    virtual Fields<3, numerical::ComplexR> &getSolidFieldsC() {
        throw std::runtime_error("PointWindow::getSolidFieldsC || requesting solid fields from a fluid window.");
    }
    
    virtual const Fields<1, numerical::ComplexR> &getFluidFieldsC() const {
        throw std::runtime_error("PointWindow::getFluidFieldsC || requesting fluid fields from a solid window.");
    }
    
    virtual const Fields<3, numerical::ComplexR> &getSolidFieldsC() const {
        throw std::runtime_error("PointWindow::getSolidFieldsC || requesting solid fields from a fluid window.");
    }
    
protected:
    // fields on a solid point
    Fields<dims, DisplType> mFields;
    std::array<std::unique_ptr<Scanning1D>, dims> mScanning;
};

template <int dims, typename DisplType>
class SolidRCPointWindow: public SFRCPointWindow<dims, DisplType> {};

template <typename DisplType>
class SolidRCPointWindow<3, DisplType>: public SFRCPointWindow<3, DisplType> {
public:
    using SFRCPointWindow<3, DisplType>::SFRCPointWindow;
    
    void applyAxialBC_SF() {
        static const numerical::ComplexR czero = 0.;
        static const numerical::ComplexR cJ = {0., 1.};
        static const numerical::Real half = .5;
        
        // alpha = 0
        this->mFields.mStiff(0, 0) = czero;
        this->mFields.mStiff(0, 1) = czero;
        
        // alpha > 0
        if (this->mNu_1 - 1 >= 1) {
            // alpha = 1
            this->mFields.mStiff(1, 0) = (this->mFields.mStiff(1, 0) - cJ * this->mFields.mStiff(1, 1)) * half;
            this->mFields.mStiff(1, 1) = cJ * this->mFields.mStiff(1, 0);
            this->mFields.mStiff(1, 2) = czero;
            
            // alpha >= 2
            this->mFields.mStiff.bottomRows(this->mNu_1 - 2).setZero();
        }
    }
    
    void addForceSource(const eigen::CMatXN3 &force, int nu_1_force, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(this->mNu_1, nu_1_force);
        this->mFields.mStiff.topRows(nu_1_min) +=
        force(Eigen::seqN(Eigen::fix<0>, nu_1_min),
              Eigen::seqN(ipnt, Eigen::fix<3>, Eigen::fix<spectral::nPEM>));
    }
    
    
    eigen::RMatX3 &getSolidStiffForWindowSum() {return this->mFields.mStiffR;};
    
    Fields<3, numerical::Real> &getSolidFieldsR() {
        throw std::runtime_error("SolidRCPointWindow::getSolidFieldsR || requesting real-type fields from complex-type window.");
    }
    
    const Fields<3, numerical::Real> &getSolidFieldsR() const {
        throw std::runtime_error("SolidRCPointWindow::getSolidFieldsR || requesting real-type fields from complex-type window.");
    }
    
    Fields<3, numerical::ComplexR> &getSolidFieldsC() {
        throw std::runtime_error("SolidRCPointWindow::getSolidFieldsC || requesting complex-type fields from real-type window.");
    }
    
    const Fields<3, numerical::ComplexR> &getSolidFieldsC() const {
        throw std::runtime_error("SolidRCPointWindow::getSolidFieldsC || requesting complex-type fields from real-type window.");
    }
};

template <int dims, typename DisplType>
class FluidRCPointWindow: public SFRCPointWindow<dims, DisplType> {};

template <typename DisplType>
class FluidRCPointWindow<1, DisplType>: public SFRCPointWindow<1, DisplType> {
public:
    using SFRCPointWindow<1, DisplType>::SFRCPointWindow;
  
    void applyAxialBC_SF() {
        this->mFields.mStiff.bottomRows(this->mNu_1 - 1).setZero();
        this->maskNyquist();
    }
    
    void storeDelta() {
        if (mFluidStore.mDeltaStore.rows() > 0) {
            if (this->mInFourier) {
                mFluidStore.mDeltaStore = this->mFields.mStiff;
            } else {
                fft::gFFT_1.computeR2C(this->mFields.mStiffR, mFluidStore.mDeltaStore, this->mNr);
            }
        }
    }
    
    void storePressure() {
        // apply pressure source
        if (mFluidStore.mPressureSource.rows() > 0) {
            this->transformToFourier();
            // add pressure to new acceleration
            this->mFields.mStiff += mFluidStore.mPressureSource;
            // zero pressure for the next time step
            mFluidStore.mPressureSource.setZero();
            this->transformToPhysical();
        }
        
        // store acceleration for pressure output
        if (mFluidStore.mPressureStore.rows() > 0) {
            mFluidStore.mPressureStore = this->mFields.mStiff;
        }
    }
    
    /////////////////////////// source ///////////////////////////
    // prepare pressure source
    void preparePressureSource() {
        mFluidStore.mPressureSource = eigen::CColX::Zero(this->mNu_1, 1);
    }
    
    // add pressure source
    void addPressureSource(const eigen::CMatXN &pressure,
                           int nu_1_pressure, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(this->mNu_1, nu_1_pressure);
        mFluidStore.mPressureSource.topRows(nu_1_min) +=
        pressure.block(0, ipnt, nu_1_min, 1);
    }
    
    /////////////////////////// wavefield output ///////////////////////////
    // prepare pressure output
    void preparePressureOutput() {
        // pressure is mAccel, which may not be allocated by the time scheme
        if (mFluidStore.mPressureStore.rows() == 0) {
            mFluidStore.mPressureStore = eigen::CColX::Zero(this->mNu_1, 1);
        }
    }
    
    // prepare delta output
    void prepareDeltaOutput() {
        // delta is mStiff, but we need to store it because mStiff is set
        // to zero after dividing by mass
        if (mFluidStore.mDeltaStore.rows() == 0) {
            mFluidStore.mDeltaStore = eigen::CColX::Zero(this->mNu_1, 1);
        }
    }
    
    // scatter pressure to element
    void scatterPressureToElementWindow(eigen::CMatXN &pressure,
                                  int nu_1_element, int ipnt) const {
        // copy lower orders
        pressure.block(0, ipnt, this->mNu_1, 1) = mFluidStore.mPressureStore;
        
        // mask higher orders
        pressure.block(this->mNu_1, ipnt, nu_1_element - this->mNu_1, 1).setZero();
    }
    
    // scatter delta to element
    void scatterDeltaToElementWindow(eigen::CMatXN &delta,
                               int nu_1_element, int ipnt) const {
        // copy lower orders
        delta.block(0, ipnt, this->mNu_1, 1) = mFluidStore.mDeltaStore;
        
        // mask higher orders
        delta.block(this->mNu_1, ipnt, nu_1_element - this->mNu_1, 1).setZero();
    }
    
    void disableScanning() {
        this->mScanning[0].reset();
        this->mScanning[0] = nullptr;
    }
    
    eigen::RColX &getFluidStiffForWindowSum() {return this->mFields.mStiffR;};
    
    Fields<1, numerical::Real> &getFluidFieldsR() {
        throw std::runtime_error("FluidRCPointWindow::getFluidFieldsR || requesting real-type fields from complex-type window.");
    }
    
    const Fields<1, numerical::Real> &getFluidFieldsR() const {
        throw std::runtime_error("FluidRCPointWindow::getFluidFieldsR || requesting real-type fields from complex-type window.");
    }
    
    Fields<1, numerical::ComplexR> &getFluidFieldsC() {
        throw std::runtime_error("FluidRCPointWindow::getFluidFieldsC || requesting complex-type fields from real-type window.");
    }
    
    const Fields<1, numerical::ComplexR> &getFluidFieldsC() const {
        throw std::runtime_error("FluidRCPointWindow::getFluidFieldsC || requesting complex-type fields from real-type window.");
    }
    
private:
    FluidStore mFluidStore;  
};

#endif
