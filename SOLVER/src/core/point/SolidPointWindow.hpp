//
//  SolidPoint.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/31/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  solid GLL point

#ifndef SolidPointWindow_hpp
#define SolidPointWindow_hpp

#include "PointWindow.hpp"
#include "eigen_element.hpp"
#include "Scanning1D.hpp"

class TimeScheme;

class SolidPointWindow: public PointWindow {
public:
    // constructor
    SolidPointWindow(const eigen::RMatX2 &windowSumPhi,
               std::unique_ptr<const Mass> &mass,
               const TimeScheme &timeScheme,
               const std::shared_ptr<Point> point);
    
public:
    bool isFluid() const {return false;};
  
    /////////////////////////// measure ///////////////////////////
    // random displ
    void randomDispl() {
        mFields.mDispl.setRandom();
        mFields.mDispl.row(0).imag().setZero();
    }
    
    // random stiff
    void randomStiff() {
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
    
    /////////////////////////// time loop ///////////////////////////
    // check stability
    bool stable() const {
        return mFields.mDispl.allFinite();
    }
    
    // stiff to accel
    void transformToPhysical();
    void computeStiffToAccel();
    void transformToFourier();
    void maskNyquist();
    void applyAxialBC();
    
    /////////////////////////// window sum  ///////////////////////////

    eigen::RColX getStiffForWindowSum(const int i) {return mFields.mStiffR.col(i);};
    void collectStiffFromWindowSum(const eigen::RColX &stiff, const int i) {
        mFields.mStiffR.col(i) = stiff.cwiseProduct(mWindowSumFrac);
    }
    
    eigen::RMatX3 getStiffForCommR() {return mFields.mStiffR;};
    eigen::CMatX3 getStiffForCommC() {return mFields.mStiff;};
    void collectStiffFromMessaging(const eigen::RMatX3 &stiff) {mFields.mStiffR = stiff;};
    void collectStiffFromMessaging(const eigen::CMatX3 &stiff) {mFields.mStiff = stiff;};
    
    /////////////////////////// element ///////////////////////////
    // scatter displ to element
    void scatterDisplToElementWindow(eigen::vec_ar3_CMatPP_RM &displ,
                               int nu_1_element, int ipol, int jpol) const {
        // copy lower orders
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            displ[alpha][0](ipol, jpol) = mFields.mDispl(alpha, 0);
            displ[alpha][1](ipol, jpol) = mFields.mDispl(alpha, 1);
            displ[alpha][2](ipol, jpol) = mFields.mDispl(alpha, 2);
        }
        
        // mask higher orders
        static const numerical::ComplexR czero = 0.;
        for (int alpha = mNu_1; alpha < nu_1_element; alpha++) {
            displ[alpha][0](ipol, jpol) = czero;
            displ[alpha][1](ipol, jpol) = czero;
            displ[alpha][2](ipol, jpol) = czero;
        }
    }
    
    // gather stiff from element
    void gatherStiffFromElementWindow(const eigen::vec_ar3_CMatPP_RM &stiff,
                                int ipol, int jpol) {
        // add lower orders only
        for (int alpha = 0; alpha < mNu_1; alpha++) {
            mFields.mStiff(alpha, 0) -= stiff[alpha][0](ipol, jpol);
            mFields.mStiff(alpha, 1) -= stiff[alpha][1](ipol, jpol);
            mFields.mStiff(alpha, 2) -= stiff[alpha][2](ipol, jpol);
        }
        if (!mD && mFields.mStiff.real().norm() > 0) mD = true;
        //if (mD) std::cout << "stiff gathered" << std::endl << mFields.mStiff << std::endl << std::endl;
    }
    
    
    /////////////////////////// source ///////////////////////////
    // add force source (external)
    void addForceSource(const eigen::CMatXN3 &force,
                        int nu_1_force, int ipnt) {
        // add minimum orders only
        int nu_1_min = std::min(mNu_1, nu_1_force);
        mFields.mStiff.topRows(nu_1_min) +=
        force(Eigen::seqN(Eigen::fix<0>, nu_1_min),
              Eigen::seqN(ipnt, Eigen::fix<3>, Eigen::fix<spectral::nPEM>));
    }
    
    
    /////////////////////////// fields ///////////////////////////

    SolidPointWindow* getSolidPointWindow() {return this;};
    
    // get
    const Fields<3> &getSolidFields() const {
        return mFields;
    }
    
    // set
    Fields<3> &getSolidFields() {
        return mFields;
    }
    
private:
    // fields on a solid point
    Fields<3> mFields;
    bool mD = false;
    
    /////////////////////////// wavefield scanning ///////////////////////////
public:
    // enable scanning
    void enableScanning() {
        if (!mScanningS) {
            mScanningS = std::make_unique<Scanning1D>();
            mScanningP = std::make_unique<Scanning1D>();
            mScanningZ = std::make_unique<Scanning1D>();
        }
    }
    
    // do scanning
    void doScanning(numerical::Real relTolFourierH2, numerical::Real relTolH2,
                    numerical::Real absTolH2, int maxNumPeaks) {
        if (mScanningS) {
            mScanningS->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(0));
            mScanningP->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(1));
            mScanningZ->doScanning(relTolFourierH2, relTolH2, absTolH2,
                                   maxNumPeaks, mFields.mDispl.col(2));
        }
    }
    
    // report scanning Nr
    int reportScanningNr() const {
        if (mScanningS) {
            int nrS = mScanningS->reportScanningNr(mNr);
            int nrP = mScanningP->reportScanningNr(mNr);
            int nrZ = mScanningZ->reportScanningNr(mNr);
            return std::max(std::max(nrS, nrP), nrZ);
        } else {
            return -1;
        }
    }
    
private:
    // scanning
    std::unique_ptr<Scanning1D> mScanningS = nullptr;
    std::unique_ptr<Scanning1D> mScanningP = nullptr;
    std::unique_ptr<Scanning1D> mScanningZ = nullptr;
};

#endif /* SolidPoint_hpp */
