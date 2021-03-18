#ifndef LocalizedNr_hpp
#define LocalizedNr_hpp

#include "numerical.hpp"

struct LocalizedNr {
    LocalizedNr(const std::vector<std::pair<double,double>> phiNrNew):
    mPhiNrWin(phiNrNew) {
      // nothing
    }
  
    double getNrUpscaled(const int i) const {
        double dphi = mPhiNrWin[(i + 1) % M()].first - mPhiNrWin[i].first;
        setBounds2Pi(dphi);
        return (2 * numerical::dPi) / dphi * mPhiNrWin[i].second;
    }
    
    void setNrUpscaled(const int i, const double nr_new) {
        double dphi = mPhiNrWin[(i + 1) % M()].first - mPhiNrWin[i].first;
        setBounds2Pi(dphi);
        mPhiNrWin[i].second = dphi / (2 * numerical::dPi) * nr_new;
    }
    
    double getScaledMaxNrInRange(const double phi1, const double phi2, const int tol) const {
        int i1 = -1;
        int i2 = -1;
        for (int m = 0; m < M(); m++) {
            if (i1 < 0) {
                if (mPhiNrWin[m].first + tol > phi1) {
                    i1 = (m == 0) ? M() - 1 : m - 1;
                    if (i2 > 0) break;
                }
            }
            if (i2 < 0) {
                if (mPhiNrWin[m].first - tol > phi2) {
                    i2 = m;
                    if (i1 > 0) break;
                }
            }
        }

        double max_nr = 0;
        double dphi;
        if (i1 <= i2) {
            for (int i = i1; i <= i2; i++) {
                max_nr = std::max(max_nr, getNrUpscaled(i));
            }
            dphi = phi2 - phi1;
        } else { // window crosses x-axis
            for (int i = i1; i < M(); i++) {
                max_nr = std::max(max_nr, getNrUpscaled(i));
            }
            for (int i = 0; i <= i2; i++) {
                max_nr = std::max(max_nr, getNrUpscaled(i));
            }
            dphi = 2 * numerical::dPi + phi2 - phi1;
        }
        return dphi / (2 * numerical::dPi) * max_nr;
    }
    
    LocalizedNr merged(int i) const {
        i = i % M();
        int i_next = (i + 1) % M();
        LocalizedNr merged_wins(mPhiNrWin);
        merged_wins.mPhiNrWin.erase(merged_wins.mPhiNrWin.begin() + i);
        double maxNr = std::max(getNrUpscaled(i_next),getNrUpscaled(i));
        merged_wins.setNrUpscaled(i,maxNr);
        return merged_wins;
    }
    
    void splitWindow(const int i, const int n) {
        double dphi = mPhiNrWin[(i + 1) % M()].first - mPhiNrWin[i].first;
        setBounds2Pi(dphi);
        
        double nr_new = mPhiNrWin[i].second / n;
        
        mPhiNrWin.insert(mPhiNrWin.begin() + n, (n - 1), mPhiNrWin[i]);
        for (int j = 0; j < n; j++) {
            mPhiNrWin[i + j].first = mPhiNrWin[i].first + j * dphi;
            mPhiNrWin[i + j].second = nr_new;
        }
    }
    
    int M() const {
        return mPhiNrWin.size();
    }
    
    double getSize(int i) {
        int i_next = (i + 1) % M();
        double dphi = mPhiNrWin[i_next].first - mPhiNrWin[i].first;
        setBounds2Pi(dphi);
        return dphi;
    }
    
    int nextSmallWindow(const int minNr) const {
        for (int i = 0; i < M(); i++) {
            if (mPhiNrWin[i].second < minNr) return i;
        } 
        return -1;
    }
    
    double sumNr() const {
        double sumNr = 0;
        for (auto it = mPhiNrWin.begin(); it != mPhiNrWin.end(); it++) {
            sumNr += it->second;
        }
        return sumNr;
    }
    
    
    std::vector<std::pair<double,double>> mPhiNrWin;
    
private:
    void setBounds2Pi(double &phi) const {
        while (phi < 0) {
            phi += 2 * numerical::dPi;
        }
        while (phi >= 2 * numerical::dPi) {
            phi -= 2 * numerical::dPi;
        }
    }
};

#endif /* LocalizedNr_hpp */
