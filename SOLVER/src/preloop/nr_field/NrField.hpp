#ifndef NrField_hpp
#define NrField_hpp

#include <memory>
#include <numeric>
#include "eigen_sem.hpp"

#include "LocalizedNr.hpp"
#include "LocalizedNrField.hpp"

class ExodusMesh;
class GeneralNrField;

class NrField {
public:
    NrField (const ExodusMesh &exodusMesh);
    
    std::vector<std::pair<double, double>> makeNodalNrAtPoint(const eigen::DCol2 &sz) const;
    
    std::pair<eigen::DMatX2, eigen::IMatX4> makeElementalNrWindows(
            const std::vector<std::vector<std::pair<double, double>>> NodalNr, const double s) const;
            
    void finalizeNrWindows(std::vector<std::unique_ptr<std::tuple<eigen::DRow4, 
            eigen::IRowN, eigen::IRowN, bool>>> &quadWins, eigen::DRow2 phi_undivided, 
            eigen::IRowN nr_undivided, double phi2_prev, double phi1_next, double s) const;
            
    int getOverlapMinNr() const {return mOverlapMinNr;};
    
private:
    // horrible but necessary window maipulation
    static void combineNrWindows(std::vector<LocalizedNr> &wfv, const double angle_tol, const bool max_only);
    void removeSmallNrWindows(std::vector<LocalizedNr> &wfv, const double maxNrUpscaled) const;
    std::vector<std::vector<LocalizedNr>> getOptions(const std::vector<LocalizedNr> wfv, 
            int j_field, int i_small, const double maxNrUpscaled) const;
    
    // utility
    static int NrWindowVectorSum(const std::vector<LocalizedNr> &wfv) {
        return std::accumulate(wfv.begin(), wfv.end(), 0,
            [] (int sum, const LocalizedNr& wf) {
                return wf.M();
            }
        );
    }
    
    // Lucky numbers
    static bool isLuckyNumber(int n);
    static void nextLuckyNumber(eigen::IRowN &PointNr);
    static int nextLuckyNumber(int n);
    
    std::shared_ptr<const GeneralNrField> mGeneralNrField = nullptr;
    std::vector<std::unique_ptr<LocalizedNrField>> mLocalizedNrFields;
    
    bool mBoundByInplane, mUseLuckyNumbers, mAlign;
    
    double mDistTol, mAveGLLSpacing, mOverlapPhi;
    
    int mOverlapMinNr, mWindowMinNr, mWindowSize;
    
    static std::vector<int> sLuckyNumbers;
};

#endif /* NrField_hpp */
