#include "NrField.hpp"
#include "inparam.hpp"
#include "LocalizedNr.hpp"
#include "LocalizedNrField.hpp"
#include "GeneralNrField.hpp"
#include "window_tools.hpp"
#include "ExodusMesh.hpp"

NrField::NrField(const ExodusMesh &exodusMesh) {
    mDistTol = exodusMesh.getGlobalVariable("dist_tolerance");
  
    mGeneralNrField = GeneralNrField::buildInparam(mDistTol);
    LocalizedNrField::buildInparam(mLocalizedNrFields, mGeneralNrField, mDistTol);

    mBoundByInplane = inparam::gInparamNr.get<bool>("bound_Nr_by_inplane");
    if (mBoundByInplane) mAveGLLSpacing = exodusMesh.computeAveGLLSpacing();
    mUseLuckyNumbers = inparam::gInparamAdvanced.get<bool>("develop:fftw_lucky_numbers");
    
    mOverlapMinNr = inparam::gInparamNr.get<int>("Nr_windows:min_points_overlap"); // including 0 point but not including 1 point, min = 2
    mOverlapPhi = inparam::gInparamNr.get<double>("Nr_windows:phi_overlap") * numerical::dDegree;
    mWindowSize = inparam::gInparamNr.get<int>("Nr_windows:azimuthal_window_size");
    mAlign = inparam::gInparamNr.get<bool>("Nr_windows:align_azimuthal_sampling_points");

    mWindowMinNr = 2 * mOverlapMinNr + 1;
}

std::vector<std::pair<double, double>> NrField::makeNodalNrAtPoint(const eigen::DCol2 &sz) const {
    // create a single window if using global nr only
    if (mLocalizedNrFields.size() == 0) {
        int nr_global = mGeneralNrField->getNrAtPoint(sz);
        if (mBoundByInplane) {
            int maxNrUpscaled = std::max((int)round(2 * numerical::dPi * sz(0) / mAveGLLSpacing), 1);
            if (sz(0) < mDistTol) maxNrUpscaled = 5;
            nr_global = std::min(nr_global, maxNrUpscaled);
        }
        std::vector<std::pair<double, double>> NrWindows(1, std::make_pair(0., nr_global));
        return NrWindows;
    }
    
    // combine window fields from all input types
    std::vector<LocalizedNr> NrWindows;
    for (auto &field: mLocalizedNrFields) {
        NrWindows.push_back(field->getWindowsAtPoint(sz));
    }
    combineNrWindows(NrWindows, mDistTol / sz(0), true);
    
    // limit by inplane resolution
    double maxNrUpscaled = -1;
    if (mBoundByInplane) {
        maxNrUpscaled = 2 * numerical::dPi * sz(0) / mAveGLLSpacing;
        if (sz(0) < mDistTol) maxNrUpscaled = 5;
        for (int i = 0; i < NrWindows[0].M(); i++) {
            if (NrWindows[0].getNrUpscaled(i) > maxNrUpscaled) {
                NrWindows[0].setNrUpscaled(i, maxNrUpscaled);
            }
        }
    }
    
    // handle windows with less than min nr
    // compare 3 options and choose the one with least total nr increase:
    // 1. increase window nr to min nr (only if possible based on inplane resolution)
    // 2. merge with previous window and take max nr
    // 3. merge with next window and take max nr
    //
    // this method is a brute-force loop through all option trees
    // it could potentially lead to memory overflow 
    // but such a situation should not occur given practical limitations of
    // the window functions
    removeSmallNrWindows(NrWindows, maxNrUpscaled);
    return NrWindows[0].mPhiNrWin;
}

std::pair<eigen::DMatX2, eigen::IMatX4> NrField::makeElementalNrWindows(
    const std::vector<std::vector<std::pair<double, double>>> NodalNr, const double s) const {
      
    std::vector<LocalizedNr> ElementalNr;
    bool local_fields = false;
    for (int inode = 0; inode < 4; inode++) {
        LocalizedNr ews(NodalNr[inode]);
        ElementalNr.push_back(ews);
        if (ews.M() > 1) local_fields = true;
    }
    
    if (!local_fields) {
        std::pair<eigen::DMatX2, eigen::IMatX4> elementalNrWindows(eigen::DMatX2::Zero(1, 2), eigen::IMatX4::Zero(1, 4));
        elementalNrWindows.first(0, 0) = 0.;
        elementalNrWindows.first(0, 1) = 0.;
        for (int i = 0; i < 4; i++) {
            elementalNrWindows.second(0, i) = (int)round(ElementalNr[i].mPhiNrWin[0].second);
        }
        return elementalNrWindows;
    }

    double maxNrUpscaled = -1;
    if (mBoundByInplane) {
        maxNrUpscaled = 2 * numerical::dPi * s / mAveGLLSpacing;
    }
    double angle_tol = mDistTol / s;
    
    combineNrWindows(ElementalNr, angle_tol, false);
    removeSmallNrWindows(ElementalNr, maxNrUpscaled);
    
    // calculate overlap
    std::vector<double> half_ol(ElementalNr[0].M());
    if (ElementalNr[0].M() <= 1) {
        // nothing
    } else if (mOverlapPhi > 0) {
        std::fill(half_ol.begin(), half_ol.end(), mOverlapPhi / 2);
    } else {
        for (int m = 0; m < ElementalNr[0].M(); m++) {
           int m_next = (m + 1) % ElementalNr[0].M();
           
           double min_nr_upscaled = std::numeric_limits<double>::max();
           for (int inode = 0; inode < 4; inode++) {
                min_nr_upscaled = std::min(min_nr_upscaled, ElementalNr[inode].getNrUpscaled(m));
                min_nr_upscaled = std::min(min_nr_upscaled, ElementalNr[inode].getNrUpscaled(m_next));
            }
            half_ol[m] = 2 * numerical::dPi / min_nr_upscaled * mOverlapMinNr / 2;
        }
    }
    
    // window size is consistent between nodes now -> add overlap and store as matrix
    std::pair<eigen::DMatX2, eigen::IMatX4> elementalNrWindows(eigen::DMatX2::Zero(ElementalNr[0].M(), 2), eigen::IMatX4::Zero(ElementalNr[0].M(), 4));
    for (int m = 0; m < ElementalNr[0].M(); m++) {
        double phi1, phi2, fracSizeIncrease;
        
        if (ElementalNr[0].M() <= 1) {
            // no overlap for single window
            phi1 = 0;
            phi2 = 0;
            
            fracSizeIncrease = 1;
        } else {
            int m_next = (m + 1) % ElementalNr[0].M();
              
            phi1 = window_tools::setBounds2Pi(ElementalNr[0].mPhiNrWin[m].first - half_ol[m]);
            phi2 = window_tools::setBounds2Pi(ElementalNr[0].mPhiNrWin[m_next].first + half_ol[m]);

            fracSizeIncrease = (1. + 2 * half_ol[m] / ElementalNr[0].getSize(m));
        }
        
        eigen::IRow4 nr;
        for (int inode = 0; inode < 4; inode++) {
            ElementalNr[inode].setNrUpscaled(m, ElementalNr[inode].getNrUpscaled(m) * fracSizeIncrease);
            nr(inode) = (int)round(ElementalNr[inode].mPhiNrWin[m].second);
        }
        
        elementalNrWindows.first(m, 0) = phi1;
        elementalNrWindows.first(m, 1) = phi2;
        elementalNrWindows.second.row(m) = nr;
    }
    
    return elementalNrWindows;
}
    
void NrField::finalizeNrWindows(
    std::vector<std::unique_ptr<std::tuple<eigen::DRow4, eigen::IRowN, int, bool>>> &quadWins,
    eigen::DRow2 phi_undivided, eigen::IRowN nr_undivided,
    double phi2_prev, double phi1_next, double s) const {
     
    bool subdivision = false;
    int nwin_prev = quadWins.size();
    
    double angle_tol = mDistTol / s;
    double win_size = phi_undivided(1) - phi_undivided(0);
    if (win_size <= angle_tol) win_size += 2 * numerical::dPi;
    
    // attempt to divide into regular windows
    if (mWindowSize > 0) {
        // calculate overlap
        int min_nr = nr_undivided.minCoeff();
        
        double overlap;
        if (mOverlapPhi > 0) {
            overlap = mOverlapPhi;
        } else {
            overlap = win_size / min_nr * mOverlapMinNr;
        }
        double regWinPhi = overlap * mWindowSize; // inner size
        
        // recalculate exact window size for regular intervals
        double win_size_corrected = window_tools::setBounds2Pi(phi2_prev - phi1_next);
        if (win_size_corrected < angle_tol) win_size_corrected += 2 * numerical::dPi;
        win_size_corrected += overlap;
        int nwin = (int)floor((win_size_corrected + angle_tol) / (regWinPhi));
        regWinPhi = win_size / nwin;

        if (nwin > 1 && min_nr / nwin >= mWindowMinNr) {
            subdivision = true;
            
            quadWins.resize(quadWins.size() + nwin);
            
            eigen::IRowN nr_sub;
            if (mAlign) {
                eigen::IRowN nr_ol = (nr_undivided.array().cast<double>() * overlap / win_size).array().round().cast<int>();
                nr_sub = (nr_ol.array() * (mWindowSize + 1)) + 1;
            } else {
                nr_sub = (nr_undivided.array().cast<double>() * (regWinPhi + overlap) / win_size).array().round().cast<int>() + 1;
                if (mUseLuckyNumbers) {
                    nr_sub = nr_sub.unaryExpr([](int nr) {return nextLuckyNumber(nr);});
                }
            }
            
            if (win_size >= 2 * numerical::dPi - angle_tol) { // single window -> add overlap
                phi_undivided(0) = window_tools::setBounds2Pi(phi_undivided(0) - overlap / 2);
                phi_undivided(1) = window_tools::setBounds2Pi(phi_undivided(1) + overlap / 2);
                phi2_prev = phi_undivided(1);
                phi1_next = phi_undivided(0);  
            }
            
            
            
            // to easily calculate regular window positions we 
            // ignore potentially different overlaps from other windows
            double phi_start = window_tools::setBounds2Pi(phi2_prev - overlap);
            for (int m = 0; m < nwin; m++) {
                eigen::DRow4 phi;
                phi << window_tools::setBounds2Pi(phi_start + m * regWinPhi),
                       window_tools::setBounds2Pi(phi_start + m * regWinPhi + overlap),
                       window_tools::setBounds2Pi(phi_start + (m + 1) * regWinPhi),
                       window_tools::setBounds2Pi(phi_start + (m + 1) * regWinPhi + overlap);
                std::tuple<eigen::DRow4, eigen::IRowN, int, bool> win_tuple = std::make_tuple(phi, nr_sub, -1, true);
                quadWins[nwin_prev + m] = std::make_unique<std::tuple<eigen::DRow4, eigen::IRowN, int, bool>>();
                *quadWins[nwin_prev + m] = std::make_tuple(phi, nr_sub, -1, true);
            }
            
            // re-introduce potentially different overlaps from start and end
            if (win_size < 2 * numerical::dPi - angle_tol) { // not needed for single window
                double fracSizeDiff;
                // first window
                (std::get<0>(*quadWins[nwin_prev]))(0) = phi_undivided(0);
                fracSizeDiff = (regWinPhi + phi2_prev - phi_undivided(0)) / (regWinPhi + overlap);
                std::get<1>(*quadWins[nwin_prev]) = ((std::get<1>(*quadWins[nwin_prev])).array().cast<double>() * fracSizeDiff).array().round().cast<int>();
                // last window
                (std::get<0>(*quadWins[nwin_prev + nwin - 1]))(3) = phi_undivided(1);
                fracSizeDiff = (regWinPhi + phi_undivided(1) - phi1_next) / (regWinPhi + overlap);
                std::get<1>(*quadWins[nwin_prev + nwin - 1]) = ((std::get<1>(*quadWins[nwin_prev + nwin - 1])).array().cast<double>() * fracSizeDiff).array().round().cast<int>();
                // apply lucky numbers
                if (mUseLuckyNumbers) {
                    std::get<1>(*quadWins[nwin_prev]) = (std::get<1>(*quadWins[nwin_prev])).unaryExpr([](int nr) {return nextLuckyNumber(nr);});
                    std::get<1>(*quadWins[nwin_prev + nwin - 1]) = (std::get<1>(*quadWins[nwin_prev + nwin - 1])).unaryExpr([](int nr) {return nextLuckyNumber(nr);});
                }
            }
        }
    }
    
    if (!subdivision) {
        // no subdivision
        if (mUseLuckyNumbers) {
            nr_undivided = nr_undivided.unaryExpr([](int nr) {return nextLuckyNumber(nr);});
        }

        eigen::DRow4 phi ;
        phi << phi_undivided(0), phi2_prev, phi1_next, phi_undivided(1);
        bool hasOverlap = (win_size < 2 * numerical::dPi - angle_tol); // no overlap for single window 
        
        quadWins.push_back(std::make_unique<std::tuple<eigen::DRow4, eigen::IRowN, int, bool>>());
        *quadWins[nwin_prev] = std::make_tuple(phi, nr_undivided, -1, hasOverlap);
    }
}

void NrField::combineNrWindows(std::vector<LocalizedNr> &wfv, const double angle_tol, const bool max_only) {
    if (wfv.size() <= 1) return;
    
    // calculate combined size for preallocation
    int nr_raw = NrWindowVectorSum(wfv);
    
    // gather all phi
    std::vector<double> phi_combined;
    for (auto &it: wfv) {
        std::vector<double> phi(it.M());
        for (int m = 0; m < it.M(); m++) phi[m] = it.mPhiNrWin[m].second;
        phi_combined.insert(phi_combined.end(), std::begin(phi), std::end(phi));
    }
    
    // remove duplicated based on dist_tol
    std::sort(phi_combined.begin(), phi_combined.end());
    auto ip = std::unique(phi_combined.begin(), phi_combined.end(),
        [angle_tol] (const double& phi1, const double& phi2) {
            return phi1 >= phi2 - angle_tol;
        }
    );
    phi_combined.resize(distance(phi_combined.begin(), ip));
    
    // calculate nr of each input window field in the new windows
    eigen::DMatXX new_nr_local = eigen::DMatXX::Zero(phi_combined.size(), wfv.size());
    for (int i = 0; i < phi_combined.size(); i++) {
        int i_next = (i + 1) % phi_combined.size();
        for (int j = 0; j < wfv.size(); j++) {
            new_nr_local(i,j) = wfv[j].getScaledMaxNrInRange(phi_combined[i], phi_combined[i_next], angle_tol);
        }
    }
    
    // make new window field(s)
    std::vector<LocalizedNr> combinedWindowFields;
    
    if (max_only) { // for nodal (create a single window field)
        std::vector<std::pair<double, double>> new_PhiNr(phi_combined.size());
        for (int i = 0; i < phi_combined.size(); i++) {
            new_PhiNr[i] = std::make_pair(phi_combined[i], new_nr_local.row(i).maxCoeff());
        }
        LocalizedNr MaxNrCombinedWindowField(new_PhiNr);
        combinedWindowFields.push_back(MaxNrCombinedWindowField);
    } else { // for elemental (nr can differ between nodes but windows must be the same)
        for (int j = 0; j != wfv.size(); j++) {
            std::vector<std::pair<double, double>> new_PhiNr(phi_combined.size());
            for (int i = 0; i < phi_combined.size(); i++) {
                new_PhiNr[i] = std::make_pair(phi_combined[i], new_nr_local(i,j));
            }
            LocalizedNr CombinedWindowField(new_PhiNr);
            combinedWindowFields.push_back(CombinedWindowField);
        }
    }
    
    wfv = combinedWindowFields;
}

void NrField::removeSmallNrWindows(std::vector<LocalizedNr> &wfv, const double maxNrUpscaled) const {
    if (wfv[0].M() < 2) {
        for (auto &wf: wfv) {
            wf.mPhiNrWin[0].second = std::max(wf.mPhiNrWin[0].second, 1.);
        }
        return;
    }
    
    int i_small, j_node;
    for (int j = 0; j < wfv.size(); j++) {
        i_small = wfv[j].nextSmallWindow(mWindowMinNr);
        if (i_small >= 0) {
            j_node = j;
            break;
        }
    }
    if (i_small < 0) return;
    
    std::vector<std::vector<LocalizedNr>> options;
    options.push_back(wfv);

    int io = 0;
    while (i_small >= 0) {
      
        wfv = options[io];
        options.erase(options.begin() + io);
        std::vector<std::vector<LocalizedNr>> new_options = getOptions(wfv, j_node, i_small, maxNrUpscaled);
        options.insert(options.end(), new_options.begin(), new_options.end());
        
        bool found = false;
        for (int io_ = 0; io_ < options.size(); io_++) {
            for (int j = 0; j < wfv.size(); j++) {
                i_small = options[io_][j].nextSmallWindow(mWindowMinNr);
                if (i_small >= 0) {
                    io = io_;
                    j_node = j;
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
    }
    
    std::vector<LocalizedNr> best_option = *std::min_element(options.begin(), options.end(),
        [] (const std::vector<LocalizedNr>& o1, const std::vector<LocalizedNr>& o2) {
            return (NrWindowVectorSum(o1) < NrWindowVectorSum(o2));
        }
    );
    
    wfv = best_option;
}

std::vector<std::vector<LocalizedNr>> NrField::getOptions(const std::vector<LocalizedNr> wfv, int j_field, int i_small, const double maxNrUpscaled) const {
        std::vector<std::vector<LocalizedNr>> options;
        
        // option 1 = increase local nr to required minimum nr
        std::vector<LocalizedNr> option = wfv;
        option[j_field].mPhiNrWin[i_small].second = (double)mWindowMinNr;
        if (maxNrUpscaled < 0 || option[j_field].getNrUpscaled(i_small) < maxNrUpscaled) options.push_back(option);
        
        // option 2 = merge with previous window
        for (int j = 0; j < wfv.size(); j++) {
            option[j] = wfv[j].merged(i_small);
        }
        options.push_back(option);
        
        // option 3 = merge with next window
        if (wfv[0].M() > 2) {
            for (int j = 0; j < wfv.size(); j++) {
                option[j] = wfv[j].merged(i_small + 1);
            }
            options.push_back(option);
        }
        
        return options;
}

// is lucky number
bool NrField::isLuckyNumber(int n) {
    // check: max prime factor <= 13
    int num = n;
    for (int i = 2; i <= num; i++) {
        while (num % i == 0) {
            if (i > 13) {
                return false;
            }
            num /= i;
        }
    }
    
    // check: pow of 11 + pow of 13 <= 1
    num = n;
    int e = 0;
    while (num % 11 == 0) {
        num /= 11;
        e++;
    }
    int f = 0;
    while (num % 13 == 0) {
        num /= 13;
        f++;
    }
    if (e + f > 1) {
        return false;
    }
    
    // true
    return true;
}

void NrField::nextLuckyNumber(eigen::IRowN &PointNr) {
    PointNr = PointNr.unaryExpr([](int nr) {
        return nextLuckyNumber(nr);
    });
}

// next lucky number
int NrField::nextLuckyNumber(int n) {
    auto it = std::lower_bound(sLuckyNumbers.begin(),
                               sLuckyNumbers.end(), n);
    if (it == sLuckyNumbers.end()) {
        // no optimization for large numbers
        return n;
    }
    return *it;
}

// lucky numbers
std::vector<int> NrField::sLuckyNumbers = {
    1    , 2    , 3    , 4    , 5    , 6    , 7    , 8    , 9    , 10   ,
    11   , 12   , 13   , 14   , 15   , 16   , 18   , 20   , 21   , 22   ,
    24   , 25   , 26   , 27   , 28   , 30   , 32   , 33   , 35   , 36   ,
    39   , 40   , 42   , 44   , 45   , 48   , 49   , 50   , 52   , 54   ,
    55   , 56   , 60   , 63   , 64   , 65   , 66   , 70   , 72   , 75   ,
    77   , 78   , 80   , 81   , 84   , 88   , 90   , 91   , 96   , 98   ,
    99   , 100  , 104  , 105  , 108  , 110  , 112  , 117  , 120  , 125  ,
    126  , 128  , 130  , 132  , 135  , 140  , 144  , 147  , 150  , 154  ,
    156  , 160  , 162  , 165  , 168  , 175  , 176  , 180  , 182  , 189  ,
    192  , 195  , 196  , 198  , 200  , 208  , 210  , 216  , 220  , 224  ,
    225  , 231  , 234  , 240  , 243  , 245  , 250  , 252  , 256  , 260  ,
    264  , 270  , 273  , 275  , 280  , 288  , 294  , 297  , 300  , 308  ,
    312  , 315  , 320  , 324  , 325  , 330  , 336  , 343  , 350  , 351  ,
    352  , 360  , 364  , 375  , 378  , 384  , 385  , 390  , 392  , 396  ,
    400  , 405  , 416  , 420  , 432  , 440  , 441  , 448  , 450  , 455  ,
    462  , 468  , 480  , 486  , 490  , 495  , 500  , 504  , 512  , 520  ,
    525  , 528  , 539  , 540  , 546  , 550  , 560  , 567  , 576  , 585  ,
    588  , 594  , 600  , 616  , 624  , 625  , 630  , 637  , 640  , 648  ,
    650  , 660  , 672  , 675  , 686  , 693  , 700  , 702  , 704  , 720  ,
    728  , 729  , 735  , 750  , 756  , 768  , 770  , 780  , 784  , 792  ,
    800  , 810  , 819  , 825  , 832  , 840  , 864  , 875  , 880  , 882  ,
    891  , 896  , 900  , 910  , 924  , 936  , 945  , 960  , 972  , 975  ,
    980  , 990  , 1000 , 1008 , 1024 , 1029 , 1040 , 1050 , 1053 , 1056 ,
    1078 , 1080 , 1092 , 1100 , 1120 , 1125 , 1134 , 1152 , 1155 , 1170 ,
    1176 , 1188 , 1200 , 1215 , 1225 , 1232 , 1248 , 1250 , 1260 , 1274 ,
    1280 , 1296 , 1300 , 1320 , 1323 , 1344 , 1350 , 1365 , 1372 , 1375 ,
    1386 , 1400 , 1404 , 1408 , 1440 , 1456 , 1458 , 1470 , 1485 , 1500 ,
    1512 , 1536 , 1540 , 1560 , 1568 , 1575 , 1584 , 1600 , 1617 , 1620 ,
    1625 , 1638 , 1650 , 1664 , 1680 , 1701 , 1715 , 1728 , 1750 , 1755 ,
    1760 , 1764 , 1782 , 1792 , 1800 , 1820 , 1848 , 1872 , 1875 , 1890 ,
    1911 , 1920 , 1925 , 1944 , 1950 , 1960 , 1980 , 2000 , 2016 , 2025 ,
    2048 , 2058 , 2079 , 2080 , 2100 , 2106 , 2112 , 2156 , 2160 , 2184 ,
    2187 , 2200 , 2205 , 2240 , 2250 , 2268 , 2275 , 2304 , 2310 , 2340 ,
    2352 , 2376 , 2400 , 2401 , 2430 , 2450 , 2457 , 2464 , 2475 , 2496 ,
    2500 , 2520 , 2548 , 2560 , 2592 , 2600 , 2625 , 2640 , 2646 , 2673 ,
    2688 , 2695 , 2700 , 2730 , 2744 , 2750 , 2772 , 2800 , 2808 , 2816 ,
    2835 , 2880 , 2912 , 2916 , 2925 , 2940 , 2970 , 3000 , 3024 , 3072 ,
    3080 , 3087 , 3120 , 3125 , 3136 , 3150 , 3159 , 3168 , 3185 , 3200 ,
    3234 , 3240 , 3250 , 3276 , 3300 , 3328 , 3360 , 3375 , 3402 , 3430 ,
    3456 , 3465 , 3500 , 3510 , 3520 , 3528 , 3564 , 3584 , 3600 , 3640 ,
    3645 , 3675 , 3696 , 3744 , 3750 , 3773 , 3780 , 3822 , 3840 , 3850 ,
    3888 , 3900 , 3920 , 3960 , 3969 , 4000 , 4032 , 4050 , 4095 , 4096 ,
    4116 , 4125 , 4158 , 4160 , 4200 , 4212 , 4224 , 4312 , 4320 , 4368 ,
    4374 , 4375 , 4400 , 4410 , 4455 , 4459 , 4480 , 4500 , 4536 , 4550 ,
    4608 , 4620 , 4680 , 4704 , 4725 , 4752 , 4800 , 4802 , 4851 , 4860 ,
    4875 , 4900 , 4914 , 4928 , 4950 , 4992 , 5000 , 5040 , 5096 , 5103 ,
    5120 , 5145 , 5184 , 5200 , 5250 , 5265 , 5280 , 5292 , 5346 , 5376 ,
    5390 , 5400 , 5460 , 5488 , 5500 , 5544 , 5600 , 5616 , 5625 , 5632 ,
    5670 , 5733 , 5760 , 5775 , 5824 , 5832 , 5850 , 5880 , 5940 , 6000 ,
    6048 , 6075 , 6125 , 6144 , 6160 , 6174 , 6237 , 6240 , 6250 , 6272 ,
    6300 , 6318 , 6336 , 6370 , 6400 , 6468 , 6480 , 6500 , 6552 , 6561 ,
    6600 , 6615 , 6656 , 6720 , 6750 , 6804 , 6825 , 6860 , 6875 , 6912 ,
    6930 , 7000 , 7020 , 7040 , 7056 , 7128 , 7168 , 7200 , 7203 , 7280 ,
    7290 , 7350 , 7371 , 7392 , 7425 , 7488 , 7500 , 7546 , 7560 , 7644 ,
    7680 , 7700 , 7776 , 7800 , 7840 , 7875 , 7920 , 7938 , 8000 , 8019 ,
    8064 , 8085 , 8100 , 8125 , 8190 , 8192 , 8232 , 8250 , 8316 , 8320 ,
    8400 , 8424 , 8448 , 8505 , 8575 , 8624 , 8640 , 8736 , 8748 , 8750 ,
    8775 , 8800 , 8820 , 8910 , 8918 , 8960 , 9000 , 9072 , 9100 , 9216 ,
    9240 , 9261 , 9360 , 9375 , 9408 , 9450 , 9477 , 9504 , 9555 , 9600 ,
    9604 , 9625 , 9702 , 9720 , 9750 , 9800 , 9828 , 9856 , 9900 , 9984 ,
    10000, 10080, 10125, 10192, 10206, 10240, 10290, 10368, 10395, 10400,
    10500, 10530, 10560, 10584, 10692, 10752, 10780, 10800, 10920, 10935,
    10976, 11000, 11025, 11088, 11200, 11232, 11250, 11264, 11319, 11340,
    11375, 11466, 11520, 11550, 11648, 11664, 11700, 11760, 11880, 11907,
    12000, 12005, 12096, 12150, 12250, 12285, 12288, 12320, 12348, 12375,
    12474, 12480, 12500, 12544, 12600, 12636, 12672, 12740, 12800, 12936,
    12960, 13000, 13104, 13122, 13125, 13200, 13230, 13312, 13365, 13377,
    13440, 13475, 13500, 13608, 13650, 13720, 13750, 13824, 13860, 14000,
    14040, 14080, 14112, 14175, 14256, 14336, 14400, 14406, 14553, 14560,
    14580, 14625, 14700, 14742, 14784, 14850, 14976, 15000, 15092, 15120,
    15288, 15309, 15360, 15400, 15435, 15552, 15600, 15625, 15680, 15750,
    15795, 15840, 15876, 15925, 16000, 16038, 16128, 16170, 16200, 16250,
    16380, 16384, 16464, 16500, 16632, 16640, 16800, 16807, 16848, 16875,
    16896, 17010, 17150, 17199, 17248, 17280, 17325, 17472, 17496, 17500,
    17550, 17600, 17640, 17820, 17836, 17920, 18000, 18144, 18200, 18225,
    18375, 18432, 18480, 18522, 18711, 18720, 18750, 18816, 18865, 18900,
    18954, 19008, 19110, 19200, 19208, 19250, 19404, 19440, 19500, 19600,
    19656, 19683, 19712, 19800, 19845, 19968, 20000, 20160, 20250, 20384,
    20412, 20475, 20480, 20580, 20625, 20736, 20790, 20800, 21000, 21060,
    21120, 21168, 21384, 21504, 21560, 21600, 21609, 21840, 21870, 21875,
    21952, 22000, 22050, 22113, 22176, 22275, 22295, 22400, 22464, 22500,
    22528, 22638, 22680, 22750, 22932, 23040, 23100, 23296, 23328, 23400,
    23520, 23625, 23760, 23814, 24000, 24010, 24057, 24192, 24255, 24300,
    24375, 24500, 24570, 24576, 24640, 24696, 24750, 24948, 24960, 25000,
    25088, 25200, 25272, 25344, 25480, 25515, 25600, 25725, 25872, 25920,
    26000, 26208, 26244, 26250, 26325, 26400, 26411, 26460, 26624, 26730,
    26754, 26880, 26950, 27000, 27216, 27300, 27440, 27500, 27648, 27720,
    27783, 28000, 28080, 28125, 28160, 28224, 28350, 28431, 28512, 28665,
    28672, 28800, 28812, 28875, 29106, 29120, 29160, 29250, 29400, 29484,
    29568, 29700, 29952, 30000, 30184, 30240, 30375, 30576, 30618, 30625,
    30720, 30800, 30870, 31104, 31185, 31200, 31213, 31250, 31360, 31500,
    31590, 31680, 31752, 31850, 32000, 32076, 32256, 32340, 32400, 32500,
    32760, 32768, 32805, 32928, 33000, 33075, 33264, 33280, 33600, 33614,
    33696, 33750, 33792, 33957, 34020, 34125, 34300, 34375, 34398, 34496,
    34560, 34650, 34944, 34992, 35000, 35100, 35200, 35280, 35640, 35672,
    35721, 35840, 36000, 36015, 36288, 36400, 36450, 36750, 36855, 36864,
    36960, 37044, 37125, 37422, 37440, 37500, 37632, 37730, 37800, 37908,
    38016, 38220, 38400, 38416, 38500, 38808, 38880, 39000, 39200, 39312,
    39366, 39375, 39424, 39600, 39690, 39936, 40000, 40095, 40131, 40320,
    40425, 40500, 40625, 40768, 40824, 40950, 40960, 41160, 41250, 41472,
    41580, 41600, 42000, 42120, 42240, 42336, 42525, 42768, 42875, 43008,
    43120, 43200, 43218, 43659, 43680, 43740, 43750, 43875, 43904, 44000,
    44100, 44226, 44352, 44550, 44590, 44800, 44928, 45000, 45056, 45276,
    45360, 45500, 45864, 45927, 46080, 46200, 46305, 46592, 46656, 46800,
    46875, 47040, 47250, 47385, 47520, 47628, 47775, 48000, 48020, 48114,
    48125, 48384, 48510, 48600, 48750, 49000, 49140, 49152, 49280, 49392,
    49500, 49896, 49920, 50000, 50176, 50400, 50421, 50544, 50625, 50688,
    50960, 51030, 51200, 51450, 51597, 51744, 51840, 51975, 52000, 52416,
    52488, 52500, 52650, 52800, 52822, 52920, 53248, 53460, 53508, 53760,
    53900, 54000, 54432, 54600, 54675, 54880, 55000, 55125, 55296, 55440,
    55566, 56000, 56133, 56160, 56250, 56320, 56448, 56595, 56700, 56862,
    56875, 57024, 57330, 57344, 57600, 57624, 57750, 58212, 58240, 58320,
    58500, 58800, 58968, 59049, 59136, 59400, 59535, 59904, 60000, 60025,
    60368, 60480, 60750, 61152, 61236, 61250, 61425, 61440, 61600, 61740,
    61875, 62208, 62370, 62400, 62426, 62500, 62720, 63000, 63180, 63360,
    63504, 63700, 64000, 64152, 64512, 64680, 64800, 64827, 65000, 65520,
    65536, 65610, 65625, 65856, 66000, 66150, 66339, 66528, 66560, 66825,
    66885, 67200, 67228, 67375, 67392, 67500, 67584, 67914, 68040, 68250,
    68600, 68750, 68796, 68992, 69120, 69300, 69888, 69984, 70000, 70200,
    70400, 70560, 70875, 71280, 71344, 71442, 71680, 72000, 72030, 72171,
    72576, 72765, 72800, 72900, 73125, 73500, 73710, 73728, 73920, 74088,
    74250, 74844, 74880, 75000, 75264, 75460, 75600, 75816, 76032, 76440,
    76545, 76800, 76832, 77000, 77175, 77616, 77760, 78000, 78125, 78400,
    78624, 78732, 78750, 78848, 78975, 79200, 79233, 79380, 79625, 79872,
    80000, 80190, 80262, 80640, 80850, 81000, 81250, 81536, 81648, 81900,
    81920, 82320, 82500, 82944, 83160, 83200, 83349, 84000, 84035, 84240,
    84375, 84480, 84672, 85050, 85293, 85536, 85750, 85995, 86016, 86240,
    86400, 86436, 86625, 87318, 87360, 87480, 87500, 87750, 87808, 88000,
    88200, 88452, 88704, 89100, 89180, 89600, 89856, 90000, 90112, 90552,
    90720, 91000, 91125, 91728, 91854, 91875, 92160, 92400, 92610, 93184,
    93312, 93555, 93600, 93639, 93750, 94080, 94325, 94500, 94770, 95040,
    95256, 95550, 96000, 96040, 96228, 96250, 96768, 97020, 97200, 97500,
    98000, 98280, 98304, 98415, 98560, 98784, 99000, 99225, 99792, 99840,
    100000};
