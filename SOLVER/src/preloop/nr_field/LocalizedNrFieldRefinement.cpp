
#include "LocalizedNrFieldRefinement.hpp"
#include "geodesy.hpp"
#include "LocalizedNr.hpp"
#include "GeneralNrField.hpp"

LocalizedNrFieldRefinement::LocalizedNrFieldRefinement(const std::shared_ptr<const GeneralNrField> &gnr, 
    eigen::DMatX3 &coords, double &value, bool relative, bool geographic, bool ellipticity):
mGeneralNrField(gnr), mCoords(coords), mValue(value), mIsRelative(relative),
mGeographic(geographic), mEllipticity(ellipticity) {
    // nothing
}

LocalizedNr LocalizedNrFieldRefinement::getWindowsAtPoint(const eigen::DCol2 &sz) const {
    double nr_bg = (double)mGeneralNrField->getNrAtPoint(sz);
    
    eigen::DCol2 crds = sz;
    if (!geodesy::isCartesian()) {
        geodesy::sz2rtheta(crds, false);
        if (mGeographic) {
            eigen::DRow3 tpr;
            tpr << crds(0), 0., crds(1);
            geodesy::spz2llr(tpr, mEllipticity, false);
            crds(0) = tpr(0);
            crds(1) = tpr(2);
        }
    }
    
    if (mCoords(0, 0) < crds(0) && crds(0) < mCoords(1, 0) && mCoords(0, 2) < crds(1) && crds(1) < mCoords(1, 2)) {
        // s,z in windowed range
        double phi_start = mCoords(0, 1);
        double phi_end = mCoords(1, 1);
        if (mGeographic) { // not source centered -> need to calculate angle based on coords
            eigen::DRow3 crds_win_start, crds_win_end;
            crds_win_start << crds(0), mCoords(0, 1), crds(1);
            crds_win_end << crds(0), mCoords(1, 1), crds(1);
            geodesy::llr2spz(crds_win_start, mEllipticity);
            geodesy::llr2spz(crds_win_end, mEllipticity);
            phi_start = crds_win_start(1);
            phi_end = crds_win_end(1);
        }
        
        std::vector<std::pair<double, double>> phiNr{std::make_pair(phi_start, mValue), std::make_pair(phi_end, 0.)};
        LocalizedNr wins(phiNr);

        if (mIsRelative) {
            int nr = mValue * nr_bg;
            wins.setNrUpscaled(0, nr);
        }
        
        wins.setNrUpscaled(1, nr_bg);
        
        return wins;
    } else {
        // s,z not in windowed range
        std::vector<std::pair<double, double>> phiNr(1, std::make_pair(0., nr_bg));
        return LocalizedNr(phiNr);
    }
}

std::string LocalizedNrFieldRefinement::verbose() const {
    return "todo: refinement verbose.";
}
