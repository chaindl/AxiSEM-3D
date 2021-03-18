#include "LocalizedNrFieldPointwise.hpp"
#include "LocalizedNr.hpp"

LocalizedNrFieldPointwise::LocalizedNrFieldPointwise(const std::string &fname, double distTol) {
    // nothing
};
  
// get nr by (s, z)
LocalizedNr LocalizedNrFieldPointwise::getWindowsAtPoint(const eigen::DCol2 &sz) const {
    std::vector<std::pair<double, double>> phiNr(1, std::make_pair(0., 0.));
    return LocalizedNr(phiNr);
};
    
// verbose
std::string LocalizedNrFieldPointwise::verbose() const {
    return "Pointwise local not yet implemented.";
};
