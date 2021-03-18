#include "LocalizedNrField.hpp"
#include "LocalizedNrFieldPointwise.hpp"
#include "LocalizedNrFieldRefinement.hpp"
#include "GeneralNrField.hpp"
#include "geodesy.hpp"
#include "numerical.hpp"
#include "inparam.hpp"

void LocalizedNrField::buildInparam(std::vector<std::unique_ptr<LocalizedNrField>> &localFields,
    const std::shared_ptr<const GeneralNrField> &generalField, double distTol) {
    
    if (!generalField) { // pointwise_local
        // need to think about how to structure this data...
        const std::string &fname =
        inparam::gInparamNr.get<std::string>("pointwise_local:nc_data_file");
        localFields.push_back(std::make_unique<LocalizedNrFieldPointwise>(fname, distTol));
        
    } else { // refinement windows or nothing
        int winCount = inparam::gInparamNr.get<int>("Nr_windows:list_of_refinement_windows:[?]");
    
        // loop over list
        for (int windex = 0; windex < winCount; windex++) {
            // window name and key in inparam
            std::string keyInparam = "Nr_windows:list_of_refinement_windows";
            const std::string &winName = inparam::gInparamNr.
            get<std::string>(keyInparam + ":{" + bstring::toString(windex) + "}");
            keyInparam += ":[" + bstring::toString(windex) + "]:" + winName;
            
            // activated or not
            if (!inparam::gInparamNr.get<bool>(keyInparam + ":activated")) {
                continue;
            }
            
            std::string cs = keyInparam + ":coordinates:";
            
            const std::vector<double> &lat =inparam::gInparamNr.getVector<double>(cs + "radial_range");
            const std::vector<double> &lon =inparam::gInparamNr.getVector<double>(cs + "azimuthal_range");
            const std::vector<double> &depth =inparam::gInparamNr.getVector<double>(cs + "vertical_range");
            
            if (lat.size() != 2 || lon.size() != 2 || depth.size() != 2) {
                throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                         "coordinate ranges for refinement windows ||"
                                         "must be given as vectors of 2. ||"
                                         "Refinement name:" + keyInparam);
            }
            
            eigen::DMatX3 crds;
            crds << lat[0], lon[0], depth[0],
                    lat[1], lon[1], depth[1];
            
            // depth unit conversion
            const std::string &lu = inparam::gInparamNr.get<std::string>(cs + "length_unit");
            if (lu == "km") {
                crds.col(2) *= 1000;
            } else if (lu == "m") {
                // nothing
            } else {
                throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                         "unknown length unit:" + lu + " ||"
                                         "Refinement name:" + keyInparam);
            }
            
            // depth to radius
            const std::string &typev = inparam::gInparamNr.get<std::string>(cs + "vertical_type");
            if (typev == "DEPTH") {
                bool solidSurf = false;
                if (inparam::gInparamNr.contains(cs + "depth_below_solid_surface")) {
                    solidSurf = inparam::gInparamNr.get<bool>(cs + "depth_below_solid_surface");
                }
                if (solidSurf) {
                    crds.col(2) = geodesy::getOuterSolidRadius() - crds.col(2).reverse().array();
                } else {
                    crds.col(2) = geodesy::getOuterRadius() - crds.col(2).reverse().array();
                }
            } else if (typev == "RADIUS") {
                // nothing
            } else {
                throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                         "unknown vertical coordinate type:" + typev + " ||"
                                         "Refinement name:" + keyInparam);
            }
            
            // lat lon to s phi
            bool ellipticity;
            bool geographic;
            const std::string &typeh = inparam::gInparamNr.get<std::string>(cs + "horizontal_type");
            if (typeh == "LATITUDE_LONGITUDE") {
                bool geographic = true;
                bool ellipticity = false;
                if (inparam::gInparamNr.contains(cs + "ellipticity")) {
                    ellipticity = inparam::gInparamNr.get<bool>(cs + "ellipticity");
                }
            } else if (typeh == "DISTANCE_AZIMUTH") {
                // more unit conversions
                if (geodesy::isCartesian() && lu == "km") crds.col(0) *= 1000;
                
                const std::string &typea = inparam::gInparamNr.get<std::string>(cs + "angle_unit");
                if (typea == "degree") {
                    crds.col(1) *= numerical::dDegree;
                    if (!geodesy::isCartesian()) crds.col(0) *= numerical::dDegree;
                } else if (typea == "radian") {
                    // nothing
                } else {
                    throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                             "unknown angular unit:" + typea + " ||"
                                             "Refinement name:" + keyInparam);
                }
            } else {
                throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                         "unknown horizontal coordinate type:" + typeh + " ||"
                                         "Refinement name:" + keyInparam);
            }
            
            const std::string &type_Nr = inparam::gInparamNr.get<std::string>(keyInparam + ":type");
            
            double nr;
            bool relative = false;
            if (type_Nr == "RELATIVE") {
                nr = inparam::gInparamNr.get<double>(keyInparam + ":value");
                relative = true;
            } else {
                nr = (double)inparam::gInparamNr.get<int>(keyInparam + ":value");
                
                if (type_Nr == "GLOBAL") {
                    double dphi = crds(1, 1) - crds(0, 1);
                    if (dphi < 0) dphi += 2 * numerical::dPi;
                    nr = nr * dphi / (2 * numerical::dPi);
                    
                } else if (type_Nr == "LOCAL") {
                    // nothing
                } else {
                    throw std::runtime_error("LocalisedNrField::buildInparam ||"
                                             "unknown local Nr type:" + type_Nr + " ||"
                                             "Refinement name:" + keyInparam);
                }
            }
            
            localFields.push_back(std::make_unique<LocalizedNrFieldRefinement>(generalField, crds, nr, relative, geographic, ellipticity));
        }
    }
}
