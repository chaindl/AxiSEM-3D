//
//  GeneralNrField.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/14/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of Nr(s,z)

#include "inparam.hpp"
#include "GeneralNrFieldConstant.hpp"
#include "GeneralNrFieldAnalytical.hpp"
#include "GeneralNrFieldPointwise.hpp"
#include "GeneralNrFieldStructured.hpp"
#include "timer.hpp"

// build from inparam
std::shared_ptr<const GeneralNrField> GeneralNrField::
buildInparam(const double distTolerance) {
    // read type and lucky number
    const std::string &type = inparam::gInparamNr.
    getWithLimits<std::string>("type_Nr", {
        {"CONSTANT", "CONSTANT"},
        {"ANALYTICAL", "ANALYTICAL"},
        {"POINTWISE", "POINTWISE"},
        {"STRUCTURED", "STRUCTURED"},
        {"POINTWISE_LOCAL", "POINTWISE_LOCAL"}});
    
    if (type == "POINTWISE_LOCAL") return nullptr;
    
    // type-dependent
    timer::gPreloopTimer.begin("Building Nr(s,z) of type " + type);
    std::unique_ptr<const GeneralNrField> GeneralNrField;
    if (type == "CONSTANT") {
        int nr = inparam::gInparamNr.getWithBounds("constant", 1);
        GeneralNrField = std::make_unique<const GeneralNrFieldConstant>(nr);
        
    } else if (type == "ANALYTICAL") {
        GeneralNrField = std::make_unique<const GeneralNrFieldAnalytical>();
        
    } else if (type == "POINTWISE") {
        const std::string &fname =
        inparam::gInparamNr.get<std::string>("pointwise:nc_data_file");
        const double &factor =
        inparam::gInparamNr.get<double>("pointwise:multip_factor");
        GeneralNrField = std::make_unique<const GeneralNrFieldPointwise>
        (fname, factor, distTolerance);
        
    } else if (type == "STRUCTURED") {
        const std::string &fname =
        inparam::gInparamNr.get<std::string>("structured:nc_data_file");
        int valOOR =
        inparam::gInparamNr.getWithBounds("structured:value_out_of_range", 0);
        GeneralNrField = std::make_unique<const GeneralNrFieldStructured>(fname, valOOR);
    } else {
        // impossible
        throw std::runtime_error("NrField::buildInparam || Impossible.");
    }
    timer::gPreloopTimer.ended("Building Nr(s,z) of type " + type);
    return GeneralNrField;
}
