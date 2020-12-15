//
//  ABC.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/30/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  parameters for absorbing boundary condition

#include "ABC.hpp"
#include "inparam.hpp"
#include "vector_tools.hpp"
#include "ExodusMesh.hpp"

// build from inparam
std::unique_ptr<ABC> ABC::buildInparam(const ExodusMesh &exodusMesh) {
    const InparamYAML &gm = inparam::gInparamModel;
    std::string root = "absorbing_boundary";
    
    // create
    std::unique_ptr<ABC> abc = std::make_unique<ABC>();
    
    // keys
    abc->mUserKeys = gm.getVector<std::string>(root + ":boundaries");
    if (abc->mUserKeys.size() == 0) {
        return abc;
    }
    
    // check keys
    if (abc->mUserKeys.size() > 3) {
        throw std::runtime_error("ABC::buildInparam || "
                                 "Parameter absorbing_boundary:boundaries || "
                                 "cannot contain more than three entities.");
    }
    std::vector<std::string> uniqueKeys = abc->mUserKeys;
    vector_tools::sortUnique(uniqueKeys);
    if (uniqueKeys.size() != abc->mUserKeys.size()) {
        throw std::runtime_error("ABC::buildInparam || "
                                 "Parameter absorbing_boundary:boundaries || "
                                 "cannot contain duplicated entities.");
    }
    for (auto it = abc->mUserKeys.begin(); it != abc->mUserKeys.end(); ++it) {
        if (*it != "RIGHT" && *it != "BOTTOM" && *it != "TOP") {
            throw
            std::runtime_error("ABC::buildInparam || "
                               "Parameter absorbing_boundary:boundaries || "
                               "contains an invalid entity: " + *it);
        }
    }
    
    // Clayton
    abc->mClayton = gm.get<bool>(root + ":enable_Clayton_Enquist");
    
    // sponge
    root += ":Kosloff_Kosloff";
    abc->mSponge = gm.get<bool>(root + ":enable");
    if (!abc->mClayton && !abc->mSponge) {
        throw std::runtime_error("ABC::buildInparam || "
                                 "Both Clayton-Enquist and Kosloff-Kosloff "
                                 "are disabled.");
    }
    
    // get parameters for sponge
    std::vector<double> relSpan;
    if (abc->mSponge) {
        relSpan = gm.getVector<double>(root + ":relative_spans");
        // check sizes
        if (abc->mUserKeys.size() != relSpan.size()) {
            throw std::runtime_error
            ("ABC::buildInparam || "
             "The number of Kosloff_Kosloff::relative_spans || "
             "does not match the number of boundaries.");
        }
        // check values
        if (*std::min_element(relSpan.begin(), relSpan.end()) < .01 ||
            *std::max_element(relSpan.begin(), relSpan.end()) > .25) {
            throw std::runtime_error("ABC::buildInparam || "
                                     "Kosloff_Kosloff:relative_spans "
                                     "must range between 0.01 and 0.25.");
        }
        // read expr string
        abc->mGammaExprSolidStr =
        gm.get<std::string>(root + ":gamma_expr_solid");
        abc->mGammaExprFluidStr =
        gm.get<std::string>(root + ":gamma_expr_fluid");
    }
    
    // setup in mesh
    for (int ikey = 0; ikey < abc->mUserKeys.size(); ikey++) {
        const std::string &key = abc->mUserKeys[ikey];
        // get boundary info
        double outer = 0., meshSpan = 0.;
        if (exodusMesh.boundaryInfoABC(key, outer, meshSpan)) {
            abc->mBoundaryKeys.push_back(key);
            if (abc->mSponge) {
                double span = relSpan[ikey] * meshSpan;
                abc->mSpongeOuterSpan.insert({key, {outer, span}});
            }
        }
    }
    
    // mesh and expressions for sponge
    if (!abc->sponge()) {
        return abc;
    }
    
    // setup expressions
    // constant
    abc->mT0 = exodusMesh.getGlobalVariable("minimum_period");
    // symbol table
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("VP", sVP);
    symbol_table.add_variable("RHO", sRHO);
    symbol_table.add_variable("SPAN", sSPAN);
    symbol_table.add_constant("T0", abc->mT0);
    // fluid
    abc->mGammaExprFluid.register_symbol_table(symbol_table);
    // solid
    symbol_table.add_variable("VS", sVS);
    abc->mGammaExprSolid.register_symbol_table(symbol_table);
    // parse
    exprtk::parser<double> parser;
    if (!parser.compile(abc->mGammaExprSolidStr, abc->mGammaExprSolid)) {
        throw std::runtime_error("ABC::buildInparam || Error parsing "
                                 "Kosloff_Kosloff::gamma_expr_solid.");
    }
    if (!parser.compile(abc->mGammaExprFluidStr, abc->mGammaExprFluid)) {
        throw std::runtime_error("ABC::buildInparam || Error parsing "
                                 "Kosloff_Kosloff::gamma_expr_fluid.");
    }
    return abc;
}

// verbose
std::string ABC::verbose() const {
    using namespace bstring;
    std::stringstream ss;
    ss << boxTitle("Absorbing Boundary");
    
    // user keys
    if (mUserKeys.size() == 0) {
        ss << "* Absorbing boundary has been disabled.\n";
        ss << boxBaseline() << "\n\n";
        return ss.str();
    }
    ss << boxEquals(0, 25, "user-specified boundaries", mUserKeys, "=", true);
    
    // mesh keys
    if (mBoundaryKeys.size() == 0) {
        ss << "* The mesh contains none of these boundaries.\n";
        ss << boxBaseline() << "\n\n";
        return ss.str();
    }
    ss << boxEquals(0, 25, "those contained in mesh", mBoundaryKeys, "=", true);
    
    // methods
    ss << boxEquals(0, 25, "Clayton-Enquist enabled", mClayton);
    ss << boxEquals(0, 25, "Kosloff-Kosloff enabled", mSponge);
    
    // sponge
    if (mSponge) {
        ss << boxSubTitle(0, "Parameters for Kosloff-Kosloff");
        for (const std::string &key: mBoundaryKeys) {
            // data
            double outer = std::get<0>(mSpongeOuterSpan.at(key));
            double span = std::get<1>(mSpongeOuterSpan.at(key));
            // verbose
            ss << "  * " << key << ":\n";
            ss << boxEquals(4, 21, "boundary location", outer);
            ss << boxEquals(4, 21, "span of sponge layer", std::abs(span));
            ss << boxEquals(4, 21, "range of sponge layer",
                            range(outer - span, outer));
        }
        ss << "  * Expression for γ-factor" << ":\n";
        ss << "    in solid: " + mGammaExprSolidStr + "\n";
        ss << "    in fluid: " + mGammaExprFluidStr + "\n";
    }
    ss << boxBaseline() << "\n\n";
    return ss.str();
}

// get U0 solid
eigen::DColX ABC::getU0Solid(double span,
                             const eigen::DColX &vp,
                             const eigen::DColX &vs,
                             const eigen::DColX &rho) const {
    // initialize
    int nr = (int)vp.rows();
    eigen::DColX U0 = eigen::DColX::Zero(nr, 1);
    // compute by expr
    sSPAN = span;
    for (int alpha = 0; alpha < nr; alpha++) {
        sVP = vp(alpha);
        sVS = vs(alpha);
        sRHO = rho(alpha);
        U0(alpha) = mGammaExprSolid.value();
    }
    return U0;
}

// get U0 fluid
eigen::DColX ABC::getU0Fluid(double span,
                             const eigen::DColX &vp,
                             const eigen::DColX &rho) const {
    // initialize
    int nr = (int)vp.rows();
    eigen::DColX U0 = eigen::DColX::Zero(nr, 1);
    // compute by expr
    sSPAN = span;
    for (int alpha = 0; alpha < nr; alpha++) {
        sVP = vp(alpha);
        sRHO = rho(alpha);
        U0(alpha) = mGammaExprFluid.value();
    }
    return U0;
}
