//
//  ABC.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/30/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  parameters for absorbing boundary condition

#ifndef ABC_hpp
#define ABC_hpp

#include <string>
#include <vector>
#include <map>
#include <memory>

#include "ExodusMesh.hpp"
#include "vector_tools.hpp"

// exprtk
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "exprtk/exprtk.hpp"
#pragma clang diagnostic pop

class ABC {
public:
    // build from inparam
    static std::unique_ptr<ABC> buildInparam(const ExodusMesh &exodusMesh);
    
    // verbose
    std::string verbose() const;
    
    // get
    const std::vector<std::string> &getBoundaryKeys() const {
        return mBoundaryKeys;
    }
    
    // get
    bool clayton() const {
        return mClayton;
    }
    
    // get
    bool sponge() const {
        return mSponge;
    }
    
    // get sponge outer and span
    const std::tuple<double, double> &
    getSpongeOuterSpan(const std::string &key) const {
        return mSpongeOuterSpan.at(key);
    }
    
    // get U0 solid
    eigen::DColX getU0Solid(double span,
                            const eigen::DColX &vp,
                            const eigen::DColX &vs,
                            const eigen::DColX &rho) const;
    
    // get U0 fluid
    eigen::DColX getU0Fluid(double span,
                            const eigen::DColX &vp,
                            const eigen::DColX &rho) const;
    
private:
    /////// general data ///////
    std::vector<std::string> mUserKeys;
    std::vector<std::string> mBoundaryKeys;
    bool mClayton = false;
    bool mSponge = false;
    
    /////// data for sponge ///////
    // geometry
    std::map<std::string, std::tuple<double, double>> mSpongeOuterSpan;
    
    // gamma expressions
    // using static for variables to keep const modifier
    inline static double sVP = 0.;
    inline static double sVS = 0.;
    inline static double sRHO = 0.;
    inline static double sSPAN = 0.;
    double mT0 = 0.;
    std::string mGammaExprSolidStr;
    std::string mGammaExprFluidStr;
    exprtk::expression<double> mGammaExprSolid;
    exprtk::expression<double> mGammaExprFluid;
};

#endif /* ABC_hpp */
