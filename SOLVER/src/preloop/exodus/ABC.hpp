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
class ExodusMesh;

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
    
    // get sponge data
    std::tuple<double, double, double>
    getSpongeData(const std::string &key) const {
        return mSpongeData.at(key);
    }
    
private:
    // data from inparam
    std::vector<std::string> mUserKeys;
    bool mClayton = false;
    bool mSponge = false;
    
    // data based on mesh
    std::vector<std::string> mBoundaryKeys;
    std::map<std::string, std::tuple<double, double, double>> mSpongeData;
};

#endif /* ABC_hpp */
