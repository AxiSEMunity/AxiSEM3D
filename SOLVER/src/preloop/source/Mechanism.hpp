//
//  Mechanism.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 5/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  base class of source mechanism

#ifndef Mechanism_hpp
#define Mechanism_hpp

#include "eigen_sem.hpp"
#include "bstring.hpp"

class Quad;
class STF;
class Domain;

namespace eigen {
    using numerical::ComplexD;
    // coords
    typedef Eigen::Matrix<double, 1, 3> DRow3;
    typedef Eigen::Matrix<double, 3, 3> DMat33;
    // source patterns
    typedef Eigen::Matrix<ComplexD, Eigen::Dynamic, nPEM> ZMatXN;
    typedef Eigen::Matrix<ComplexD, Eigen::Dynamic, nPEM * 3> ZMatXN3;
    typedef Eigen::Matrix<ComplexD, Eigen::Dynamic, nPEM * 6> ZMatXN6;
}

class Mechanism {
public:
    // supported source mechanisms
    enum class SM_Type {MomentTensor, ForceVector, FluidPressure};
    inline static const std::map<SM_Type, std::string> sSM_TypeStr = {
        {SM_Type::MomentTensor, "moment tensor"},
        {SM_Type::ForceVector, "force vector"},
        {SM_Type::FluidPressure, "fluid pressure"}};
    inline static const std::map<SM_Type, int> sSM_TypeDim = {
        {SM_Type::MomentTensor, 6},
        {SM_Type::ForceVector, 3},
        {SM_Type::FluidPressure, 1}};
    
    // constructor
    Mechanism(const SM_Type &type, const eigen::DColX &data):
    mType(type), mData(data) {
        // verify
        if (sSM_TypeDim.at(mType) != mData.size()) {
            throw std::runtime_error
            ("Mechanism::Mechanism || "
             "Invalid data size for " + sSM_TypeStr.at(mType) + ": || "
             "Required = " + bstring::toString(sSM_TypeDim.at(mType)) + " || "
             "Provided = " + bstring::toString(mData.size()));
        }
    }
    
    // solid or fluid
    bool inFluid() const {
        return mType == SM_Type::FluidPressure;
    }
    
    // release element source
    void
    release(const eigen::DMat33 &Qzsp, bool sourceOnAxis,
            const eigen::DRowN &inplaneFactor, double phi,
            const Quad &quad, std::unique_ptr<STF> &stf, Domain &domain) const;
    
    // verbose
    std::string verbose(int keyWidth) const {
        std::stringstream ss;
        ss << bstring::boxSubTitle(2, "Source mechanism");
        ss << bstring::boxEquals(4, keyWidth, sSM_TypeStr.at(mType),
                                 std::vector<double>
                                 (mData.data(), mData.data() + mData.size()));
        return ss.str();
    }
    
private:
    // type
    const SM_Type mType;
    // data
    const eigen::DColX mData;
    
    ////////////////////////////// static //////////////////////////////
public:
    // build from inparam
    static std::unique_ptr<const Mechanism>
    buildInparam(int sindex, const std::string &sourceName);
};

#endif /* Mechanism_hpp */
