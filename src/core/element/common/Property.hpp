//
//  Property.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 8/14/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  elemental properties
//  * properties such as lambda and mu for elasticity
//  * handling 1D and 3D properties in the same manner (via FieldArithmetic)

#ifndef Property_hpp
#define Property_hpp

#include "numerical.hpp"
#include "eigen.hpp"
#include <memory>

template <int P, int N = P * P>
class Property {
    typedef Eigen::Matrix<double, P, P, Eigen::RowMajor> DMatPP_RM;
    typedef Eigen::Matrix<numerical::Real, P, P, Eigen::RowMajor> RMatPP_RM;
    typedef Eigen::Matrix<double, Eigen::Dynamic, N> DMatXN;
    typedef Eigen::Matrix<numerical::Real, Eigen::Dynamic, N> RMatXN;
    
public:
    // default constructor (workspace)
    Property() = default;
    
    // 1D constructor
    Property(const DMatPP_RM &data1D):
    mData1D(std::make_unique<RMatPP_RM>
            (data1D.template cast<numerical::Real>())) {
        // nothing
    }
    
    // 3D constructor
    Property(const DMatXN &data3D):
    mData3D(data3D.template cast<numerical::Real>()) {
        // nothing
    }
    
    // copy constructor
    Property(const Property &other):
    mData1D((other.mData1D == nullptr) ? nullptr :
            std::make_unique<RMatPP_RM>(*(other.mData1D))),
    mData3D(other.mData3D) {
        // nothing
    }
    
    // check compatibility
    void checkCompatibility(bool isProp1D, int nr, bool elemInFourier,
                            const std::string &propName,
                            const std::string &ownerName) const {
        if (isProp1D) {
            // 1D/3D
            if (mData1D == nullptr) {
                throw std::runtime_error("Property::checkCompatibility || "
                                         "Incompatible 1D/3D flags. || "
                                         "Property name: " + propName + " || "
                                         "Owner class: " + ownerName);
            }
        } else {
            // 1D/3D
            if (mData3D.rows() == 0 || elemInFourier) {
                throw std::runtime_error("Property::checkCompatibility || "
                                         "Incompatible 1D/3D flags. || "
                                         "Property name: " + propName + " || "
                                         "Owner class: " + ownerName);
            }
            // size
            if (mData3D.rows() != nr) {
                throw std::runtime_error("Property::checkCompatibility || "
                                         "Incompatible sizes in 3D. || "
                                         "Property name: " + propName + " || "
                                         "Owner class: " + ownerName);
            }
        }
    }
    
    // get 1D
    inline const RMatPP_RM &get1D() const {
        return *mData1D;
    }
    
    // get 3D
    inline const RMatXN &get3D() const {
        return mData3D;
    }
    
    ////////////////// for workspace //////////////////
    // expand workspace
    void expandWorkspace(bool isProp1D, int nr) {
        if (isProp1D) {
            if (mData1D == nullptr) {
                mData1D = std::make_unique<RMatPP_RM>();
            }
        } else {
            if (mData3D.rows() < nr) {
                mData3D.resize(nr, N);
            }
        }
    }
    
    // this = other * factor
    void set(bool isProp1D, const Property &other, numerical::Real factor) {
        if (isProp1D) {
            *mData1D = *(other.mData1D) * factor;
        } else {
            mData3D.topRows(other.mData3D.rows()) = other.mData3D * factor;
        }
    }
    
    // this = other * factor
    void set(bool isProp1D, const Property &other, const Property &factor) {
        if (isProp1D) {
            *mData1D = other.mData1D->cwiseProduct(*(factor.mData1D));
        } else {
            mData3D.topRows(other.mData3D.rows()) =
            other.mData3D.cwiseProduct(factor.mData3D);
        }
    }
    
private:
    // 1D data
    std::unique_ptr<RMatPP_RM> mData1D = nullptr;
    
    // 3D data
    RMatXN mData3D = RMatXN(0, N);
};

#endif /* Property_hpp */
