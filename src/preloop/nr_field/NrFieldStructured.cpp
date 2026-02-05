//
//  NrFieldStructured.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  structured Nr(s,z)

#include "NrFieldStructured.hpp"
#include "geodesy.hpp"

// constructor
NrFieldStructured::NrFieldStructured(const std::string &fname,
                                     int valOutOfRange):
mFilename(fname), mValueOutOfRange(valOutOfRange) {
    // coords
    std::array<std::string, 2> coordNames;
    if (geodesy::isCartesian()) {
        coordNames = {"structured_s", "structured_z"};
    } else {
        coordNames = {"structured_r", "structured_t"};
    }
    // data
    std::vector<std::pair<std::string, double>> varInfo = {
        std::make_pair("structured_Nr", 1.)};
    mGrid = std::make_unique<const StructuredGrid<2, int>>
    (fname, coordNames, varInfo, std::array<int, 2>({0, 1}));
}

// get nr by (s, z)
eigen::IColX NrFieldStructured::
getNrAtPoints(const eigen::DMatX2_RM &sz) const {
    eigen::IColX nr(sz.rows());
    if (geodesy::isCartesian()) {
        // Cartesian
        for (int ip = 0; ip < sz.rows(); ip++) {
            nr(ip) = mGrid->compute(sz.row(ip), mValueOutOfRange);
        }
    } else {
        // spherical
        const auto &rt = geodesy::sz2rtheta(sz, true);
        for (int ip = 0; ip < sz.rows(); ip++) {
            nr(ip) = mGrid->compute(rt.row(ip), mValueOutOfRange);
        }
    }
    return nr;
}

// verbose
std::string NrFieldStructured::verbose() const {
    using namespace bstring;
    
    // title
    std::stringstream ss;
    ss << boxTitle("Nr(s,z)");
    ss << boxEquals(0, 18, "type", "STRUCTURED");
    ss << boxEquals(0, 18, "NetCDF file", mFilename);
    ss << boxEquals(0, 18, "out-of-range value", mValueOutOfRange);
    
    // range
    const auto &crds = mGrid->getGridCoords();
    const std::string &rangeX1 = range(crds[0].front(), crds[0].back());
    const std::string &rangeX2 = range(crds[1].front(), crds[1].back());
    if (geodesy::isCartesian()) {
        ss << boxEquals(0, 18, "range of s-axis", rangeX1);
        ss << boxEquals(0, 18, "range of z-axis", rangeX2);
    } else {
        ss << boxEquals(0, 18, "range of r-axis", rangeX1);
        ss << boxEquals(0, 19, "range of θ-axis", rangeX2);
    }
    
    // value
    const auto &data = mGrid->getGridData();
    ss << boxEquals(0, 18, "min Nr", data.minimum());
    ss << boxEquals(0, 18, "max Nr", data.maximum());
    long sum = ((Eigen::Tensor<long, 0, Eigen::RowMajor>)
                data.cast<long>().sum())(0);
    int mean = (int)round(1. * sum / data.size());
    ss << boxEquals(0, 18, "ave Nr", mean);
    ss << boxEquals(0, 18, "sum Nr", sum);
    ss << boxBaseline() << "\n\n";
    return ss.str();
}
