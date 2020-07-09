//
//  NrFieldPointwise.cpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/15/20.
//  Copyright © 2020 Kuangdai Leng. All rights reserved.
//

//  pointwise Nr(s,z)

#include "NrFieldPointwise.hpp"
#include "NetCDF_Reader.hpp"

// constructor
NrFieldPointwise::NrFieldPointwise(const std::string &fname, double factor,
                                   double distTolExact):
mFilename(fname), mFactor(factor), mDistTolExact(distTolExact) {
    // rtree
    std::array<std::pair<std::string, double>, 1> varInfo = {
        std::make_pair("pointwise_Nr", 1.)};
    mRTree = std::make_unique<RTreeND<2, 1, int>>
    (fname, "pointwise_sz", varInfo);
    
    // read starting Nr used by wavefield scanning
    NetCDF_Reader reader;
    reader.open(io::popInputDir(fname));
    try {
        eigen::IColX nrStart;
        reader.readMatrix("starting_Nr_for_scanning", nrStart);
        mSumNrStart = nrStart.sum();
    } catch (...) {
        // not created by wavefield scanning
        mSumNrStart = -1;
    }
}

// get nr by (s, z)
eigen::IColX NrFieldPointwise::getNrAtPoints(const eigen::DMatX2_RM &sz) const {
    eigen::IColX result(sz.rows());
    for (int ip = 0; ip < sz.rows(); ip++) {
        // get from r-tree without distance limit
        double max = std::numeric_limits<double>::max();
        int nr = mRTree->compute(sz.row(ip), 4, max, 0, mDistTolExact);
        // apply multiplication factor
        result(ip) = (int)round(nr * mFactor);
    }
    return result;
}

// verbose
std::string NrFieldPointwise::verbose() const {
    std::stringstream ss;
    // statistics for verbose
    const eigen::IColX &controlNr = mRTree->getAllValues();
    long sumNr = controlNr.cast<long>().sum();
    int aveNr = (int)round(1. * sumNr / controlNr.size());
    
    // verbose
    ss << bstring::boxTitle("Nr(s,z)");
    ss << bstring::boxEquals(0, 21, "type", "POINTWISE");
    ss << bstring::boxEquals(0, 21, "NetCDF file", mFilename);
    ss << bstring::boxEquals(0, 21, "multiplication factor", mFactor);
    ss << bstring::boxEquals(0, 21, "# control points", mRTree->size());
    ss << bstring::boxEquals(0, 21, "min Nr", controlNr.minCoeff());
    ss << bstring::boxEquals(0, 21, "max Nr", controlNr.maxCoeff());
    ss << bstring::boxEquals(0, 21, "ave Nr", aveNr);
    ss << bstring::boxEquals(0, 21, "sum Nr", sumNr);
    if (mSumNrStart != -1) {
        ss << "* This Nr(s,z) was created by wavefield scanning:\n";
        ss << bstring::boxEquals(2, 19, "compression ratio",
                                 1. * sumNr / mSumNrStart);
    }
    ss << bstring::boxBaseline() << "\n\n";
    return ss.str();
}
