//
//  RTreeND.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 9/7/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  boost RTree for searching nearest points

#ifndef RTreeND_hpp
#define RTreeND_hpp

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include <boost/geometry.hpp>
#pragma clang diagnostic pop

#include "NetCDF_Reader.hpp"
#include "mpi.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "eigen_tools.hpp"

// D: dimensions, must be 1 or 2 or 3
// V: number of variables
// T: type of variables
template <int D, int V, typename T>
class RTreeND {
    ///////////// typedef /////////////
    // coordinates
    typedef typename std::conditional<D == 1,
    Eigen::Matrix<double, Eigen::Dynamic, D>,
    Eigen::Matrix<double, Eigen::Dynamic, D, Eigen::RowMajor>>::type DMatXD_RM;
    // data
    typedef typename std::conditional<V == 1,
    Eigen::Matrix<T, Eigen::Dynamic, V>,
    Eigen::Matrix<T, Eigen::Dynamic, V, Eigen::RowMajor>>::type TMatXV_RM;
    // data unit
    typedef Eigen::Matrix<T, 1, V> TRowV;
    typedef Eigen::Matrix<double, 1, V> DRowV;
    // location
    typedef boost::geometry::model::point<double, D,
    boost::geometry::cs::cartesian> RTreeLoc;
    // leaf
    typedef std::pair<RTreeLoc, TRowV> RTreeLeaf;
    
    
    ///////////// cs array to RTreeLoc /////////////
    // 1D
    template <int R = D, typename ArrayCS>
    static typename std::enable_if<R == 1,
    RTreeLoc>::type toRTreeLoc(const ArrayCS &loc) {
        return RTreeLoc(loc[0]);
    }
    // 2D
    template <int R = D, typename ArrayCS>
    static typename std::enable_if<R == 2,
    RTreeLoc>::type toRTreeLoc(const ArrayCS &loc) {
        return RTreeLoc(loc[0], loc[1]);
    }
    // 3D
    template <int R = D, typename ArrayCS>
    static typename std::enable_if<R == 3,
    RTreeLoc>::type toRTreeLoc(const ArrayCS &loc) {
        return RTreeLoc(loc[0], loc[1], loc[2]);
    }
    
    
    ///////////// public /////////////
public:
    // constructor
    RTreeND() = default;
    
    // constructor
    RTreeND(const std::string &ncFile, const std::string &coordVarName,
            const std::array<std::pair<std::string, double>, V> &varInfo) {
        // read data
        DMatXD_RM ctrlCrds;
        TMatXV_RM ctrlVals;
        timer::gPreloopTimer.begin("Reading control-point data");
        timer::gPreloopTimer.message("data file: " + io::popInputDir(ncFile));
        if (mpi::root()) {
            // open
            NetCDF_Reader reader(io::popInputDir(ncFile));
            // read coords
            reader.readMatrixDouble(coordVarName, ctrlCrds);
            // read data
            ctrlVals.resize(ctrlCrds.rows(), V);
            for (int ivar = 0; ivar < V; ivar++) {
                // read to double
                eigen::DColX ctrlVal;
                reader.readMatrixDouble(varInfo[ivar].first, ctrlVal);
                ctrlVal *= varInfo[ivar].second;
                // round
                if (std::is_integral<T>::value) {
                    ctrlVal = ctrlVal.array().round();
                }
                // cast to Ts
                ctrlVals.col(ivar) = ctrlVal.template cast<T>();
            }
        }
        timer::gPreloopTimer.ended("Reading control-point data");
        
        // broadcast
        timer::gPreloopTimer.begin("Broadcasting control-point data");
        mpi::bcastEigen(ctrlCrds);
        mpi::bcastEigen(ctrlVals);
        
        // memory info
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(ctrlCrds, "coordinates at control points"));
        timer::gPreloopTimer.message
        (eigen_tools::memoryInfo(ctrlVals, "values at control points"));
        timer::gPreloopTimer.ended("Broadcasting control-point data");
        
        // build RTree
        timer::gPreloopTimer.begin("Building R-tree");
        addLeafs(ctrlCrds, ctrlVals);
        timer::gPreloopTimer.ended("Building R-tree");
    }
    
    // add a leaf
    template <typename ArrayCS>
    void addLeaf(const ArrayCS &loc, const TRowV &val) {
        mRTree.insert({toRTreeLoc(loc), val});
    }
    
    // add scalar, only for V = 1
    template <int VN = V, typename ArrayCS>
    typename std::enable_if<VN == 1, void>::type
    addLeaf(const ArrayCS &loc, T val) {
        static TRowV oneVal;
        oneVal(0) = val;
        addLeaf(loc, oneVal);
    }
    
    // add leafs
    void addLeafs(const DMatXD_RM &locs, const TMatXV_RM &vals) {
        for (int ir = 0; ir < locs.rows(); ir++) {
            mRTree.insert({toRTreeLoc(locs.row(ir)), vals.row(ir)});
        }
    }
    
    // query
    // difficult to make this collective
    template <typename ArrayCS>
    void query(const ArrayCS &loc, int count,
               std::vector<double> &dists, std::vector<TRowV> &vals) const {
        // location
        const RTreeLoc &rloc = toRTreeLoc(loc);
        // KNN query
        std::vector<RTreeLeaf> leafs;
        mRTree.query(boost::geometry::index::nearest(rloc, count),
                     std::back_inserter(leafs));
        // get distance
        dists.clear();
        dists.reserve(count);
        std::transform(leafs.begin(), leafs.end(), std::back_inserter(dists),
                       [&rloc](const auto &leaf) {
            return boost::geometry::distance(leaf.first, rloc);});
        // get values
        vals.clear();
        vals.reserve(count);
        std::transform(leafs.begin(), leafs.end(), std::back_inserter(vals),
                       [](const auto &leaf) {return leaf.second;});
    }
    
    // compute
    template <typename ArrayCS>
    TRowV compute(const ArrayCS &loc, int count, double maxDistInRange,
                  const TRowV &valOutOfRange, double distTolExact) const {
        // query
        static std::vector<double> dists;
        static std::vector<TRowV> vals;
        query(loc, count, dists, vals);
        
        // average vaules by inverse distance
        double invDistSum = 0.;
        static DRowV valTarget;
        valTarget.setZero();
        int numInRange = 0;
        for (int ip = 0; ip < dists.size(); ip++) {
            // exactly on a leaf
            if (dists[ip] < distTolExact) {
                return vals[ip];
            }
            // out of range
            if (dists[ip] > maxDistInRange) {
                continue;
            }
            // weight by inverse distance
            invDistSum += 1. / dists[ip];
            valTarget += vals[ip].template cast<double>() / dists[ip];
            numInRange++;
        }
        
        // out of range
        if (numInRange == 0) {
            return valOutOfRange;
        }
        
        // average and round
        valTarget /= invDistSum;
        if (std::is_integral<T>::value) {
            valTarget = valTarget.array().round();
        }
        return valTarget.template cast<T>();
    }
    
    // compute scalar, only for V = 1
    template <int VN = V, typename ArrayCS>
    typename std::enable_if<VN == 1, T>::type
    compute(const ArrayCS &loc, int count, double maxDistInRange,
            T valOutOfRange, double distTolExact) const {
        static TRowV oneVal;
        oneVal(0) = valOutOfRange;
        return compute(loc, count, maxDistInRange, oneVal, distTolExact)(0);
    }
    
    // size
    int size() const {
        return (int)mRTree.size();
    }
    
    // get all values
    TMatXV_RM getAllValues() const {
        // get values
        std::vector<TRowV> vals;
        std::transform(mRTree.begin(), mRTree.end(), std::back_inserter(vals),
                       [](const auto &leaf) {return leaf.second;});
        // cast to matrix
        TMatXV_RM mat(vals.size(), V);
        for (int ip = 0; ip < vals.size(); ip++) {
            mat.row(ip) = vals[ip];
        }
        return mat;
    }
    
private:
    // rtree
    boost::geometry::index::rtree<RTreeLeaf,
    boost::geometry::index::quadratic<16>> mRTree;
};

#endif /* RTreeND_hpp */
