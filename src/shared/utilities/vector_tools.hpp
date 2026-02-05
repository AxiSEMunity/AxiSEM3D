//
//  vector_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 3/5/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  std::vector tools

#ifndef vector_tools_hpp
#define vector_tools_hpp

#include <numeric>
#include <vector>
#include <map>
#include <algorithm>
#include "numerical.hpp"

namespace vector_tools {
    // sorted interpolation
    template <typename T>
    void linearInterpSorted(const std::vector<T> &points, T value,
                            int &index0, int &index1,
                            T &factor0, T &factor1,
                            int begin = 0, int end = -1) {
        // truncate points to specified range
        if (end < 0) {
            end = (int)points.size();
        }
        auto itBegin = points.begin() + begin;
        auto itEnd = points.begin() + end;
        
        // check size
        if (itEnd - itBegin < 2) {
            throw std::runtime_error("vector_tools::linearInterpSorted || "
                                     "At least two points must be provided.");
        }
        
        // check value
        T tolerance = ((*(itEnd - 1) - *(itBegin)) * numerical::epsilon<T>());
        if (value < *(itBegin) - tolerance ||
            value > *(itEnd - 1) + tolerance) {
            throw std::runtime_error("vector_tools::linearInterpSorted || "
                                     "Value is out of range.");
        }
        
        // find
        index1 = (int)(std::upper_bound(itBegin, itEnd, value) - itBegin);
        if (index1 == 0) {
            // value slightly smaller than begin
            index1++;
        } else if (index1 == itEnd - itBegin) {
            // value slightly bigger than end
            index1--;
        }
        index0 = index1 - 1;
        
        // shift to absolute location
        index0 += begin;
        index1 += begin;
        
        // compute value using at() to check bounds
        factor1 = ((value - points.at(index0)) /
                   (points.at(index1) - points.at(index0)));
        factor0 = (T)1. - factor1;
    }
    
    // check sorted unique
    template <typename T>
    bool isSortedUnique(const std::vector<T> &vec) {
        auto adjIter = std::adjacent_find(vec.begin(), vec.end(),
                                          std::greater_equal<T>());
        return adjIter == vec.end();
    }
    
    // sort and unique std::vector
    template <typename T>
    void sortUnique(std::vector<T> &vec) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }
    
    // find in a sorted unique std::vector
    template <typename T>
    bool findSortedUnique(const std::vector<T> &vec, const T &val) {
        return std::binary_search(vec.begin(), vec.end(), val);
    }
    
    // sum of std::vector
    template <typename T>
    T sum(const std::vector<T> &vec) {
        return std::accumulate(vec.begin(), vec.end(), (T)0);
    }
    
    // total size in nested std::vector
    template <typename T>
    size_t totalSize(const std::vector<T> &vec) {
        size_t total = 0;
        for (const T &v: vec) {
            total += v.size();
        }
        return total;
    }
    
    // max length in std::vector<std::string>
    template <typename Vec>
    int maxLength(const Vec &vec) {
        int mkl = 0;
        for (auto it = vec.begin(); it != vec.end(); ++it) {
            mkl = std::max(mkl, (int)(it->length()));
        }
        return mkl;
    }
    
    // max key length in std::map
    template <typename T>
    int maxKeyLength(const std::map<std::string, T> &map) {
        int mkl = 0;
        for (auto it = map.begin(); it != map.end(); ++it) {
            mkl = std::max(mkl, (int)(it->first.length()));
        }
        return mkl;
    }
    
    // sum of values in std::map
    template <typename T>
    T sumValues(const std::map<std::string, T> &map) {
        T sum = (T)0;
        for (auto it = map.begin(); it != map.end(); ++it) {
            sum += it->second;
        }
        return sum;
    }
    
    // aggregate value by key in std::map
    template <typename T>
    void aggregate(std::map<std::string, T> &map,
                   const std::string &key, T value) {
        map.insert({key, (T)0});
        map.at(key) += value;
    }
}

#endif /* vector_tools_hpp */
