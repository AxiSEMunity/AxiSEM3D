//
//  c1_tools.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 16/10/2024.
//  Copyright © 2024 Kuangdai Leng. All rights reserved.
//

#ifndef c1_tools_h
#define c1_tools_h

#include <algorithm>
#include <cmath>

namespace c1_tools {
  inline int
  columnSurf(bool includeSediment, bool includeIce) {
    int colSurf = 5; // no ice, no sediment
    if (includeIce) {
      colSurf = 1; // ice
    } else if (includeSediment) {
      colSurf = 2; // sediment
    }
    return colSurf;
  };

  inline void
  interpLinear(double target, const eigen::DColX& bases, int& loc, double& weight) {
    if (target < bases(0) || target > bases(bases.size() - 1)) {
      loc = -1;
      weight = 0.;
      return;
    }

    for (int i = 1; i < bases.size(); i++) {
      if (target <= bases(i)) {
        loc = i - 1;
        weight = 1. - 1. / (bases(loc + 1) - bases(loc)) * (target - bases(loc));
        return;
      }
    }
  };

  inline void
  gaussianSmoothing(eigen::DColX& data, int order, double dev, bool period) {
    if (data.size() == 0)
      return;
    order = std::min(order, ((int)data.size() + 1) / 2 - 1);
    if (order == 0)
      return;
    dev *= order;

    // gaussian kernel
    eigen::DColX gaussian(order * 2 + 1);
    for (int i = 0; i <= order; i++) {
      gaussian(order + i) = exp(-.5 * i * i / (dev * dev));
      gaussian(order - i) = gaussian(order + i);
    }
    gaussian /= gaussian.sum();

    // convolve
    eigen::DColX result = eigen::DColX::Zero(data.size());
    for (int i = 0; i < data.size(); i++) {
      for (int j = -order; j <= order; j++) {
        int k = i + j;
        if (period) {
          while (k < 0) {
            k += data.size();
          }
          while (k > data.size() - 1) {
            k -= data.size();
          }
        } else {
          // using fixed padding
          if (k < 0) {
            k = 0;
          }
          if (k > data.size() - 1) {
            k = (int)data.size() - 1;
          }
        }
        result(i) += gaussian(j + order) * data(k);
      }
    }

    // assign
    data = result;
  };

  inline void
  gaussianSmoothing(eigen::DMatXX& data,
      eigen::IColX orderRow,
      eigen::DColX devRow,
      bool periodRow,
      eigen::IColX orderCol,
      eigen::DColX devCol,
      bool periodCol) {
    for (int i = 0; i < data.rows(); i++) {
      eigen::DColX temp = data.row(i).transpose();
      gaussianSmoothing(temp, orderRow(i), devRow(i), periodRow);
      data.row(i) = temp.transpose();
    }
    for (int i = 0; i < data.cols(); i++) {
      eigen::DColX temp = data.col(i);
      gaussianSmoothing(temp, orderCol(i), devCol(i), periodCol);
      data.col(i) = temp;
    }
  };
} // namespace c1_tools

#endif /* c1_tools_h */
