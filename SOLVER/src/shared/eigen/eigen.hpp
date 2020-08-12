//
//  eigen.hpp
//  AxiSEM3D
//
//  Created by Kuangdai Leng on 1/7/19.
//  Copyright © 2019 Kuangdai Leng. All rights reserved.
//

//  Eigen3 header

#ifndef eigen_hpp
#define eigen_hpp

///////////////// preprocessors /////////////////
// no auto resizing
#define EIGEN_NO_AUTOMATIC_RESIZING

// no dynamic allocation once time loop begins
#define EIGEN_RUNTIME_NO_MALLOC

// disable static alignment to enhance code portability
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0

// initialize matrix with NAN
#define EIGEN_INITIALIZE_MATRICES_BY_NAN

///////////////// include /////////////////
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include <Eigen/Dense>
#pragma clang diagnostic pop

#endif /* eigen_hpp */
