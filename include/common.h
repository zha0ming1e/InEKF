//
// Created by zhaomingle on 4/4/23.
//

#ifndef INEKF_COMMON_H
#define INEKF_COMMON_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace inekf {

    // global variable of numerical tolerance for small rotation or other increment
    const double NUMERICAL_TOLERANCE = 1.0e-10;

    // state type: world centric and robot body centric
    enum class StateType {
        // world centric
        WorldCentric,
        // robot body centric
        BodyCentric
    };

    // error type: left invariant and right invariant
    enum class ErrorType {
        // left invariant
        LeftInvariant,
        // right invariant
        RightInvariant
    };

    // role type of IMU measurement: input and output
    enum class IMURoleType {
        // input type: driven input of the system
        INPUT,
        // output type: output measurement of the system
        OUTPUT
    };

    // remove a particular row and column with a same index
    inline void removeMatrixRowAndColumn(unsigned int index, Eigen::MatrixXd& mat) {
        unsigned int dim = mat.cols();
        mat.block(index, 0, dim - index - 1, dim) = mat.bottomRows(dim - index - 1).eval();
        mat.block(0, index, dim, dim - index - 1) = mat.rightCols(dim - index - 1).eval();
        mat.conservativeResize(dim - 1, dim - 1);
    }

}

#endif //INEKF_COMMON_H
