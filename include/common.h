#ifndef INEKF_COMMON_H
#define INEKF_COMMON_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <mutex>

//#include <omp.h>

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
        // output type: measurement output of the system
        OUTPUT
    };

    // remove num_to_remove rows and columns starting from a given index idx in a square matrix
    inline void removeMatrixRowAndColumn(unsigned int idx, Eigen::MatrixXd& mat, unsigned int num_to_remove = 1) {
        unsigned int dim = mat.cols();
        // check reasonable parameters
        assert((idx >= 0) && (idx < dim) && (num_to_remove > 0) && (num_to_remove < (dim - idx + 1)));

        mat.block(idx, 0, dim - idx - num_to_remove, dim) = mat.bottomRows(dim - idx - num_to_remove).eval();
        mat.block(0, idx, dim, dim - idx - num_to_remove) = mat.rightCols(dim - idx - num_to_remove).eval();
        mat.conservativeResize(dim - num_to_remove, dim - num_to_remove);
    }

}

#endif //INEKF_COMMON_H
