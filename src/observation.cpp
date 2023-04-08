//
// Created by zhaomingle on 4/5/23.
//

#include "observation.h"

namespace inekf {

    // class Observation
    Observation::Observation(Eigen::VectorXd& Y, Eigen::VectorXd& b,
                             Eigen::MatrixXd& H, Eigen::MatrixXd& N,
                             Eigen::MatrixXd& Pi) : Y_(Y), b_(b), H_(H), N_(N), Pi_(Pi) {
        //
    }

    // class Kinematics
    Kinematics::Kinematics(int id_in, Eigen::Matrix4d pose_in,
                           Eigen::Matrix<double,6,6> covariance_in)
                           : id(id_in), pose(pose_in), covariance(covariance_in) {
        //
    }

    // class Landmark
    Landmark::Landmark(int id_in, Eigen::Vector3d position_in,
                       Eigen::Matrix3d covariance_in)
                       : id(id_in), position(position_in), covariance(covariance_in) {
        //
    }

}
