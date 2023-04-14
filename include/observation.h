#ifndef INEKF_OBSERVATION_H
#define INEKF_OBSERVATION_H

#include "common.h"

namespace inekf {

    // class Observation: general observation models
    class Observation {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        Observation() = default;

        Observation(Eigen::VectorXd& Y, Eigen::VectorXd& b, Eigen::MatrixXd& H, Eigen::MatrixXd& N, Eigen::MatrixXd& Pi);

        // destructor
        ~Observation() = default;

    private:
        // measurement vector Y
        Eigen::VectorXd Y_;
        // constant vector b
        Eigen::VectorXd b_;
        // observation model Jacobian matrix H
        Eigen::MatrixXd H_;
        // observation model covariance
        Eigen::MatrixXd N_;
        // auxiliary selection matrix \Pi
        Eigen::MatrixXd Pi_;
    };

    // class Kinematics: kinematics measurement with index, pose and covariance
    class Kinematics {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        Kinematics() = default;

        Kinematics(int id_in, Eigen::Matrix4d pose_in, Eigen::Matrix<double, 6, 6> covariance_in);

        // destructor
        ~Kinematics() = default;

    public:
        // timestamp
        double timestamp = 0;
        // measurement index
        int id = 0;
        // measurement pose
        Eigen::Matrix4d pose;
        // measurement covariance
        Eigen::Matrix<double, 6, 6> covariance;
    };

    // class Landmark: landmark measurement with index, position and covariance
    class Landmark {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        Landmark() = default;

        Landmark(int id_in, Eigen::Vector3d position_in, Eigen::Matrix3d covariance_in);

        // destructor
        ~Landmark() = default;

    public:
        // timestamp
        double timestamp = 0;
        // measurement index
        int id = 0;
        // measurement position
        Eigen::Vector3d position;
        // measurement covariance
        Eigen::Matrix3d covariance;
    };

    // types for measurements and measurement vectors
    // a map with an integer as key and a 3-dimensional vector as value
    typedef std::map<unsigned int, Eigen::Vector3d, std::less<unsigned int>, Eigen::aligned_allocator<std::pair<const unsigned int, Eigen::Vector3d>>> mapIntVector3d;
    typedef std::map<unsigned int, Eigen::Vector3d, std::less<unsigned int>, Eigen::aligned_allocator<std::pair<const unsigned int, Eigen::Vector3d>>>::const_iterator mapIntVector3dIterator;
    // a vector of kinematics measurements
    typedef std::vector<Kinematics, Eigen::aligned_allocator<Kinematics>> vectorKinematics;
    typedef std::vector<Kinematics, Eigen::aligned_allocator<Kinematics>>::const_iterator vectorKinematicsIterator;
    // a vector of landmark measurements
    typedef std::vector<Landmark, Eigen::aligned_allocator<Landmark>> vectorLandmarks;
    typedef std::vector<Landmark, Eigen::aligned_allocator<Landmark>>::const_iterator vectorLandmarksIterator;

}

#endif //INEKF_OBSERVATION_H
