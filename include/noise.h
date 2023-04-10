#ifndef INEKF_NOISE_H
#define INEKF_NOISE_H

#include "common.h"

namespace inekf {

    // class Noise: noise covariance matrices for angular velocity & translational velocity measurements, gyroscope & accelerometer biases, contact & landmark measurements
    class Noise {

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructor
        Noise();

        // destructor
        ~Noise() = default;

        // functions
        // set parameters
        // gyroscope measurement
        void setGyroscopeNoise(double std);

        void setGyroscopeNoise(const Eigen::Vector3d& std_vec);

        void setGyroscopeNoise(const Eigen::Matrix3d& cov);
        // accelerometer measurement
        void setAccelerometerNoise(double std);

        void setAccelerometerNoise(const Eigen::Vector3d& std_vec);

        void setAccelerometerNoise(const Eigen::Matrix3d& cov);
        // gyroscope bias
        void setGyroscopeBiasNoise(double std);

        void setGyroscopeBiasNoise(const Eigen::Vector3d& std_vec);

        void setGyroscopeBiasNoise(const Eigen::Matrix3d& cov);
        // accelerometer bias
        void setAccelerometerBiasNoise(double std);

        void setAccelerometerBiasNoise(const Eigen::Vector3d& std_vec);

        void setAccelerometerBiasNoise(const Eigen::Matrix3d& cov);
        // contact measurement
        void setContactNoise(double std);

        void setContactNoise(const Eigen::Vector3d& std_vec);

        void setContactNoise(const Eigen::Matrix3d& cov);
        // landmark measurement
        void setLandmarkNoise(double std);

        void setLandmarkNoise(const Eigen::Vector3d& std_vec);

        void setLandmarkNoise(const Eigen::Matrix3d& cov);

        // get parameters
        Eigen::Matrix3d getGyroscopeCov();

        Eigen::Matrix3d getAccelerometerCov();

        Eigen::Matrix3d getGyroscopeBiasCov();

        Eigen::Matrix3d getAccelerometerBiasCov();

        Eigen::Matrix3d getContactCov();

        Eigen::Matrix3d getLandmarkCov();

        // output stream
        friend std::ostream& operator<<(std::ostream& os, const Noise& p);

    private:
        // covariance for angular velocity measurements
        Eigen::Matrix3d Q_g_;
        // covariance for linear acceleration measurements
        Eigen::Matrix3d Q_a_;

        // covariance for gyroscope bias
        Eigen::Matrix3d Q_bg_;
        // covariance for accelerometer bias
        Eigen::Matrix3d Q_ba_;

        // covariance for contact measurements
        Eigen::Matrix3d Q_c_;
        // covariance for landmark measurements
        Eigen::Matrix3d Q_l_;
    };

}

#endif //INEKF_NOISE_H
