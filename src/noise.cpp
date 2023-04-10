#include "noise.h"

namespace inekf {

    // class Noise
    Noise::Noise() {
        // set default values
        setGyroscopeNoise(0.01);
        setAccelerometerNoise(0.1);
        setGyroscopeBiasNoise(0.00001);
        setAccelerometerBiasNoise(0.0001);
        setContactNoise(0.1);
    }

    void Noise::setGyroscopeNoise(double std) { Q_g_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setGyroscopeNoise(const Eigen::Vector3d &std_vec) {
        Q_g_ << std_vec(0) * std_vec(0) , 0, 0,
                0, std_vec(1) * std_vec(1), 0,
                0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setGyroscopeNoise(const Eigen::Matrix3d &cov) { Q_g_ = cov; }

    void Noise::setAccelerometerNoise(double std) { Q_a_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setAccelerometerNoise(const Eigen::Vector3d &std_vec) {
        Q_a_ << std_vec(0) * std_vec(0) , 0, 0,
                0, std_vec(1) * std_vec(1), 0,
                0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setAccelerometerNoise(const Eigen::Matrix3d &cov) { Q_a_ = cov; }

    void Noise::setGyroscopeBiasNoise(double std) { Q_bg_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setGyroscopeBiasNoise(const Eigen::Vector3d &std_vec) {
        Q_bg_ << std_vec(0) * std_vec(0) , 0, 0,
                 0, std_vec(1) * std_vec(1), 0,
                 0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setGyroscopeBiasNoise(const Eigen::Matrix3d &cov) { Q_bg_ = cov; }

    void Noise::setAccelerometerBiasNoise(double std) { Q_ba_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setAccelerometerBiasNoise(const Eigen::Vector3d &std_vec) {
        Q_ba_ << std_vec(0) * std_vec(0) , 0, 0,
                 0, std_vec(1) * std_vec(1), 0,
                 0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setAccelerometerBiasNoise(const Eigen::Matrix3d &cov) { Q_ba_ = cov; }

    void Noise::setContactNoise(double std) { Q_c_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setContactNoise(const Eigen::Vector3d& std_vec) {
        Q_c_ << std_vec(0) * std_vec(0) , 0, 0,
                0, std_vec(1) * std_vec(1), 0,
                0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setContactNoise(const Eigen::Matrix3d &cov) { Q_c_ = cov; }

    void Noise::setLandmarkNoise(double std) { Q_l_ = std * std * Eigen::Matrix3d::Identity(); }

    void Noise::setLandmarkNoise(const Eigen::Vector3d &std_vec) {
        Q_l_ << std_vec(0) * std_vec(0) , 0, 0,
                0, std_vec(1) * std_vec(1), 0,
                0, 0, std_vec(2) * std_vec(2);
    }

    void Noise::setLandmarkNoise(const Eigen::Matrix3d &cov) { Q_l_ = cov; }

    Eigen::Matrix3d Noise::getGyroscopeCov() { return Q_g_; }

    Eigen::Matrix3d Noise::getAccelerometerCov() { return Q_a_; }

    Eigen::Matrix3d Noise::getGyroscopeBiasCov() { return Q_bg_; }

    Eigen::Matrix3d Noise::getAccelerometerBiasCov() { return Q_ba_; }

    Eigen::Matrix3d Noise::getContactCov() { return Q_c_; }

    Eigen::Matrix3d Noise::getLandmarkCov() { return Q_l_; }

    std::ostream &operator<<(std::ostream &os, const Noise &p) {
        os << "---------- Noise ---------- \n" << std::endl;
        os << "Gyroscope Measurement Covariance: \n" << p.Q_g_ << std::endl << std::endl;
        os << "Accelerometer Measurement Covariance: \n" << p.Q_a_ << std::endl << std::endl;
        os << "Gyroscope Bias Covariance: \n" << p.Q_bg_ << std::endl << std::endl;
        os << "Accelerometer Bias Covariance: \n" << p.Q_ba_ << std::endl << std::endl;
        os << "Contact Measurement Covariance: \n" << p.Q_c_ << std::endl << std::endl;
        os << "Landmark Measurement Covariance: \n" << p.Q_l_ << std::endl << std::endl;
        os << "--------------------------- \n" << std::endl;

        return os;
    }

}
