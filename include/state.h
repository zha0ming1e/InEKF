#ifndef INEKF_STATE_H
#define INEKF_STATE_H

#include "matrix_lie_group.h"

namespace inekf {

    // class BaseState: base system state on matrix Lie group SE_K(3)
    class BaseState {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        BaseState();

        explicit BaseState(unsigned int K);

        explicit BaseState(const Eigen::MatrixXd& X);

        explicit BaseState(const SEK3& X);

        BaseState(const Eigen::MatrixXd& X, const Eigen::MatrixXd& P);

        BaseState(const SEK3& X, const Eigen::MatrixXd& P);

        BaseState(double timestamp, unsigned int K);

        BaseState(double timestamp, const Eigen::MatrixXd& X);

        BaseState(double timestamp, const SEK3& X);

        BaseState(double timestamp, const Eigen::MatrixXd& X, const Eigen::MatrixXd& P);

        BaseState(double timestamp, const SEK3& X, const Eigen::MatrixXd& P);

        // destructor
        ~BaseState() = default;

        // functions
        // get
        double getTimestamp() const { return timestamp_; }

        SEK3 getX() const { return X_; }

        Eigen::MatrixXd getXMatrix() const { return X_.getX(); }

        Eigen::MatrixXd getP() const { return P_; }

        Eigen::Vector3d getVector(unsigned int id) const { return X_.getVector(id); }
        // dimension of state matrix X: columns
        unsigned int getXDim() const { return X_.getXDim(); }
        // degree of freedom of state
        unsigned int getXDoF() const { return X_.getXDoF(); }
        // dimension of state covariance matrix P: columns
        unsigned int getPDim() const { return P_.cols(); }

        StateType getStateType() const { return state_type_; }

        SEK3 getXInverse() const { return X_.getInverse(); }

        Eigen::MatrixXd getXInverseMatrix() const { return X_.getInverse().getX(); }

        SEK3 getWorldX() const;

        Eigen::MatrixXd getWorldXMatrix() const;

        SEK3 getBodyX() const;

        Eigen::MatrixXd getBodyXMatrix() const;
        // set
        void setT(const double timestamp) { timestamp_ = timestamp; }

        void setX(const SEK3& X) { X_ = X; }

        void setX(const Eigen::MatrixXd& X) { X_ = SEK3(X); }

        void setP(const Eigen::MatrixXd& P) { P_ = P; }

        void copyDiagXMatrix(int n, Eigen::MatrixXd& BigX) const;

        void copyDiagXInverseMatrix(int n, Eigen::MatrixXd& BigXinv) const;

    protected:
        // timestamp
        double timestamp_ = 0;
        // state type: world-centric or robot body centric
        StateType state_type_ = StateType::WorldCentric;
        // R (rotation), v (velocity), p (position), ... states on the matrix Lie group SE_K(3)
        // R v p ...
        // 0 1 0 ...
        // 0 0 1 ...
        // . . . ...
        SEK3 X_;
        // state covariance P
        Eigen::MatrixXd P_;
    };

    // class State: system state class including IMU bias states on Euclidean vector space
    class State : public BaseState {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        State();

        explicit State(const Eigen::MatrixXd& X);

        explicit State(const SEK3& X);

        State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta);

        State(const SEK3& X, const Eigen::VectorXd& Theta);

        State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P);

        State(const SEK3& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P);

        State(double timestamp, const Eigen::MatrixXd& X);

        State(double timestamp, const SEK3& X);

        State(double timestamp, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta);

        State(double timestamp, const SEK3& X, const Eigen::VectorXd& Theta);

        State(double timestamp, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P);

        State(double timestamp, const SEK3& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P);

        // destructor
        ~State() = default;

        // functions
        Eigen::VectorXd getTheta() const { return Theta_; }

        Eigen::Matrix3d getRotation() const { return X_.getRotationMatrix(); }

        Eigen::Vector3d getVelocity() const { return X_.getVector(3); }

        Eigen::Vector3d getPosition() const { return X_.getVector(4); }

        Eigen::Vector3d getGyroscopeBias() const { return Theta_.head(3); }

        Eigen::Vector3d getAccelerometerBias() const { return Theta_.tail(3); }

        Eigen::Matrix3d getRotationCovariance() const { return P_.block<3, 3>(0, 0); }

        Eigen::Matrix3d getVelocityCovariance() const { return P_.block<3, 3>(3, 3); }

        Eigen::Matrix3d getPositionCovariance() const { return P_.block<3, 3>(6, 6); }

        Eigen::Matrix3d getGyroscopeBiasCovariance() const { return P_.block<3, 3>(P_.rows() - 6, P_.rows() - 6); }

        Eigen::Matrix3d getAccelerometerBiasCovariance() const { return P_.block<3, 3>(P_.rows() - 3, P_.rows() - 3); }
        // dimension of bias state vector \Theta: rows
        unsigned int getThetaDim() const { return Theta_.rows(); }

        Eigen::Matrix3d getWorldRotation() const;

        Eigen::Vector3d getWorldVelocity() const;

        Eigen::Vector3d getWorldPosition() const;

        Eigen::Matrix3d getBodyRotation() const;

        Eigen::Vector3d getBodyVelocity() const;

        Eigen::Vector3d getBodyPosition() const;

        void setTheta(const Eigen::VectorXd& Theta) { Theta_ = Theta; }

        void setRotation(const Eigen::Matrix3d& R) { X_.setRotation(R); }

        void setRotation(const SEK3& sek3) { X_.setRotation(sek3); }

        void setVelocity(const Eigen::Vector3d& v) { X_.setVector(3, v); }

        void setPosition(const Eigen::Vector3d& p) { X_.setVector(4, p); }

        void setGyroscopeBias(const Eigen::Vector3d& bg) { Theta_.head(3) = bg; }

        void setAccelerometerBias(const Eigen::Vector3d& ba) { Theta_.tail(3) = ba; }

        void setRotationCovariance(const Eigen::Matrix3d& cov) { P_.block<3, 3>(0,0) = cov; }

        void setVelocityCovariance(const Eigen::Matrix3d& cov) { P_.block<3, 3>(3,3) = cov; }

        void setPositionCovariance(const Eigen::Matrix3d& cov) { P_.block<3, 3>(6,6) = cov; }

        void setGyroscopeBiasCovariance(const Eigen::Matrix3d& cov) { P_.block<3, 3>(P_.rows() - 6,P_.rows() - 6) = cov; }

        void setAccelerometerBiasCovariance(const Eigen::Matrix3d& cov) { P_.block<3, 3>(P_.rows() - 3,P_.rows() - 3) = cov; }

        friend std::ostream& operator<<(std::ostream& os, const State& s);

    private:
        // b_g (gyroscope bias), b_a (accelerometer bias) states on the Euclidean vector space
        Eigen::VectorXd Theta_;
        // b_g (gyroscope bias), b_a (accelerometer bias) states on the translation matrix Lie group T(N)
        //TN Theta_TN_;
    };

}

#endif //INEKF_STATE_H
