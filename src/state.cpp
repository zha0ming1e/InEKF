#include "state.h"
#include "matrix_lie_group.h"

namespace inekf {

    // class BaseState
    BaseState::BaseState() : X_(2), P_(Eigen::MatrixXd::Identity(9, 9)) {
        //
    }

    BaseState::BaseState(unsigned int K) : X_(K) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(const Eigen::MatrixXd& X) : X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(const SEK3& X) : X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(const Eigen::MatrixXd& X, const Eigen::MatrixXd& P) : X_(X), P_(P) {
        //
    }

    BaseState::BaseState(const SEK3& X, const Eigen::MatrixXd& P) : X_(X), P_(P) {
        //
    }

    BaseState::BaseState(double t, unsigned int K) : t_(t), X_(K) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(double t, const Eigen::MatrixXd& X) : t_(t), X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(double t, const SEK3& X) : t_(t), X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    BaseState::BaseState(double t, const Eigen::MatrixXd& X, const Eigen::MatrixXd& P) : t_(t), X_(X), P_(P) {
        //
    }

    BaseState::BaseState(double t, const SEK3& X, const Eigen::MatrixXd& P) : t_(t), X_(X), P_(P) {
        //
    }

    SEK3 BaseState::getWorldX() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getX();
        } else {
            return this->getXInverse();
        }
    }

    Eigen::MatrixXd BaseState::getWorldXMatrix() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getXMatrix();
        } else {
            return this->getXInverseMatrix();
        }
    }

    SEK3 BaseState::getBodyX() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getX();
        } else {
            return this->getXInverse();
        }
    }

    Eigen::MatrixXd BaseState::getBodyXMatrix() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getXMatrix();
        } else {
            return this->getXInverseMatrix();
        }
    }

    void BaseState::copyDiagXMatrix(int n, Eigen::MatrixXd& BigX) const {
        unsigned int dimX = this->getXDim();
        auto Xmat = X_.getX();
        for(int i = 0; i < n; ++i) {
            unsigned int startIndex = BigX.rows();
            BigX.conservativeResize(startIndex + dimX, startIndex + dimX);
            BigX.block(startIndex, 0, dimX, startIndex) = Eigen::MatrixXd::Zero(dimX, startIndex);
            BigX.block(0, startIndex, startIndex, dimX) = Eigen::MatrixXd::Zero(startIndex, dimX);
            BigX.block(startIndex, startIndex, dimX, dimX) = Xmat;
        }
    }

    void BaseState::copyDiagXInverseMatrix(int n, Eigen::MatrixXd& BigXinv) const {
        unsigned int dimX = this->getXDim();
        auto Xinv = this->getXInverseMatrix();
        for(int i = 0; i < n; ++i) {
            unsigned int startIndex = BigXinv.rows();
            BigXinv.conservativeResize(startIndex + dimX, startIndex + dimX);
            BigXinv.block(startIndex, 0, dimX, startIndex) = Eigen::MatrixXd::Zero(dimX, startIndex);
            BigXinv.block(0, startIndex, startIndex, dimX) = Eigen::MatrixXd::Zero(startIndex, dimX);
            BigXinv.block(startIndex, startIndex, dimX, dimX) = Xinv;
        }
    }

    // class State
    State::State() : BaseState(2), Theta_(Eigen::VectorXd::Zero(6)) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(const Eigen::MatrixXd& X) : BaseState(X), Theta_(Eigen::VectorXd::Zero(6)) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(const SEK3& X) : BaseState(X), Theta_(Eigen::VectorXd::Zero(6)) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta) : BaseState(X), Theta_(Theta) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(const SEK3& X, const Eigen::VectorXd& Theta) : BaseState(X), Theta_(Theta) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(X, P), Theta_(Theta) {
        //
    }

    State::State(const SEK3& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(X, P), Theta_(Theta) {
        //
    }

    State::State(double t, const Eigen::MatrixXd& X) : BaseState(t, X), Theta_(Eigen::VectorXd::Zero(6)) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(double t, const SEK3& X) : BaseState(t, X), Theta_(Eigen::VectorXd::Zero(6)) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(double t, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta) : BaseState(t, X), Theta_(Theta) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(double t, const SEK3& X, const Eigen::VectorXd& Theta) : BaseState(t, X), Theta_(Theta) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    State::State(double t, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(t, X, P), Theta_(Theta) {
        //
    }

    State::State(double t, const SEK3& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(t, X, P), Theta_(Theta) {
        //
    }

    Eigen::Matrix3d State::getWorldRotation() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getRotation();
        } else {
            return this->getRotation().transpose();
        }
    }

    Eigen::Vector3d State::getWorldVelocity() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getVelocity();
        } else {
            return -this->getRotation().transpose() * this->getVelocity();
        }
    }

    Eigen::Vector3d State::getWorldPosition() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getPosition();
        } else {
            return -this->getRotation().transpose() * this->getPosition();
        }
    }

    Eigen::Matrix3d State::getBodyRotation() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getRotation();
        } else {
            return this->getRotation().transpose();
        }
    }

    Eigen::Vector3d State::getBodyVelocity() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getVelocity();
        } else {
            return -this->getRotation().transpose() * this->getVelocity();
        }
    }

    Eigen::Vector3d State::getBodyPosition() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getPosition();
        } else {
            return -this->getRotation().transpose() * this->getPosition();
        }
    }

    std::ostream& operator<<(std::ostream& os, const State& s) {
        os << "---------- State ---------- \n" << std::endl;
        os << "t: \n" << s.getT() << std::endl << std::endl;
        os << "X: \n" << s.getXMatrix() << std::endl << std::endl;
        os << "Theta: \n" << s.getTheta() << std::endl << std::endl;
        os << "P: \n" << s.getP() << std::endl;
        os << "--------------------------- \n" << std::endl;
        return os;
    }

}
