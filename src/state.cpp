//
// Created by zhaomingle on 4/5/23.
//

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
        for(int i = 0; i < n; ++i) {
            unsigned int startIndex = BigX.rows();
            BigX.conservativeResize(startIndex + dimX, startIndex + dimX);
            BigX.block(startIndex, 0, dimX, startIndex) = Eigen::MatrixXd::Zero(dimX, startIndex);
            BigX.block(0, startIndex, startIndex, dimX) = Eigen::MatrixXd::Zero(startIndex, dimX);
            BigX.block(startIndex, startIndex, dimX, dimX) = X_.getX();
        }
    }

    void BaseState::copyDiagXInverseMatrix(int n, Eigen::MatrixXd& BigXinv) const {
        unsigned int dimX = this->getXDim();
        Eigen::MatrixXd Xinv = this->getXInverseMatrix();
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

    State::State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(X), Theta_(Theta) {
        P_ = P;
    }

    State::State(const SEK3& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState(X), Theta_(Theta) {
        P_ = P;
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
        os << "---------- State ----------" << std::endl;
        os << "X: \n" << s.getXMatrix() << std::endl << std::endl;
        os << "Theta: \n" << s.getTheta() << std::endl << std::endl;
        os << "P: \n" << s.getP() << std::endl;
        os << "-----------------------------------";
        return os;
    }

}
