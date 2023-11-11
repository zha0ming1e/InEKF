#include "state_template.h"
#include "matrix_lie_group_template.h"

namespace inekf {

    // class BaseState
    template <unsigned int K>
    BaseState<K>::BaseState() : X_(SEK3<K>()), P_(Eigen::MatrixXd::Identity(9, 9)) {
        //
    }

    template <unsigned int K>
    BaseState<K>::BaseState(const Eigen::MatrixXd& X) : X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    template <unsigned int K>
    BaseState<K>::BaseState(const SEK3<K>& X) : X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    template <unsigned int K>
    BaseState<K>::BaseState(const Eigen::MatrixXd& X, const Eigen::MatrixXd& P) : X_(X), P_(P) {
        //
    }

    template <unsigned int K>
    BaseState<K>::BaseState(const SEK3<K>& X, const Eigen::MatrixXd& P) : X_(X), P_(P) {
        //
    }

    template <unsigned int K>
    BaseState<K>::BaseState(double timestamp) : timestamp_(timestamp), X_(SEK3<K>()) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    template <unsigned int K>
    BaseState<K>::BaseState(double timestamp, const Eigen::MatrixXd& X) : timestamp_(timestamp), X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    template <unsigned int K>
    BaseState<K>::BaseState(double timestamp, const SEK3<K>& X) : timestamp_(timestamp), X_(X) {
        P_ = Eigen::MatrixXd::Identity(this->getXDoF(), this->getXDoF());
    }

    template <unsigned int K>
    BaseState<K>::BaseState(double timestamp, const Eigen::MatrixXd& X, const Eigen::MatrixXd& P) : timestamp_(timestamp), X_(X), P_(P) {
        //
    }

    template <unsigned int K>
    BaseState<K>::BaseState(double timestamp, const SEK3<K>& X, const Eigen::MatrixXd& P) : timestamp_(timestamp), X_(X), P_(P) {
        //
    }

    template <unsigned int K>
    SEK3<K> BaseState<K>::getWorldX() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getX();
        } else {
            return this->getXInverse();
        }
    }

    template <unsigned int K>
    Eigen::MatrixXd BaseState<K>::getWorldXMatrix() const {
        if (state_type_ == StateType::WorldCentric) {
            return this->getXMatrix();
        } else {
            return this->getXInverseMatrix();
        }
    }

    template <unsigned int K>
    SEK3<K> BaseState<K>::getBodyX() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getX();
        } else {
            return this->getXInverse();
        }
    }

    template <unsigned int K>
    Eigen::MatrixXd BaseState<K>::getBodyXMatrix() const {
        if (state_type_ == StateType::BodyCentric) {
            return this->getXMatrix();
        } else {
            return this->getXInverseMatrix();
        }
    }

    template <unsigned int K>
    void BaseState<K>::copyDiagXMatrix(int n, Eigen::MatrixXd& BigX) const {
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

    template <unsigned int K>
    void BaseState<K>::copyDiagXInverseMatrix(int n, Eigen::MatrixXd& BigXinv) const {
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
    template <unsigned int K>
    State<K>::State() : BaseState<K>(), Theta_(Eigen::VectorXd::Zero(6)) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(const Eigen::MatrixXd& X) : BaseState<K>(X), Theta_(Eigen::VectorXd::Zero(6)) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(const SEK3<K>& X) : BaseState<K>(X), Theta_(Eigen::VectorXd::Zero(6)) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta) : BaseState<K>(X), Theta_(Theta) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(const SEK3<K>& X, const Eigen::VectorXd& Theta) : BaseState<K>(X), Theta_(Theta) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState<K>(X, P), Theta_(Theta) {
        //
    }

    template <unsigned int K>
    State<K>::State(const SEK3<K>& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState<K>(X, P), Theta_(Theta) {
        //
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const Eigen::MatrixXd& X) : BaseState<K>(timestamp, X), Theta_(Eigen::VectorXd::Zero(6)) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const SEK3<K>& X) : BaseState<K>(timestamp, X), Theta_(Eigen::VectorXd::Zero(6)) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta) : BaseState<K>(timestamp, X), Theta_(Theta) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const SEK3<K>& X, const Eigen::VectorXd& Theta) : BaseState<K>(timestamp, X), Theta_(Theta) {
        this->P_ = Eigen::MatrixXd::Identity(this->getXDoF() + this->getThetaDim(), this->getXDoF() + this->getThetaDim());
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const Eigen::MatrixXd& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState<K>(timestamp, X, P), Theta_(Theta) {
        //
    }

    template <unsigned int K>
    State<K>::State(double timestamp, const SEK3<K>& X, const Eigen::VectorXd& Theta, const Eigen::MatrixXd& P) : BaseState<K>(timestamp, X, P), Theta_(Theta) {
        //
    }

    template <unsigned int K>
    Eigen::Matrix3d State<K>::getWorldRotation() const {
        if (this->state_type_ == StateType::WorldCentric) {
            return this->getRotation();
        } else {
            return this->getRotation().transpose();
        }
    }

    template <unsigned int K>
    Eigen::Vector3d State<K>::getWorldVelocity() const {
        if (this->state_type_ == StateType::WorldCentric) {
            return this->getVelocity();
        } else {
            return -this->getRotation().transpose() * this->getVelocity();
        }
    }

    template <unsigned int K>
    Eigen::Vector3d State<K>::getWorldPosition() const {
        if (this->state_type_ == StateType::WorldCentric) {
            return this->getPosition();
        } else {
            return -this->getRotation().transpose() * this->getPosition();
        }
    }

    template <unsigned int K>
    Eigen::Matrix3d State<K>::getBodyRotation() const {
        if (this->state_type_ == StateType::BodyCentric) {
            return this->getRotation();
        } else {
            return this->getRotation().transpose();
        }
    }

    template <unsigned int K>
    Eigen::Vector3d State<K>::getBodyVelocity() const {
        if (this->state_type_ == StateType::BodyCentric) {
            return this->getVelocity();
        } else {
            return -this->getRotation().transpose() * this->getVelocity();
        }
    }

    template <unsigned int K>
    Eigen::Vector3d State<K>::getBodyPosition() const {
        if (this->state_type_ == StateType::BodyCentric) {
            return this->getPosition();
        } else {
            return -this->getRotation().transpose() * this->getPosition();
        }
    }

    template <unsigned int K>
    std::ostream& operator<<(std::ostream& os, const State<K>& s) {
        os << "---------- State ---------- \n" << std::endl;
        os << "Timestamp: \n" << s.getTimestamp() << std::endl << std::endl;
        os << "X: \n" << s.getXMatrix() << std::endl << std::endl;
        os << "Theta: \n" << s.getTheta() << std::endl << std::endl;
        os << "P: \n" << s.getP() << std::endl;
        os << "--------------------------- \n" << std::endl;
        return os;
    }

}
