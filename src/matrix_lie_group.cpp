//
// Created by zhaomingle on 4/4/23.
//

#include "matrix_lie_group.h"

namespace inekf {

    long factorial(int n) {
        return (n == 0 || n == 1) ? 1 : (factorial(n - 1) * n);
    }

    // 3D Lie algebra
    Eigen::Matrix3d skew(const Eigen::Vector3d& phi) {
        Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
        M << 0, -phi[2], phi[1],
             phi[2], 0, -phi[0],
             -phi[1], phi[0], 0;

        return M;
    }

    Eigen::Matrix3d Gamma_m(const Eigen::Vector3d& phi, const int m) {
        // check m >= 0
        assert(m >= 0);

        // m-th order Gamma series: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (phi_^\wedge)^n
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d phix = skew(phi);
        double phi_norm = phi.norm();
        double phi_norm2 = phi_norm * phi_norm;

        // check if phi_norm_ is small enough
        if (phi_norm < NUMERICAL_TOLERANCE) {
            // ignore higher order n >= 2
            return ((1.0 / static_cast<double>(factorial(m)) * I) + (1.0 / static_cast<double>(factorial(1 + m)) * phix));
        }

        // closed form for Gamma_0, Gamma_1 and Gamma_2
        switch (m) {
            case 0: // Gamma_0: exponential map of SO(3)
                return (I + (sin(phi_norm) / phi_norm) * phix + ((1.0 - cos(phi_norm)) / phi_norm2) * phix * phix);
            case 1: // Gamma_1: left Jacobian of SO(3)
                return (I + ((1.0 - cos(phi_norm)) / phi_norm2) * phix + ((phi_norm - sin(phi_norm)) / (phi_norm * phi_norm2)) * phix * phix);
            case 2: // Gamma_2: Hartley et al. 2020
                return (0.5 * I + ((phi_norm - sin(phi_norm)) / (phi_norm * phi_norm2)) * phix + ((phi_norm2 + 2.0 * cos(phi_norm) - 2.0) / (2.0 * phi_norm2 * phi_norm2)) * phix * phix);
            default: { // general cases
                Eigen::Matrix3d R = I + (sin(phi_norm) / phi_norm) * phix + ((1.0 - cos(phi_norm)) / phi_norm2) * phix * phix;
                Eigen::Matrix3d S = I;
                Eigen::Matrix3d phixk = I;
                long k_factorial = 1;
                for (int k = 1; k <= m; ++k) {
                    k_factorial = k_factorial * k;
                    phixk = (phixk * phix).eval();
                    S = (S + (1.0 / static_cast<double>(k_factorial)) * phixk).eval();
                }
                // m is odd (&& m != 0) or even
                if (m % 2) {
                    // m is odd
                    return ((1.0 / static_cast<double>(k_factorial)) * I + (pow(-1.0, (m + 1.0) / 2.0) / pow(phi_norm, m + 1.0)) * phix * (R - S));
                } else {
                    // m is even
                    return ((1.0 / static_cast<double>(k_factorial)) * I + (pow(-1.0,m / 2.0) / pow(phi_norm, m)) * (R - S));
                }
            }
        }
    }

    // SO(3)
    Eigen::MatrixXd Hat_SO3(const Eigen::Vector3d& phi) {
        return skew(phi);
    }

    Eigen::Matrix3d Exp_SO3(const Eigen::Vector3d& phi) {
        return Gamma_m(phi, 0);
    }

    Eigen::Matrix3d LeftJacobian_SO3(const Eigen::Vector3d& phi) {
        return Gamma_m(phi, 1);
    }

    Eigen::Matrix3d RightJacobian_SO3(const Eigen::Vector3d& phi) {
        return Gamma_m(-phi, 1);
    }

    // SE_K(3)
    Eigen::MatrixXd Hat_SEK3(const Eigen::VectorXd& xi) {
        // the group of K direct isometries, K vectors
        int K = (static_cast<int>(xi.size()) - 3) / 3;
        // first 3D angular vector
        Eigen::Vector3d phi = xi.head(3);
        // skew of \phi
        Eigen::Matrix3d phix = skew(phi);

        // xi^\hat is the hat martrix of SE_K(3)
        Eigen::MatrixXd xi_hat = Eigen::MatrixXd::Zero(3 + K, 3 + K);
        // rotation skew block
        xi_hat.block<3, 3>(0, 0) = phix;
        // vector hat block
        for (int i = 0; i < K; ++i) {
            xi_hat.block<3, 1>(0, 3 + i) = xi.segment<3>(3 + 3 * i);
        }

        return xi_hat;
    }

    Eigen::MatrixXd Exp_SEK3(const Eigen::VectorXd& xi) {
        // the group of K direct isometries, K vectors
        int K = (static_cast<int>(xi.size()) - 3) / 3;
        // first 3D angular vector
        Eigen::Vector3d phi = xi.head(3);
        // class Gamma for \Gamma_m
        Gamma gamma = Gamma(phi);
        // rotation matrix R
        Eigen::Matrix3d R = gamma.getGamma(0);
        // left Jacobian J_l
        Eigen::Matrix3d Jl = gamma.getGamma(1);

        // matrix Lie group X
        Eigen::MatrixXd X = Eigen::MatrixXd::Identity(3 + K, 3 + K);
        // rotation matrix block
        X.block<3, 3>(0, 0) = R;
        // vector block
        for (int i = 0; i < K; ++i) {
            X.block<3, 1>(0, 3 + i) = Jl * xi.segment<3>(3 + 3 * i);
        }

        return X;
    }

    Eigen::MatrixXd Adjoint_SEK3(const Eigen::MatrixXd& X) {
        // the group of K direct isometries, K vectors
        int K = static_cast<int>(X.cols()) - 3;
        // adjoint matrix of X
        Eigen::MatrixXd Ad_X = Eigen::MatrixXd::Zero(3 + 3 * K, 3 + 3 * K);
        // rotation matrix
        Eigen::Matrix3d R = X.block<3, 3>(0, 0);

        // assembling
        Ad_X.block<3, 3>(0, 0) = R;
        for (int i = 0; i < K; ++i) {
            Ad_X.block<3, 3>(3 + 3 * i, 0) = skew(X.block<3, 1>(0, 3 + i)) * R;
            Ad_X.block<3, 3>(3 + 3 * i, 3 + 3 * i) = R;
        }

        return Ad_X;
    }

    // class Gamma
    Gamma::Gamma(Eigen::Vector3d phi) : phi_(phi) {
        phix_ = skew(phi_);
        phi_norm_ = phi_.norm();
        phi_norm2_ = phi_norm_ * phi_norm_;
        sin_phi_norm_ = sin(phi_norm_);
        cos_phi_norm_ = cos(phi_norm_);
    }

    Gamma::Gamma(double theta1, double theta2, double theta3) {
        phi_ = Eigen::Vector3d(theta1, theta2, theta3);
        phix_ = skew(phi_);
        phi_norm_ = phi_.norm();
        phi_norm2_ = phi_norm_ * phi_norm_;
        sin_phi_norm_ = sin(phi_norm_);
        cos_phi_norm_ = cos(phi_norm_);
    }

    Gamma::Gamma(double* theta) {
        phi_ = Eigen::Vector3d(theta[0], theta[1], theta[2]);
        phix_ = skew(phi_);
        phi_norm_ = phi_.norm();
        phi_norm2_ = phi_norm_ * phi_norm_;
        sin_phi_norm_ = sin(phi_norm_);
        cos_phi_norm_ = cos(phi_norm_);
    }

    Gamma::Gamma(const Gamma& gm) {
        phi_ = gm.getPhi();
        phix_ = gm.getPhix();
        phi_norm_ = gm.getPhiNorm();
        phi_norm2_ = gm.getPhiNorm2();
        sin_phi_norm_ = gm.getSinPhiNorm();
        cos_phi_norm_ = gm.getCosPhiNorm();
    }

    Eigen::Matrix3d Gamma::getGamma(int m) const {
        // check m >= 0
        assert(m >= 0);

        // m-th order Gamma series: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (phi_^\wedge)^n
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

        // check if phi_norm_ is small enough
        if (phi_norm_ < NUMERICAL_TOLERANCE) {
            // ignore higher order n >= 2
            return ((1.0 / static_cast<double>(factorial(m)) * I) + (1.0 / static_cast<double>(factorial(1 + m)) * phix_));
        }

        // closed form for \Gamma_0, \Gamma_1 and \Gamma_2
        switch (m) {
            case 0: // \Gamma_0: exponential map of SO(3)
                return (I + (sin_phi_norm_ / phi_norm_) * phix_ + ((1.0 - cos_phi_norm_) / phi_norm2_) * phix_ * phix_);
            case 1: // \Gamma_1: left Jacobian of SO(3)
                return (I + ((1.0 - cos_phi_norm_) / phi_norm2_) * phix_ + ((phi_norm_ - sin_phi_norm_) / (phi_norm_ * phi_norm2_)) * phix_ * phix_);
            case 2: // \Gamma_2: Hartley et al. 2020
                return (0.5 * I + ((phi_norm_ - sin_phi_norm_) / (phi_norm_ * phi_norm2_)) * phix_ + ((phi_norm2_ + 2.0 * cos_phi_norm_ - 2.0) / (2.0 * phi_norm2_ * phi_norm2_)) * phix_ * phix_);
            default: { // general cases
                Eigen::Matrix3d R = I + (sin_phi_norm_ / phi_norm_) * phix_ + ((1.0 - cos_phi_norm_) / phi_norm2_) * phix_ * phix_;
                Eigen::Matrix3d S = I;
                Eigen::Matrix3d phixk = I;
                long k_factorial = 1;
                for (int k = 1; k <= m; ++k) {
                    k_factorial = k_factorial * k;
                    phixk = (phixk * phix_).eval();
                    S = (S + (1.0 / static_cast<double>(k_factorial)) * phixk).eval();
                }
                // m is odd (&& m != 0) or even
                if (m % 2) {
                    // m is odd
                    return ((1.0 / static_cast<double>(k_factorial)) * I + (pow(-1.0, (m + 1.0) / 2.0) / pow(phi_norm_, m + 1.0)) * phix_ * (R - S));
                } else {
                    // m is even
                    return ((1.0 / static_cast<double>(k_factorial)) * I + (pow(-1.0,m / 2.0) / pow(phi_norm_, m)) * (R - S));
                }
            }
        }
    }

    // class SEKN
    SEKN::SEKN(unsigned int N, unsigned int K) : N_(N), K_(K) {
        // initialize \xi vector as zero vector
        xi_ = Eigen::VectorXd::Zero(N_ * (1 + K_));
        // initialize X matrix as identity matrix
        X_ = Eigen::MatrixXd::Identity(N_ + K_, N_ + K_);
    }

    SEKN::SEKN(unsigned int N) : N_(N) {
        //
    }

    SEKN::SEKN(const SEKN& sekn) {
        N_ = sekn.N_;
        K_ = sekn.K_;
        xi_ = sekn.xi_;
        X_ = sekn.X_;
    }

    // class SEK3
    SEK3::SEK3(unsigned int K) : SEKN(3, K) {
        // initialize gamma using phi, i.e. xi_.head(3)
        phi2Gamma();
    }

    SEK3::SEK3(Eigen::VectorXd xi) : SEKN(3) {
        xi_ = xi;
        K_ = (xi_.rows() - 3) / 3;
        // initialize gamma using phi, i.e. xi_.head(3)
        phi2Gamma();
        // initialize matrix Lie group X_ using initialized Lie algebra xi_
        expSEK3();
    }

    SEK3::SEK3(Eigen::MatrixXd X) : SEKN(3) {
        X_ = X;
        K_ = X_.cols() - 3;
        // initialize Lie algebra xi_ and gamma_ using initialized matrix Lie group X_
        // set gamma in the logSEK3()
        logSEK3();
    }

    SEK3::SEK3(const SEK3& sek3)  : SEKN(sek3) {
        gamma_ = sek3.gamma_;
    }

    Eigen::MatrixXd SEK3::getHat() const {
        // xi_hat is the hat martrix of SE_K(3)
        Eigen::MatrixXd xi_hat = Eigen::MatrixXd::Zero(3 + K_, 3 + K_);
        // rotation skew block
        xi_hat.block<3, 3>(0, 0) = gamma_.getPhix();
        // check K_ > 0, i.e. X_ is not SO(3)
        if (K_ > 0) {
            // vector hat block
            for (int i = 0; i < K_; ++i) {
                xi_hat.block<3, 1>(0, 3 + i) = xi_.segment<3>(3 + 3 * i);
            }
        }

        return xi_hat;
    }

    Eigen::MatrixXd SEK3::getAdjoint() const {
        // adjoint matrix of X_
        Eigen::MatrixXd Ad_X = Eigen::MatrixXd::Zero(3 + 3 * K_, 3 + 3 * K_);
        // rotation matrix
        Eigen::Matrix3d R = X_.block<3, 3>(0, 0);

        // assembling
        Ad_X.block<3, 3>(0, 0) = R;
        // check K_ > 0, i.e. X_ is not SO(3)
        if (K_ > 0) {
            for (int i = 0; i < K_; ++i) {
                Ad_X.block<3, 3>(3 + 3 * i, 3 + 3 * i) = R;
                Ad_X.block<3, 3>(3 + 3 * i, 0) = skew(X_.block<3, 1>(0, 3 + i)) * R;
            }
        }

        return Ad_X;
    }

    Eigen::Matrix3d SEK3::getLeftJacobianSO3() const {
        return gamma_.getGamma(1);
    }

    Eigen::Matrix3d SEK3::getRightJacobianSO3() const {
        return Gamma(-gamma_.getPhi()).getGamma(1);
    }

    SEK3 SEK3::getInverse() const {
        // inversion of X_
        Eigen::MatrixXd X_inv = Eigen::MatrixXd::Identity(X_.rows(), X_.cols());

        Eigen::Matrix3d RT = X_.block<3, 3>(0, 0).transpose();
        // rotation matrix block of X_inv
        X_inv.block<3, 3>(0, 0) = RT;
        // check K_ > 0, i.e. X_ is not SO(3)
        if (K_ > 0) {
            for (int i = 0; i < K_; ++i) {
                X_inv.block<3, 1>(0, 3 + i) = -RT * X_.block<3, 1>(0, 3 + i);
            }
        }

        return SEK3(X_inv);
    }

    void SEK3::setRotation(const Eigen::Matrix3d& R) {
        X_.block<3, 3>(0, 0) = R;
        // reset \xi and \Gamma
        // set gamma in the logSEK3()
        logSEK3();
    }

    void SEK3::setRotation(const SEK3& sek3) {
        this->setRotation(sek3.getRotationMatrix());
    }

    void SEK3::setRotationDelay(const Eigen::Matrix3d& R) {
        X_.block<3, 3>(0, 0) = R;
    }

    void SEK3::setRotationDelay(const SEK3& sek3) {
        this->setRotationDelay(sek3.getRotationMatrix());
    }

    void SEK3::setVector(unsigned int idx, const Eigen::Vector3d& vec) {
        // check index >= 3
        assert(idx >= 3);

        X_.block<3, 1>(0, idx) = vec;
        // reset a partial segment of \xi, \phi is not changed so \Gamma should not be changed
        Eigen::MatrixXd Jl_inv = this->getRightJacobianSO3().transpose();
        xi_.segment(3 * idx - 6, 3) = Jl_inv * X_.block<3, 1>(0, idx);
    }

    void SEK3::setVector(unsigned int idx, const SEK3& sek3) {
        this->setVector(idx, sek3.getVector(idx));
    }

    void SEK3::setVectorDelay(unsigned int idx, const Eigen::Vector3d& vec) {
        // check index >= 3
        assert(idx >= 3);

        X_.block<3, 1>(0, idx) = vec;
    }

    void SEK3::setVectorDelay(unsigned int idx, const SEK3& sek3) {
        this->setVectorDelay(idx, sek3.getVector(idx));
    }

    void SEK3::updateXiAndGamma() {
        // recall the private function logSEK3()
        this->logSEK3();
    }

    SEK3 SEK3::operator*(const SEK3& y) {
        // product on matrix Lie group
        Eigen::MatrixXd X_prod = this->X_ * y.X_;

        return SEK3(X_prod);
    }

    SEK3& SEK3::operator=(const SEK3& y) {
        K_ = y.K_;
        xi_ = y.xi_;
        X_ = y.X_;
        gamma_ = y.gamma_;
    }

    void SEK3::expSEK3() {
        X_ = Eigen::MatrixXd::Identity(3 + K_, 3 + K_);

        // rotation matrix block
        X_.block<3, 3>(0, 0) = gamma_.getGamma(0);
        // check K_ > 0, i.e. X_ is not SO(3)
        if (K_ > 0) {
            // left Jacobian J_l
            //Eigen::Matrix3d Jl = gamma_.getGamma(1);
            Eigen::Matrix3d Jl = this->getLeftJacobianSO3();
            // vector block
            for (int i = 0; i < K_; ++i) {
                X_.block<3, 1>(0, 3 + i) = Jl * xi_.segment<3>(3 + 3 * i);
            }
        }
    }

    void SEK3::logSEK3() {
        xi_ = Eigen::VectorXd::Zero(3 + 3 * K_);

        // angular vector segment
        Eigen::AngleAxisd angleaxis(X_.block<3, 3>(0, 0));
        xi_.head(3) = angleaxis.angle() * angleaxis.axis();

        // set \Gamma: set gamma using phi, i.e. xi_.head(3)
        phi2Gamma();

        // check K_ > 0, i.e. X_ is not SO(3)
        if (K_ > 0) {
            // inversion of left Jacobian = transposition of right Jacobian
            //Eigen::Matrix3d Jl_inv = gamma_.getGamma(1).inverse();
            //Eigen::Matrix3d Jl_inv = Gamma(-gamma_.getPhi()).getGamma(1).transpose();
            //Eigen::Matrix3d Jl_inv = RightJacobian_SO3(gamma_.getPhi()).transpose();
            Eigen::Matrix3d Jl_inv = this->getRightJacobianSO3().transpose();
            // vector segment
            for (int i = 0; i < K_; ++i) {
                xi_.segment<3>(3 + 3 * i) = Jl_inv * X_.block<3, 1>(0, 3 + i);
            }
        }
    }

    // class TN
    TN::TN(unsigned int N) : SEKN(N) {
        //
    }

    TN::TN(Eigen::VectorXd t) : SEKN(t.rows(), 1) {
        xi_.tail(N_) = t;
        X_.block(0, N_, N_, 1) = t;
    }

    TN::TN(Eigen::MatrixXd T) : SEKN(T.cols() - 1, 1) {
        xi_.tail(N_) = T.block(0, N_, N_, 1);
        X_.block(0, N_, N_, 1) = T.block(0, N_, N_, 1);
    }

    TN::TN(const TN& tn) : SEKN(tn) {
        //
    }

    Eigen::MatrixXd TN::getHat() const {
        Eigen::MatrixXd xi_hat = Eigen::MatrixXd::Zero(1 + N_, 1 + N_);
        xi_hat.block(0, N_, N_, 1) = xi_.tail(N_);

        return xi_hat;
    }

    TN TN::getInverse() const {
        Eigen::VectorXd t_inv = -xi_.tail(N_);

        return TN(t_inv);
    }

    void TN::setTVec(const Eigen::VectorXd& vec) {
        N_ = vec.rows();
        xi_ = Eigen::VectorXd::Zero(N_ + N_);
        xi_.tail(N_) = vec;
        X_ = Eigen::MatrixXd::Identity(1 + N_, 1 + N_);
        X_.block(0, N_, N_, 1) = vec;
    }

    void TN::setTMat(const Eigen::MatrixXd& mat) {
        N_ = mat.cols() - 1;
        X_ = mat;
        xi_ = Eigen::VectorXd::Zero(N_ + N_);
        xi_.tail(N_) = mat.block(0, N_, N_, 1);
    }

    TN TN::operator+(const TN& y) const {
        Eigen::VectorXd t_plus = xi_.tail(N_) + y.xi_.tail(N_);

        return TN(t_plus);
    }

    TN TN::operator-(const TN& y) const {
        Eigen::VectorXd t_minus = xi_.tail(N_) - y.xi_.tail(N_);

        return TN(t_minus);
    }

    TN TN::operator*(const double& k) const {
        Eigen::VectorXd t_mul = xi_.tail(N_) * k;

        return TN(t_mul);
    }

    TN& TN::operator=(const TN& y) {
        N_ = y.N_;
        xi_ = y.xi_;
        X_ = y.X_;
    }

}
