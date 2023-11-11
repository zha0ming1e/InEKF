#include "matrix_lie_group_template.h"

namespace inekf {

    long factorial(unsigned int n) {
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

        // m-th order Gamma series: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (\phi^{\wedge})^n
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
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
        for (int i = 0; i < K; ++i) {
            xi_hat.block<3, 1>(0, 3 + i) = xi.segment<3>(3 + 3 * i);
        }
//}

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
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
        for (int i = 0; i < K; ++i) {
            X.block<3, 1>(0, 3 + i) = Jl * xi.segment<3>(3 + 3 * i);
        }
//}

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
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
        for (int i = 0; i < K; ++i) {
            Ad_X.block<3, 3>(3 + 3 * i, 0) = skew(X.block<3, 1>(0, 3 + i)) * R;
            Ad_X.block<3, 3>(3 + 3 * i, 3 + 3 * i) = R;
        }
//}

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

    Eigen::Matrix3d Gamma::getGamma(unsigned int m) const {
        // check m >= 0
        assert(m >= 0);

        // m-th order Gamma series: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (\phi^{\wedge})^n
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

    Eigen::Matrix3d Gamma::getGammaNegativePhim1() const {
        return (Eigen::Matrix3d::Identity() + ((cos_phi_norm_ - 1.0) / phi_norm2_) * phix_ + ((phi_norm_ - sin_phi_norm_) / (phi_norm_ * phi_norm2_)) * phix_ * phix_);
    }

    // class SEKN
    template <unsigned int N, unsigned int K>
    SEKN<N, K>::SEKN() {
        // N dimensional space, N = 3 for SO(3) and SE_K(3)
        N_ = N;
        // K direct isometry, K = 0 for SO(3)
        K_ = K;
        // initialize \xi vector as zero vector
        xi_ = Eigen::VectorXd::Zero(N_ * (1 + K_));
        // initialize X matrix as identity matrix
        X_ = Eigen::MatrixXd::Identity(N_ + K_, N_ + K_);
    }

    template <unsigned int N, unsigned int K>
    SEKN<N, K>::SEKN(const SEKN<N, K>& sekn) {
        N_ = sekn.N_;
        K_ = sekn.K_;
        xi_ = sekn.xi_;
        X_ = sekn.X_;
    }

    // class SEK3
    template <unsigned int K>
    SEK3<K>::SEK3() : SEKN<3, K>() {
        // initialize gamma using phi, i.e. xi_.head(3)
        phi2Gamma();
    }

    template <unsigned int K>
    SEK3<K>::SEK3(Eigen::VectorXd xi) : SEKN<3, K>() {
        this->xi_ = xi;
        //this->K_ = (this->xi_.rows() - 3) / 3;  // unnecessary
        // initialize gamma using phi, i.e. xi_.head(3)
        phi2Gamma();
        // initialize matrix Lie group X_ using initialized Lie algebra xi_
        expSEK3();
    }

    template <unsigned int K>
    SEK3<K>::SEK3(Eigen::MatrixXd X) : SEKN<3, K>() {
        this->X_ = X;
        //this->K_ = this->X_.cols() - 3;  // unnecessary
        // initialize Lie algebra xi_ and gamma_ using initialized matrix Lie group X_
        // set gamma in the logSEK3()
        logSEK3();
    }

    template <unsigned int K>
    SEK3<K>::SEK3(const SEK3<K>& sek3) : SEKN<3, K>(sek3) {
        gamma_ = sek3.gamma_;
    }

    template <unsigned int K>
    Eigen::MatrixXd SEK3<K>::getHat() const {
        // xi_hat is the hat martrix of SE_K(3)
        Eigen::MatrixXd xi_hat = Eigen::MatrixXd::Zero(3 + this->K_, 3 + this->K_);
        // rotation skew block
        xi_hat.template block<3, 3>(0, 0) = gamma_.getPhix();
        // check K_ > 0, i.e. X_ is not SO(3)
        if (this->K_ > 0) {
            // vector hat block
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
            for (int i = 0; i < this->K_; ++i) {
                xi_hat.template block<3, 1>(0, 3 + i) = this->xi_.template segment<3>(3 + 3 * i);
            }
//}
        }

        return xi_hat;
    }

    template <unsigned int K>
    Eigen::MatrixXd SEK3<K>::getAdjoint() const {
        // adjoint matrix of X_
        Eigen::MatrixXd Ad_X = Eigen::MatrixXd::Zero(3 + 3 * this->K_, 3 + 3 * this->K_);
        // rotation matrix
        Eigen::Matrix3d R = this->X_.template block<3, 3>(0, 0);

        // assembling
        Ad_X.template block<3, 3>(0, 0) = R;
        // check K_ > 0, i.e. X_ is not SO(3)
        if (this->K_ > 0) {
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
            for (int i = 0; i < this->K_; ++i) {
                Ad_X.template block<3, 3>(3 + 3 * i, 3 + 3 * i) = R;
                Ad_X.template block<3, 3>(3 + 3 * i, 0) = skew(this->X_.template block<3, 1>(0, 3 + i)) * R;
            }
//}
        }

        return Ad_X;
    }

    template <unsigned int K>
    Eigen::Matrix3d SEK3<K>::getLeftJacobianSO3() const {
        return gamma_.getGamma(1);
    }

    template <unsigned int K>
    Eigen::Matrix3d SEK3<K>::getRightJacobianSO3() const {
        return gamma_.getGammaNegativePhim1();
    }

    template <unsigned int K>
    SEK3<K> SEK3<K>::getInverse() const {
        // inversion of X_
        Eigen::MatrixXd X_inv = Eigen::MatrixXd::Identity(this->X_.rows(), this->X_.cols());

        Eigen::Matrix3d RT = this->X_.template block<3, 3>(0, 0).transpose();
        // rotation matrix block of X_inv
        X_inv.template block<3, 3>(0, 0) = RT;
        // check K_ > 0, i.e. X_ is not SO(3)
        if (this->K_ > 0) {
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
            for (int i = 0; i < this->K_; ++i) {
                X_inv.template block<3, 1>(0, 3 + i) = -RT * this->X_.template block<3, 1>(0, 3 + i);
            }
//}
        }

        return SEK3(X_inv);
    }

    template <unsigned int K>
    void SEK3<K>::setRotation(const Eigen::Matrix3d& R) {
        this->X_.template block<3, 3>(0, 0) = R;
        // reset \xi and \Gamma
        // set gamma in the logSEK3()
        logSEK3();
    }

    template <unsigned int K>
    void SEK3<K>::setRotation(const SEK3<K>& sek3) {
        this->setRotation(sek3.getRotationMatrix());
    }

    template <unsigned int K>
    void SEK3<K>::setRotationDelay(const Eigen::Matrix3d& R) {
        this->X_.template block<3, 3>(0, 0) = R;
    }

    template <unsigned int K>
    void SEK3<K>::setRotationDelay(const SEK3<K>& sek3) {
        this->setRotationDelay(sek3.getRotationMatrix());
    }

    template <unsigned int K>
    void SEK3<K>::setVector(unsigned int idx, const Eigen::Vector3d& vec) {
        // check index >= 3
        assert(idx >= 3);

        this->X_.template block<3, 1>(0, idx) = vec;
        // reset a partial segment of \xi, \phi is not changed so \Gamma should not be changed
        Eigen::MatrixXd Jl_inv = this->getRightJacobianSO3().transpose();
        this->xi_.template segment<3>(3 * idx - 6) = Jl_inv * this->X_.template block<3, 1>(0, idx);
    }

    template <unsigned int K>
    void SEK3<K>::setVector(unsigned int idx, const SEK3<K>& sek3) {
        this->setVector(idx, sek3.getVector(idx));
    }

    template <unsigned int K>
    void SEK3<K>::setVectorDelay(unsigned int idx, const Eigen::Vector3d& vec) {
        // check index >= 3
        assert(idx >= 3);

        this->X_.template block<3, 1>(0, idx) = vec;
    }

    template <unsigned int K>
    void SEK3<K>::setVectorDelay(unsigned int idx, const SEK3<K>& sek3) {
        this->setVectorDelay(idx, sek3.getVector(idx));
    }

    template <unsigned int K>
    void SEK3<K>::updateXiAndGamma() {
        // recall the private function logSEK3()
        this->logSEK3();
    }

    template <unsigned int K>
    SEK3<K> SEK3<K>::operator*(const SEK3<K>& y) {
        // product on matrix Lie group
        Eigen::MatrixXd X_prod = this->X_ * y.X_;

        return SEK3<K>(X_prod);
    }

    template <unsigned int K>
    SEK3<K>& SEK3<K>::operator=(const SEK3<K>& y) {
        this->K_ = y.K_;
        this->xi_ = y.xi_;
        this->X_ = y.X_;
        gamma_ = y.gamma_;

        return *this;
    }

    template <unsigned int K>
    void SEK3<K>::expSEK3() {
        this->X_ = Eigen::MatrixXd::Identity(3 + this->K_, 3 + this->K_);

        // rotation matrix block
        this->X_.template block<3, 3>(0, 0) = gamma_.getGamma(0);
        // check K_ > 0, i.e. X_ is not SO(3)
        if (this->K_ > 0) {
            // left Jacobian J_l
            //Eigen::Matrix3d Jl = gamma_.getGamma(1);
            Eigen::Matrix3d Jl = this->getLeftJacobianSO3();
            // vector block
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
            for (int i = 0; i < this->K_; ++i) {
                this->X_.template block<3, 1>(0, 3 + i) = Jl * this->xi_.template segment<3>(3 + 3 * i);
            }
//}
        }
    }

    template <unsigned int K>
    void SEK3<K>::logSEK3() {
        this->xi_ = Eigen::VectorXd::Zero(3 + 3 * this->K_);

        // angular vector segment
        Eigen::AngleAxisd angleaxis(this->X_.template block<3, 3>(0, 0));
        this->xi_.template head<3>() = angleaxis.angle() * angleaxis.axis();

        // set \Gamma: set gamma using phi, i.e. xi_.head(3)
        phi2Gamma();

        // check K_ > 0, i.e. X_ is not SO(3)
        if (this->K_ > 0) {
            // inversion of left Jacobian = transposition of right Jacobian
            //Eigen::Matrix3d Jl_inv = gamma_.getGamma(1).inverse();
            //Eigen::Matrix3d Jl_inv = Gamma(-gamma_.getPhi()).getGamma(1).transpose();
            //Eigen::Matrix3d Jl_inv = RightJacobian_SO3(gamma_.getPhi()).transpose();
            Eigen::Matrix3d Jl_inv = this->getRightJacobianSO3().transpose();
            // vector segment
//{
//#pragma omp parallel for num_threads(8)
//#pragma omp parallel for
            for (int i = 0; i < this->K_; ++i) {
                this->xi_.template segment<3>(3 + 3 * i) = Jl_inv * this->X_.template block<3, 1>(0, 3 + i);
            }
//}
        }
    }

    // class TN
    template <unsigned int N>
    TN<N>::TN() : SEKN<N, 1>() {
        //
    }

    template <unsigned int N>
    TN<N>::TN(Eigen::VectorXd t) : SEKN<N, 1>() {
        this->xi_.template tail<this->N_>() = t;
        this->X_.template block<this->N_, 1>(0, this->N_) = t;
    }

    template <unsigned int N>
    TN<N>::TN(Eigen::MatrixXd T) : SEKN<N, 1>() {
        this->xi_.template tail<this->N_>() = T.template block<this->N_, 1>(0, this->N_);
        this->X_.template block<this->N_, 1>(0, this->N_) = T.template block<this->N_, 1>(0, this->N_);
    }

    template <unsigned int N>
    TN<N>::TN(const TN<N>& tn) : SEKN<N, 1>(tn) {
        //
    }

    template <unsigned int N>
    Eigen::MatrixXd TN<N>::getHat() const {
        Eigen::MatrixXd xi_hat = Eigen::MatrixXd::Zero(1 + this->N_, 1 + this->N_);
        xi_hat.template block<this->N_, 1>(0, this->N_) = this->xi_.template tail<this->N_>();

        return xi_hat;
    }

    template <unsigned int N>
    TN<N> TN<N>::getInverse() const {
        Eigen::VectorXd t_inv = -this->xi_.template tail<this->N_>();

        return TN(t_inv);
    }

    template <unsigned int N>
    void TN<N>::setTVec(const Eigen::VectorXd& vec) {
        // check dimensions
        assert(this->N_ == vec.rows());
        this->xi_ = Eigen::VectorXd::Zero(this->N_ + this->N_);
        this->xi_.template tail<this->N_>() = vec;
        this->X_ = Eigen::MatrixXd::Identity(1 + this->N_, 1 + this->N_);
        this->X_.template block<this->N_, 1>(0, this->N_) = vec;
    }

    template <unsigned int N>
    void TN<N>::setTMat(const Eigen::MatrixXd& mat) {
        // check dimensions
        assert(this->N_ == mat.cols() - 1);
        this->X_ = mat;
        this->xi_ = Eigen::VectorXd::Zero(this->N_ + this->N_);
        this->xi_.template tail<this->N_>() = mat.template block<this->N_, 1>(0, this->N_);
    }

    template <unsigned int N>
    TN<N> TN<N>::operator+(const TN<N>& y) const {
        Eigen::VectorXd t_plus = this->xi_.template tail<this->N_>() + y.xi_.template tail<this->N_>();

        return TN<N>(t_plus);
    }

    template <unsigned int N>
    TN<N> TN<N>::operator-(const TN<N>& y) const {
        Eigen::VectorXd t_minus = this->xi_.template tail<this->N_>() - y.xi_.template tail<this->N_>();

        return TN<N>(t_minus);
    }

    template <unsigned int N>
    TN<N> TN<N>::operator*(const double& k) const {
        Eigen::VectorXd t_mul = this->xi_.template tail<this->N_>() * k;

        return TN<N>(t_mul);
    }

    template <unsigned int N>
    TN<N>& TN<N>::operator=(const TN<N>& y) {
        // check dimensions
        assert(this->N_ == y.N_);
        this->xi_ = y.xi_;
        this->X_ = y.X_;

        return *this;
    }

}
