#ifndef INEKF_MATRIX_LIE_GROUP_H
#define INEKF_MATRIX_LIE_GROUP_H

#include "common.h"

namespace inekf {

    // factorial
    long factorial(unsigned int n);

    // 3D Lie algebra
    // skew-symmetric matrix of 3D Lie algebra
    Eigen::Matrix3d skew(const Eigen::Vector3d& phi);
    // Gamma series of 3D Lie algebra, Bloesch et al., 2012
    Eigen::Matrix3d Gamma_m(const Eigen::Vector3d& phi, int m);

    // SO(3)
    // hat map of SO(3)
    Eigen::MatrixXd Hat_SO3(const Eigen::Vector3d& phi);
    // exponential map of SO(3)
    Eigen::Matrix3d Exp_SO3(const Eigen::Vector3d& phi);
    // left Jacobian matrix of SO(3)
    Eigen::Matrix3d LeftJacobian_SO3(const Eigen::Vector3d& phi);
    // right Jacobian matrix of SO(3)
    Eigen::Matrix3d RightJacobian_SO3(const Eigen::Vector3d& phi);

    // SE_K(3)
    // hat map of SE_K(3)
    Eigen::MatrixXd Hat_SEK3(const Eigen::VectorXd& xi);
    // exponential map of SE_K(3)
    Eigen::MatrixXd Exp_SEK3(const Eigen::VectorXd& xi);
    // adjoint matrix of SE_K(3)
    Eigen::MatrixXd Adjoint_SEK3(const Eigen::MatrixXd& X);

    // class Gamma: \Gamma_m series
    // m-th order Gamma series: \Gamma_m = \sum_{n=0}^{\infty} \dfrac{1}{(n+m)!} (\phi^{\wedge})^n
    class Gamma {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        Gamma() = default;

        explicit Gamma(Eigen::Vector3d phi);

        Gamma(double theta1, double theta2, double theta3);

        explicit Gamma(double* theta);

        Gamma(const Gamma& gm);

        // destructor
        ~Gamma() = default;

        // functions
        Eigen::Matrix3d getGamma(unsigned int m) const;

        Eigen::Matrix3d getGammaNegativePhim1() const;

        Eigen::Vector3d getPhi() const { return phi_; }

        Eigen::Matrix3d getPhix() const { return phix_; }

        double getPhiNorm() const { return phi_norm_; }

        double getPhiNorm2() const { return phi_norm2_; }

        double getSinPhiNorm() const { return sin_phi_norm_; }

        double getCosPhiNorm() const { return cos_phi_norm_; }

    private:
        // angular vector \phi
        Eigen::Vector3d phi_;
        // skew of vector \phi
        Eigen::Matrix3d phix_;
        // norm of \phi
        double phi_norm_ = 0;
        // square norm of \phi
        double phi_norm2_ = 0;
        // sin of phi_norm
        double sin_phi_norm_ = 0;
        // cos of phi_norm
        double cos_phi_norm_ = 1;
    };

    // class SEKN: basic class of N dimensional K direct isometry matrix Lie group for SE_K(N)
    // including: SO(3) -> K = 0, N = 3
    // including: SE_K(3) -> K = K, N = 3
    // including: T(N) -> K = 1, N = N and top-left "rotation matrix" is N * N dimensional Identity matrix, i.e. I(N)
    // general representation: R_(N * N) is the general rotation matrix block
    // | R_(N * N) p1_(N * 1) ... pK_(N * 1) |
    // | 0_(1 * N)          1 ...          0 |
    // | 0_(1 * N)          0 ...          0 |
    // |         .          .   .          . |
    // |         .          .   .          . |
    // |         .          .   .          . |
    // | 0_(1 * N)          0 ...          1 |
    class SEKN {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        SEKN() = default;

        explicit SEKN(unsigned int N);

        SEKN(unsigned int N, unsigned int K);

        SEKN(const SEKN& sekn);

        // destructor
        ~SEKN() = default;

        // functions
        unsigned getN() const { return N_; }

        unsigned int getK() const { return K_; }

        Eigen::VectorXd getXi() const { return xi_; }

        Eigen::MatrixXd getX() const { return X_; }
        // get dimension of matrix Lie group X: X matrix columns
        unsigned int getXDim() const { return X_.cols(); }
        // get degree of freedom of matrix Lie group X: dimension of \xi or N * (1 + K)
        unsigned int getXDoF() const { return (N_ + N_ * K_); }
        // get rotation matrix block
        Eigen::Matrix3d getRotationMatrix() const { return X_.block(0, 0, N_, N_); }
        // get vector block from column index
        Eigen::Vector3d getVector(unsigned int idx) const { return X_.block(0, idx, N_, 1); }
        // get hat matrix of Lie algebra \xi
        virtual Eigen::MatrixXd getHat() const = 0;

    protected:
        // N dimensional space, N = 3 for SO(3) and SE_K(3)
        unsigned int N_ = 1;
        // K direct isometry, K = 0 for SO(3)
        unsigned int K_ = 0;
        // N * (1 + K) dimensional vector for Lie algebra
        Eigen::VectorXd xi_;
        // (N + K) * (N + K) dimensional matrix Lie group
        Eigen::MatrixXd X_;
    };

    // class SEK3: K direct isometry matrix Lie group SE_K(3), SO(3) is a special case for K = 0, i.e. SE_0(3)
    // general representation: R_(3 * 3) is the 3-dimensional rotation matrix block
    // | R_(3 * 3) p1_(3 * 1) ... pK_(3 * 1) |
    // | 0_(1 * 3)          1 ...          0 |
    // | 0_(1 * 3)          0 ...          0 |
    // |         .          .   .          . |
    // |         .          .   .          . |
    // |         .          .   .          . |
    // | 0_(1 * 3)          0 ...          1 |
    class SEK3 : public SEKN {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        // initialization sequence: K_ & xi_ -> gamma_ -> X_ or X_ -> K_ & -> gamma_ -> xi_
        explicit SEK3(unsigned int K);

        explicit SEK3(Eigen::VectorXd xi);

        explicit SEK3(Eigen::MatrixXd X);

        SEK3(const SEK3& sek3);

        // destructor
        ~SEK3() = default;

        // functions
        Gamma getGamma() const { return gamma_; }
        // get hat matrix of Lie algebra \xi
        Eigen::MatrixXd getHat() const override;
        // get adjoint matrix of matrix Lie group X
        Eigen::MatrixXd getAdjoint() const;
        // get left Jacobian of SO(3)
        Eigen::Matrix3d getLeftJacobianSO3() const;
        // get right Jacobian of SO(3)
        Eigen::Matrix3d getRightJacobianSO3() const;
        // get the inversion of SE_K(3)
        SEK3 getInverse() const;
        // set rotation matrix block from a rotation matrix
        void setRotation(const Eigen::Matrix3d& R);
        // set rotation matrix block from a SE_K(3) matrix Lie group
        void setRotation(const SEK3& sek3);
        // set rotation matrix block from a rotation matrix without updating \xi vector and \Gamma
        void setRotationDelay(const Eigen::Matrix3d& R);
        // set rotation matrix block from a SE_K(3) matrix Lie group without updating \xi vector and \Gamma
        void setRotationDelay(const SEK3& sek3);
        // set vector from a column index and a vector
        void setVector(unsigned int idx, const Eigen::Vector3d& vec);
        // set vector from a column index and a SE_K(3) matrix Lie group
        void setVector(unsigned int idx, const SEK3& sek3);
        // set vector from a column index and a vector without updating \xi vector and \Gamma
        void setVectorDelay(unsigned int idx, const Eigen::Vector3d& vec);
        // set vector from a column index and a SE_K(3) matrix Lie group without updating \xi vector and \Gamma
        void setVectorDelay(unsigned int idx, const SEK3& sek3);
        // update \xi vector and \Gamma using X's top-left 3 * 3 rotation matrix block -> \phi -> gamma -> \xi
        void updateXiAndGamma();
        // override operator *
        SEK3 operator*(const SEK3& y);
        // override operator =
        SEK3& operator=(const SEK3& y);

    private:
        // initialize class member \gamma using first 3 elements of Lie algebra \xi
        void phi2Gamma() { gamma_ = Gamma(xi_.head(3)); }
        // exponential map from Lie algebra vector \xi to matrix Lie group X
        void expSEK3();
        // logarithm map from matrix Lie group X to Lie algebra vector \xi
        void logSEK3();

    private:
        // fixed N = 3 for SE_K(3): N dimensional space
        const unsigned int N_ = 3;
        // \Gamma series
        Gamma gamma_;
    };

    // class TN: translation matrix Lie group T(N)
    // general representation:
    // | I_(N * N) t_(N * 1) |
    // | 0_(1 * N)         1 |
    class TN : public SEKN {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        explicit TN(unsigned int N);

        explicit TN(Eigen::VectorXd t);

        explicit TN(Eigen::MatrixXd T);

        TN(const TN& tn);

        // destructor
        ~TN() = default;

        // functions
        Eigen::VectorXd getTVec() const { return this->getXi().tail(N_); }

        Eigen::MatrixXd getTMat() const { return this->getX(); }
        // get hat matrix of Lie algebra \xi
        Eigen::MatrixXd getHat() const override;
        // get the inversion of T(N)
        TN getInverse() const;
        // set vector of T(N)
        void setTVec(const Eigen::VectorXd& vec);
        // set matrix of T(N)
        void setTMat(const Eigen::MatrixXd& mat);
        // override operator +
        TN operator+(const TN& y) const;
        // override operator -
        TN operator-(const TN& y) const;
        // override operator *
        TN operator*(const double& k) const;
        // override operator =
        TN& operator=(const TN& y);

    private:
        // K direct isometry, K = 1 for T(N)
        const unsigned int K_ = 1;
    };

}

#endif //INEKF_MATRIX_LIE_GROUP_H
