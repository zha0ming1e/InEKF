#ifndef INEKF_INEKF_H
#define INEKF_INEKF_H

#include "noise.h"
#include "matrix_lie_group.h"
#include "observation.h"
#include "state.h"

namespace inekf {

    // class Filter: basic class for filters
    class Filter {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        Filter() = default;

        explicit Filter(const Noise& noise);

        explicit Filter(const State& state);

        Filter(const State& state, const Noise& noise);

        Filter(const State& state, const Noise& noise, ErrorType error_type);

        Filter(const State& state, const Noise& noise, ErrorType error_type, IMURoleType imu_role_type);

        // destructor
        ~Filter() = default;

        // functions
        ErrorType getErrorType() const { return error_type_; }

        IMURoleType getIMURoleType() const { return imu_role_type_; }

        State getState() const { return state_; }

        Noise getNoise() const { return noise_; }

        void setErrorType(const ErrorType& error_type) { error_type_ = error_type; }

        void setIMURoleType(const IMURoleType& imu_role_type) { imu_role_type_ = imu_role_type; }

        void setState(const State& state) { state_ = state; }

        void setNoise(const Noise& noise) { noise_ = noise; }
        // reset and clear observation matrices: Z, H, N
        void resetObservationMatrices();

        // purly virtual functions
        // filter propagation and correction functions
        // propagate the state estimate and covariance in right invariant type using inertial measurements as input
        virtual void propagate(double dt, const Eigen::MatrixXd& input) = 0;
        // correct the state estimate using measurements and observation models with formulated Z, H, N matrices
        virtual void correct(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N) = 0;
        // discrete-time state transition matrix using input
        virtual Eigen::MatrixXd getStateTransitionMatrix(double dt, const Eigen::MatrixXd& input) = 0;
        // discrete-time noise covariance matrix using discrete-time state transition matrix
        virtual Eigen::MatrixXd getDiscreteNoiseMatrix(double dt, const Eigen::MatrixXd& Phi) = 0;

    protected:
        // error type: left invariant or right invariant
        ErrorType error_type_ = ErrorType::RightInvariant;
        // IMU role type: input or output
        IMURoleType imu_role_type_ = IMURoleType::INPUT;
        // state
        State state_;
        // noise
        Noise noise_;

    protected:
        // filter matrices
        // discrete-time state transition matrix: \Phi
        Eigen::MatrixXd Phi_;
        // discrete-time process noise covariance matrix: Q_d
        Eigen::MatrixXd Q_d_;
        // innovation matrix or vector: Z
        Eigen::MatrixXd Z_;
        // observation Jacobian matrix: H
        Eigen::MatrixXd H_;
        // observation noise covariance matrix: N
        Eigen::MatrixXd N_;
    };

    // class InEKF: invariant extended Kalman filter, i.e. Invariant-EKF
    class InEKF : public Filter {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        // constructors
        InEKF();

        explicit InEKF(const Noise& noise);

        explicit InEKF(const State& state);

        InEKF(const State& state, const Noise& noise);

        InEKF(const State& state, const Noise& noise, ErrorType error_type);

        InEKF(const State& state, const Noise& noise, ErrorType error_type, IMURoleType imu_role_type);

        // destructor
        ~InEKF() = default;

        // functions
        // get and set
        std::map<unsigned int, bool> getContacts() const { return contacts_; }

        std::map<unsigned int, unsigned int> getEstimatedContactPositions() const { return estimated_contact_positions_; }

        mapIntVector3d getPriorLandmarks() const { return prior_landmarks_; }

        std::map<unsigned int, unsigned int> getEstimatedLandmarks() const { return estimated_landmarks_; }

        Eigen::Vector3d getMagneticField() const { return magnetic_field_; }

        void setContacts(const std::vector<std::pair<unsigned int, bool>>& contacts);
        // set prior landmarks: static landmarks
        void setPriorLandmarks(const mapIntVector3d& prior_landmarks) { prior_landmarks_ = prior_landmarks; }

        void setMagneticField(const Eigen::Vector3d& true_magnetic_field) { magnetic_field_ = true_magnetic_field; }

        // utilities
        // reset and clear filter
        void clear();
        // remove a particular prior landmark using landmark index
        void removePriorLandmark(unsigned int landmark_id);
        // remove some prior landmarks using a vector of landmark indices
        void removePriorLandmarks(const std::vector<unsigned int>& landmark_ids);
        // remove a particular estimated landmark using landmark index
        void removeEstimatedLandmark(unsigned int landmark_id);
        // remove some estimated landmarks using a vector of landmark indices
        void removeEstimatedLandmarks(const std::vector<unsigned int>& landmark_ids);
        // keep some estimated landmarks using a vector of landmark indices
        void keepEstimatedLandmarks(const std::vector<unsigned int>& landmark_ids);

        // filter propagation and correction functions
        // propagate the state estimate and covariance in right invariant type using inertial measurements as input
        void propagate(double dt, const Eigen::MatrixXd& imu) override;
        // correct the state estimate using measurements and observation models with formulated Z, H, N matrices
        void correct(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N) override;
        // discrete-time state transition matrix using input angular velocity \omega and translational velocity a: \Phi
        Eigen::MatrixXd getStateTransitionMatrix(double dt, const Eigen::MatrixXd& imu) override;
        // discrete-time noise covariance matrix using discrete-time state transition matrix
        Eigen::MatrixXd getDiscreteNoiseMatrix(double dt, const Eigen::MatrixXd& Phi) override;

        // correct the state estimate and covariance in right invariant type using forward kinematics measurements
        void correctKinematics(const vectorKinematics& measured_kinematics);
        // correct the state estimate and covariance in right invariant type using relative landmark position measurements
        void correctLandmarks(const vectorLandmarks& measured_landmarks);
        // correct the state estimate and covariance in left invariant type using absolute z-position of contact point measurements
        void correctContactPosition(const unsigned int id, const Eigen::Vector3d& measured_contact_position, const Eigen::Matrix3d& covariance, const Eigen::Vector3d& indices);
        // correct the state estimate and covariance in right invariant type using magnetic field measurements
        //void correctMagnetometer(const Eigen::Vector3d& measured_magnetic_field, const Eigen::Matrix3d& covariance);
        // correct the state estimate and covariance in left invariant type using GPS absolute position measurements
        //void correctPosition(const Eigen::Vector3d& measured_position, const Eigen::Matrix3d& covariance, const Eigen::Vector3d& indices);

    public:
        // correct state using right invariant observation models
        //void correctRightInvariant(const Observation& obs);
        // correct state using right invariant observation models with formulated Z, H, N matrices
        void correctRightInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N);
        // correct state using left invariant observation models
        //void correctLeftInvariant(const Observation& obs);
        // correct state using left invariant observation models with formulated Z, H, N matrices
        void correctLeftInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N);

    private:
        // if estimate IMU biases
        bool estimate_bias_ = true;
        // gravity vector in the world frame (z-up)
        Eigen::Vector3d g_;
        // magnetic field vector
        Eigen::Vector3d magnetic_field_;
        // contact state: (index, if_contact)
        std::map<unsigned int, bool> contacts_;
        // estimated contact position: (index, column_index_in_state_matrix)
        std::map<unsigned int, unsigned int> estimated_contact_positions_;
        // prior landmark static vector: vector of (index, position)
        mapIntVector3d prior_landmarks_;
        // estimated landmark: (index, column_index_in_state_matrix)
        std::map<unsigned int, unsigned int> estimated_landmarks_;
    };

}

#endif //INEKF_INEKF_H
