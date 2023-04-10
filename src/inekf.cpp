#include "inekf.h"

namespace inekf {

    // class Filter
    Filter::Filter(const Noise& noise) : noise_(noise) {
        //
    }

    Filter::Filter(const State& state) : state_(state) {
        //
    }

    Filter::Filter(const State& state, const Noise& noise) : state_(state), noise_(noise) {
        //
    }

    Filter::Filter(const State& state, const Noise& noise, ErrorType error_type) : state_(state), noise_(noise), error_type_(error_type) {
        //
    }

    Filter::Filter(const State& state, const Noise& noise, ErrorType error_type, IMURoleType imu_role_type)
    : state_(state), noise_(noise), error_type_(error_type), imu_role_type_(imu_role_type) {
        //
    }

    void Filter::resetObservationMatrices() {
        // reset and clear observation matrices: Z, H, N
        Z_ = Eigen::MatrixXd();
        H_ = Eigen::MatrixXd();
        N_ = Eigen::MatrixXd();
    }

    // class InEKF
    InEKF::InEKF() : g_(Eigen::Vector3d(0, 0, -9.81)), magnetic_field_(Eigen::Vector3d::Zero()) {
        //
    }

    InEKF::InEKF(const Noise& noise) : Filter(noise),
    g_(Eigen::Vector3d(0, 0, -9.81)),
    magnetic_field_(Eigen::Vector3d(std::cos(1.2049), 0, std::sin(1.2049))) {
        //
    }

    InEKF::InEKF(const State& state) : Filter(state),
    g_(Eigen::Vector3d(0, 0, -9.81)),
    magnetic_field_(Eigen::Vector3d(std::cos(1.2049), 0, std::sin(1.2049))) {
        //
    }

    InEKF::InEKF(const State& state, const Noise& noise) : Filter(state, noise),
    g_(Eigen::Vector3d(0, 0, -9.81)),
    magnetic_field_(Eigen::Vector3d(std::cos(1.2049), 0, std::sin(1.2049))) {
        //
    }

    InEKF::InEKF(const State& state, const Noise& noise, ErrorType error_type) : Filter(state, noise, error_type),
    g_(Eigen::Vector3d(0, 0, -9.81)),
    magnetic_field_(Eigen::Vector3d(std::cos(1.2049), 0, std::sin(1.2049))) {
        //
    }

    InEKF::InEKF(const State& state, const Noise& noise, ErrorType error_type, IMURoleType imu_role_type)
    : Filter(state, noise, error_type, imu_role_type),
    g_(Eigen::Vector3d(0, 0, -9.81)),
    magnetic_field_(Eigen::Vector3d(std::cos(1.2049), 0, std::sin(1.2049))) {
        //
    }

    void InEKF::setContacts(const std::vector<std::pair<unsigned int, bool>>& contacts) {
        // insert new measured contact states
        for (auto& it : contacts) {
            auto return_label = contacts_.insert(it);
            // if the contact is already in the contacts_, re-assign it with the new measured state
            if (!return_label.second)
                return_label.first->second = it.second;
        }
    }

    void InEKF::clear() {
        state_ = State();
        noise_ = Noise();

        Phi_ = Eigen::MatrixXd();
        Q_d_ = Eigen::MatrixXd();

        Z_ = Eigen::MatrixXd();
        H_ = Eigen::MatrixXd();
        N_ = Eigen::MatrixXd();

        contacts_.clear();
        estimated_contact_positions_.clear();
        prior_landmarks_.clear();
        estimated_landmarks_.clear();
    }

    void InEKF::removePriorLandmark(const unsigned int landmark_id) {
        // search the landmark id in the prior landmarks
        auto it = prior_landmarks_.find(landmark_id);
        if (it != prior_landmarks_.end()) {
            // the remove from the prior landmarks
            prior_landmarks_.erase(it->first);
        }
    }

    void InEKF::removePriorLandmarks(const std::vector<unsigned int>& landmark_ids) {
        for (auto& id : landmark_ids) {
            this->removePriorLandmark(id);
        }
    }

    void InEKF::removeEstimatedLandmark(const unsigned int landmark_id) {
        // search the landmark id in the estimated landmarks
        auto it = estimated_landmarks_.find(landmark_id);
        if (it != estimated_landmarks_.end()) {
            // get current state X and covariance P
            Eigen::MatrixXd X_rem = state_.getXMatrix();
            Eigen::MatrixXd P_rem = state_.getP();
            // remove corresponding row and column from X
            removeMatrixRowAndColumn(it->second, X_rem, 1);
            // remove 3 rows and columns from P
            unsigned int startIdx = 3 + 3 * (it->second - 3);
            removeMatrixRowAndColumn(startIdx, P_rem, 3);
            //removeMatrixRowAndColumn(startIdx, P_rem);
            //removeMatrixRowAndColumn(startIdx, P_rem);
            //removeMatrixRowAndColumn(startIdx, P_rem);
            // update all indices for estimated_landmarks_ and estimated_contact_positions_
            for (auto& it2 : estimated_landmarks_) {
                if (it2.second > it->second)
                    it2.second -= 1;
            }
            for (auto& it2 : estimated_contact_positions_) {
                if (it2.second > it->second)
                    it2.second -= 1;
            }
            // remove from estimated landmark positions
            estimated_landmarks_.erase(it->first);
            // update state X and covariance P
            state_.setX(X_rem);
            state_.setP(P_rem);
        }
    }

    void InEKF::removeEstimatedLandmarks(const std::vector<unsigned int>& landmark_ids) {
        for (auto& it : landmark_ids) {
            this->removeEstimatedLandmark(it);
        }
    }

    void InEKF::keepEstimatedLandmarks(const std::vector<unsigned int>& landmark_ids) {
        std::vector<unsigned int> ids_to_erase;
        for (auto& it : estimated_landmarks_) {
            auto it_found = std::find(landmark_ids.begin(), landmark_ids.end(), it.first);
            if (it_found == landmark_ids.end()) {
                // get current state X and covariance P
                Eigen::MatrixXd X_rem = state_.getXMatrix();
                Eigen::MatrixXd P_rem = state_.getP();
                // remove corresponding row and column from X
                removeMatrixRowAndColumn(it.second, X_rem, 1);
                // remove 3 rows and columns from P
                unsigned int startIdx = 3 + 3 * (it.second - 3);
                removeMatrixRowAndColumn(startIdx, P_rem, 3);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                // update all indices for estimated_landmarks_ and estimated_contact_positions_
                for (auto& it2 : estimated_landmarks_) {
                    if (it2.second > it.second)
                        it2.second -= 1;
                }
                for (auto& it2 : estimated_contact_positions_) {
                    if (it2.second > it.second)
                        it2.second -= 1;
                }
                // add to ids_to_erase
                ids_to_erase.push_back(it.first);
                // update state X and covariance P
                state_.setX(X_rem);
                state_.setP(P_rem);
            }
        }
        // remove from estimated_landmarks_
        for (auto& idit : ids_to_erase) {
            estimated_landmarks_.erase(idit);
        }
    }

    Eigen::MatrixXd InEKF::getStateTransitionMatrix(const double dt, const Eigen::MatrixXd& imu) {
        // IMU measurements: angular velocity (w) and specific force (a) in the body frame
        Eigen::Vector3d w = imu.block<3, 1>(0, 0);
        Eigen::Vector3d a = imu.block<3, 1>(3, 0);

        Eigen::Vector3d phi = w * dt;
        // \Gamma with angular increment \phi
        Gamma gamma = Gamma(phi);
        Eigen::Matrix3d G0 = gamma.getGamma(0);
        Eigen::Matrix3d G0T = G0.transpose();
        Eigen::Matrix3d G1 = gamma.getGamma(1);
        Eigen::Matrix3d G1T = G1.transpose();
        Eigen::Matrix3d G2 = gamma.getGamma(2);
        Eigen::Matrix3d G2T = G2.transpose();
        Eigen::Matrix3d G3T = Gamma(-phi).getGamma(3);

        unsigned int dimX = state_.getXDim();
        unsigned int dimTheta = state_.getThetaDim();
        unsigned int dimP = state_.getPDim();
        Eigen::MatrixXd Phi = Eigen::MatrixXd::Identity(dimP, dimP);

        // compute the complicated bias terms (derived for the left invariant case)
        Eigen::Matrix3d ax = skew(a);
        Eigen::Matrix3d wx = skew(w);
        Eigen::Matrix3d wx2 = wx * wx;

        double dt2 = dt * dt;
        double dt3 = dt2 * dt;

        double theta = w.norm();
        double theta2 = theta * theta;
        double theta3 = theta2 * theta;
        double theta4 = theta3 * theta;
        double theta5 = theta4 * theta;
        double theta6 = theta5 * theta;
        double theta7 = theta6 * theta;

        double thetadt = theta * dt;
        double thetadt2 = thetadt * thetadt;
        double thetadt3 = thetadt2 * thetadt;
        double sinthetadt = sin(thetadt);
        double costhetadt = cos(thetadt);
        double sin2thetadt = sin(2 * thetadt);
        double cos2thetadt = cos(2 * thetadt);
        double thetadtcosthetadt = thetadt * costhetadt;
        double thetadtsinthetadt = thetadt * sinthetadt;

        Eigen::Matrix3d Phi25L = G0T * (ax * G2T * dt2
                + ((sinthetadt - thetadtcosthetadt) / theta3) * (wx * ax)
                - ((cos2thetadt - 4 * costhetadt + 3) / (4 * theta4)) * (wx * ax * wx)
                + ((4 * sinthetadt + sin2thetadt - 4 * thetadtcosthetadt - 2 * thetadt) / (4 * theta5)) * (wx * ax * wx2)
                + ((thetadt2 - 2 * thetadtsinthetadt - 2 * costhetadt + 2) / (2 * theta4)) * (wx2 * ax)
                - ((6 * thetadt - 8 * sinthetadt + sin2thetadt) / (4 * theta5)) * (wx2 * ax * wx)
                + ((2 * thetadt2 - 4 * thetadtsinthetadt - cos2thetadt + 1) / (4 * theta6)) * (wx2 * ax * wx2));

        Eigen::Matrix3d Phi35L = G0T * (ax * G3T * dt3
                - ((thetadtsinthetadt + 2 * costhetadt - 2) / theta4) * (wx * ax)
                - ((6 * thetadt - 8 * sinthetadt + sin2thetadt) / (8 * theta5)) * (wx * ax * wx)
                - ((2 * thetadt2 + 8 * thetadtsinthetadt + 16 * costhetadt + cos2thetadt - 17) / (8 * theta6)) * (wx * ax * wx2)
                + ((thetadt3 + 6 * thetadt - 12 * sinthetadt + 6 * thetadtcosthetadt) / (6 * theta5)) * (wx2 * ax)
                - ((6 * thetadt2 + 16 * costhetadt - cos2thetadt - 15) / (8 * theta6)) * (wx2 * ax * wx)
                + ((4 * thetadt3 + 6 * thetadt - 24 * sinthetadt - 3 * sin2thetadt + 24 * thetadtcosthetadt) / (24 * theta7)) * (wx2 * ax * wx2));
        // numerical tolerance for approximation
        const double tol =  1.0e-6;
        if (theta < tol) {
            Phi25L = (1 / 2) * ax * dt2;
            Phi35L = (1 / 6) * ax * dt3;
        }

        // fill the analytical discrete state transition matrices
        if ((state_.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::LeftInvariant) ||
            (state_.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::RightInvariant)) {
            // compute left-invariant state transition matrix \Phi
            // \Phi_{11}
            Phi.block<3, 3>(0, 0) = G0T;
            // \Phi_{21}
            Phi.block<3, 3>(3, 0) = -G0T * skew(G1 * a) * dt;
            // \Phi_{31}
            Phi.block<3, 3>(6, 0) = -G0T * skew(G2 * a) * dt2;
            // \Phi_{22}
            Phi.block<3, 3>(3, 3) = G0T;
            // \Phi_{32}
            Phi.block<3, 3>(6, 3) = G0T * dt;
            // \Phi_{33}
            Phi.block<3, 3>(6, 6) = G0T;
            for (int i = 5; i < dimX; ++i) {
                // \Phi_{(3 + i)(3 + i)}
                Phi.block<3, 3>((i - 2) * 3, (i - 2) * 3) = G0T;
            }
            // \Phi_{15}
            Phi.block<3, 3>(0, dimP - dimTheta) = -G1T * dt;
            // \Phi_{25}
            Phi.block<3, 3>(3, dimP - dimTheta) = Phi25L;
            // \Phi_{35}
            Phi.block<3, 3>(6, dimP - dimTheta) = Phi35L;
            // \Phi_{26}
            Phi.block<3, 3>(3, dimP - dimTheta + 3) = -G1T * dt;
            // \Phi_{36}
            Phi.block<3, 3>(6, dimP - dimTheta + 3) = -G0T * G2 * dt2;
        } else {
            // compute right-invariant state transition matrix
            Eigen::Matrix3d gx = skew(g_);
            Eigen::Matrix3d R = state_.getRotation();
            Eigen::Vector3d v = state_.getVelocity();
            Eigen::Vector3d p = state_.getPosition();
            Eigen::Matrix3d RG0 = R * G0;
            Eigen::Matrix3d RG1dt = R * G1 * dt;
            Eigen::Matrix3d RG2dt2 = R * G2 * dt2;
            // \Phi_{21}
            Phi.block<3, 3>(3, 0) = gx * dt;
            // \Phi_{31}
            Phi.block<3, 3>(6, 0) = 0.5 * gx * dt2;
            // \Phi_{32}
            Phi.block<3, 3>(6, 3) = Eigen::Matrix3d::Identity() * dt;
            // \Phi_{15}
            Phi.block<3, 3>(0, dimP - dimTheta) = -RG1dt;
            // \Phi_{25}
            Phi.block<3, 3>(3, dimP - dimTheta) = -skew(v + RG1dt * a + g_ * dt) * RG1dt + RG0 * Phi25L;
            // \Phi_{35}
            Phi.block<3, 3>(6, dimP - dimTheta) = -skew(p + v * dt + RG2dt2 * a + 0.5 * g_ * dt2) * RG1dt + RG0 * Phi35L;
            for (int i = 5; i < dimX; ++i) {
                // \Phi_{(3 + i)5}
                Phi.block<3, 3>((i - 2) * 3, dimP - dimTheta) = -skew(state_.getVector(i)) * RG1dt;
            }
            // \Phi_{26}
            Phi.block<3, 3>(3, dimP - dimTheta + 3) = -RG1dt;
            // \Phi_{36}
            Phi.block<3, 3>(6, dimP - dimTheta + 3) = -RG2dt2;
        }
        return Phi;
    }

    Eigen::MatrixXd InEKF::getDiscreteNoiseMatrix(const double dt, const Eigen::MatrixXd& Phi) {
        unsigned int dimTheta = state_.getThetaDim();
        unsigned int dimP = state_.getPDim();
        Eigen::MatrixXd G = Eigen::MatrixXd::Identity(dimP, dimP);

        // compute G using adjoint of X if needed, otherwise identity
        if ((state_.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::RightInvariant) ||
            (state_.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::LeftInvariant)) {
            G.block(0, 0, dimP - dimTheta, dimP - dimTheta) = state_.getWorldX().getAdjoint();
        }

        // continuous-time noise covariance: landmark noise terms will remain zero
        Eigen::MatrixXd Qc = Eigen::MatrixXd::Zero(dimP, dimP);
        Qc.block<3, 3>(0, 0) = noise_.getGyroscopeCov();
        Qc.block<3, 3>(3, 3) = noise_.getAccelerometerCov();
        // contact noise terms
        for (auto& it : estimated_contact_positions_) {
            Qc.block<3, 3>(3 + 3 * (it.second - 3), 3 + 3 * (it.second - 3)) = noise_.getContactCov();
        }
        // IMU bias noise terms
        Qc.block<3, 3>(dimP - dimTheta, dimP - dimTheta) = noise_.getGyroscopeBiasCov();
        Qc.block<3, 3>(dimP - dimTheta + 3, dimP - dimTheta + 3) = noise_.getAccelerometerBiasCov();

        // discrete noise covariance
        Eigen::MatrixXd PhiG = Phi * G;
        Eigen::MatrixXd Qd = PhiG * Qc * PhiG.transpose() * dt;

        return Qd;
    }

    void InEKF::propagate(double dt, const Eigen::MatrixXd &imu) {
        // un-biased IMU measurements: angular velocity (w) and specific force (a) in the body frame
        Eigen::MatrixXd imu_unbiased = imu.block<6, 1>(0, 0) - state_.getTheta();
        Eigen::Vector3d w = imu_unbiased.block<3, 1>(0, 0);
        Eigen::Vector3d a = imu_unbiased.block<3, 1>(3, 0);

        // get current dimensions in the filter
        auto dimX = state_.getXDim();
        auto dimP = state_.getPDim();
        auto dimTheta = state_.getThetaDim();

        // propagate: propagate error covariance using state transition matrix \Phi and discrete-time process noise Q_d
        Phi_ = this->getStateTransitionMatrix(dt, imu_unbiased);
        Q_d_ = this->getDiscreteNoiseMatrix(dt, Phi_);
        Eigen::MatrixXd P_pred = Phi_ * state_.getP() * Phi_.transpose() + Q_d_;

        // remove IMU bias related corrections if the system do not estimate IMU biases
        // set co-relative terms = 0, set auto-relative terms = 1
        if (!estimate_bias_) {
            P_pred.block(0, dimP - dimTheta, dimP - dimTheta, dimTheta) = Eigen::MatrixXd::Zero(dimP - dimTheta, dimTheta);
            P_pred.block(dimP - dimTheta, 0, dimTheta, dimP - dimTheta) = Eigen::MatrixXd::Zero(dimTheta, dimP - dimTheta);
            P_pred.block(dimP - dimTheta, dimP - dimTheta, dimTheta, dimTheta) = Eigen::MatrixXd::Zero(dimTheta, dimTheta);
        }

        // propagate nominal state estimates, i.e. estimate means
        // current nominal state estimates
        auto R = state_.getRotation();
        auto v = state_.getVelocity();
        auto p = state_.getPosition();
        // angular increment
        Eigen::Vector3d phi = w * dt;
        // \Gamma series
        Gamma gamma = Gamma(phi);
        Eigen::Matrix3d G0 = gamma.getGamma(0);
        Eigen::Matrix3d G1 = gamma.getGamma(1);
        Eigen::Matrix3d G2 = gamma.getGamma(2);
        // propagate nominal states
        auto X_pred = state_.getX();
        if (state_.getStateType() == StateType::WorldCentric) {
            // propagate world-centric state estimates
            // nominal rotation state
            //X_pred.setRotation(R * G0);
            // nominal velocity state
            //X_pred.setVector(3, v + (R * G1 * a + g_) * dt);
            // nominal position state
            //X_pred.setVector(4, p + v * dt + (R * G2 * a + 0.5 * g_) * dt * dt);

            // or: set rotation and vectors with delay and update Lie algebra vector \xi manually
            // avoiding multi-times updating computation
            X_pred.setRotationDelay(R * G0);
            X_pred.setVectorDelay(3, v + (R * G1 * a + g_) * dt);
            X_pred.setVectorDelay(4, p + v * dt + (R * G2 * a + 0.5 * g_) * dt * dt);
            X_pred.updateXiAndGamma();
        } else {
            // propagate robot body-centric state estimates
            Eigen::Matrix3d G0T = G0.transpose();
            // nominal rotation state
            //X_pred.setRotation(G0T * R);
            // nominal velocity state
            //X_pred.setVector(3, G0T * (v - (G1 * a + R * g_) * dt));
            // nominal position state
            //X_pred.setVector(4, G0T * (p + v * dt - (G2 * a + 0.5 * R * g_) * dt * dt));
            //for (int i = 5; i < dimX; ++i) {
            //    X_pred.setVector(i, G0T * state_.getX().getVector(i));
            //}

            // or: set rotation and vectors with delay and update Lie algebra vector \xi manually
            // avoiding multi-times updating computation
            X_pred.setRotationDelay(G0T * R);
            X_pred.setVectorDelay(3, G0T * (v - (G1 * a + R * g_) * dt));
            X_pred.setVectorDelay(4, G0T * (p + v * dt - (G2 * a + 0.5 * R * g_) * dt * dt));
            for (int i = 5; i < dimX; ++i) {
                X_pred.setVectorDelay(i, G0T * state_.getX().getVector(i));
            }
            X_pred.updateXiAndGamma();
        }

        // update nominal state estimates
        state_.setX(X_pred);
        state_.setP(P_pred);
    }

    void InEKF::correct(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N) {
        // depending on the error type: left invariant or right invariant
        if (error_type_ == ErrorType::RightInvariant) {
            this->correctRightInvariant(Z, H, N);
        } else {
            this->correctLeftInvariant(Z, H, N);
        }
    }

    void InEKF::correctRightInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N) {
        // get current dimensions in the filter
        auto dimTheta = state_.getThetaDim();
        auto dimP = state_.getPDim();
        // get current state covariance
        auto P = state_.getP();

        // remove IMU biases \theta from the matrix
        Eigen::MatrixXd Theta = Eigen::MatrixXd::Zero(6, 1);
        P.block<6, 6>(dimP - dimTheta, dimP - dimTheta) = 1.0e-4 * Eigen::MatrixXd::Identity(6, 6);
        P.block(0, dimP - dimTheta, dimP - dimTheta, dimTheta) = Eigen::MatrixXd::Zero(dimP - dimTheta, dimTheta);
        P.block(dimP - dimTheta, 0, dimTheta, dimP - dimTheta) = Eigen::MatrixXd::Zero(dimTheta, dimP - dimTheta);

        // transform left invariant error to right invariant error using adjoint map
        if (error_type_ == ErrorType::LeftInvariant) {
            // augment adjoint matrix including \theta rows and columns
            // | Ad_X 0 |
            // |    0 I |
            Eigen::MatrixXd Ad_aug = Eigen::MatrixXd::Identity(dimP, dimP);
            Ad_aug.block(0, 0, dimP - dimTheta, dimP - dimTheta) = state_.getX().getAdjoint();
            P = (Ad_aug * P * Ad_aug.transpose()).eval();
        }

        // Kalman gain
        Eigen::MatrixXd PHT = P * H.transpose();
        Eigen::MatrixXd S = H * PHT + N;
        Eigen::MatrixXd K = PHT * S.inverse();

        // state correction vector using Kalman gain and innovation vector
        Eigen::VectorXd delta = K * Z;
        SEK3 dX = SEK3(Eigen::VectorXd(delta.segment(0, delta.rows() - dimTheta)));
        Eigen::VectorXd dTheta = delta.segment(delta.rows() - dimTheta, dimTheta);

        // correct state estimates using right invariant error type
        SEK3 X_new = dX * state_.getX();
        Eigen::VectorXd Theta_new = Theta + dTheta;

        // update state estimates and error covariance
        // state estimates
        state_.setX(X_new);
        state_.setTheta(Theta_new);
        // covariance: Joseph update form
        Eigen::MatrixXd IKH = Eigen::MatrixXd::Identity(dimP, dimP) - K * H;
        Eigen::MatrixXd P_new = IKH * P * IKH.transpose() + K * N * K.transpose();
        // transform right invariant error to left invariant error using adjoint map
        if (error_type_ == ErrorType::LeftInvariant) {
            Eigen::MatrixXd Ad_aug_Inv = Eigen::MatrixXd::Identity(dimP, dimP);
            Ad_aug_Inv.block(0, 0, dimP - dimTheta, dimP - dimTheta) = state_.getXInverse().getAdjoint();
            P_new = (Ad_aug_Inv * P_new * Ad_aug_Inv.transpose()).eval();
        }
        state_.setP(P_new);
    }

    void InEKF::correctLeftInvariant(const Eigen::MatrixXd& Z, const Eigen::MatrixXd& H, const Eigen::MatrixXd& N) {
        // get current dimensions in the filter
        auto dimTheta = state_.getThetaDim();
        auto dimP = state_.getPDim();
        // get current state covariance
        auto P = state_.getP();

        // transform right invariant error to left invariant error using adjoint map
        if (error_type_ == ErrorType::RightInvariant) {
            // augment adjoint matrix including \theta rows and columns
            // | Ad_X  0 |
            // |    0  I |
            Eigen::MatrixXd Ad_aug_Inv = Eigen::MatrixXd::Identity(dimP, dimP);
            Ad_aug_Inv.block(0, 0, dimP - dimTheta, dimP - dimTheta) = state_.getXInverse().getAdjoint();
            P = (Ad_aug_Inv * P * Ad_aug_Inv.transpose()).eval();
        }

        // Kalman gain
        Eigen::MatrixXd PHT = P * H.transpose();
        Eigen::MatrixXd S = H * PHT + N;
        Eigen::MatrixXd K = PHT * S.inverse();

        // state correction vector using Kalman gain and innovation vector
        Eigen::VectorXd delta = K * Z;
        SEK3 dX = SEK3(Eigen::VectorXd(delta.segment(0, delta.rows() - dimTheta)));
        Eigen::VectorXd dTheta = delta.segment(delta.rows() - dimTheta, dimTheta);

        // correct state estimates using left invariant error type
        SEK3 X_new = state_.getX() * dX;
        Eigen::VectorXd Theta_new = state_.getTheta() + dTheta;

        // update state estimates and error covariance
        // state estimates
        state_.setX(X_new);
        state_.setTheta(Theta_new);
        // covariance: Joseph update form
        Eigen::MatrixXd IKH = Eigen::MatrixXd::Identity(dimP, dimP) - K * H;
        Eigen::MatrixXd P_new = IKH * P * IKH.transpose() + K * N * K.transpose();
        // transform left invariant error to right invariant error using adjoint map
        if (error_type_ == ErrorType::RightInvariant) {
            Eigen::MatrixXd Ad_aug = Eigen::MatrixXd::Identity(dimP, dimP);
            // adjoint matrix at the new state X_new
            Ad_aug.block(0, 0, dimP - dimTheta, dimP - dimTheta) = X_new.getAdjoint();
            P_new = (Ad_aug * P_new * Ad_aug.transpose()).eval();
        }
        state_.setP(P_new);
    }

    void InEKF::correctKinematics(const vectorKinematics& measured_kinematics) {
        // reset and clear observation matrices: Z, H, N
        this->resetObservationMatrices();

        // contacts need to be removed
        std::vector<std::pair<unsigned int, unsigned int>> remove_contacts;
        // new contacts
        vectorKinematics new_contacts;
        // used contact indices
        std::vector<unsigned int> used_contact_ids;

        for (auto& it : measured_kinematics) {
            // detect and skip if an ID is not unique
            if (std::find(used_contact_ids.begin(), used_contact_ids.end(), it.id) != used_contact_ids.end()) {
                std::cout << "Duplicate contact ID detected! Skipping measurement... \n";
                continue;
            } else {
                used_contact_ids.push_back(it.id);
            }

            // find contact indicator for kinematics measurements
            auto it_contact = contacts_.find(it.id);
            // skip if contact state is unknown
            if (it_contact == contacts_.end()) {
                continue;
            }
            bool contact_indicated = it_contact->second;

            // if the id in the estimated_contact_positions_
            auto it_estimated = estimated_contact_positions_.find(it.id);
            bool found = it_estimated != estimated_contact_positions_.end();

            // if contact is not indicated and id is found in estimated_contacts_, then remove state
            if (!contact_indicated && found) {
                remove_contacts.emplace_back(*it_estimated);
            } else if (contact_indicated && !found) {
                // if contact is indicated and id is not found in the estimated_contacts_, then augment state
                new_contacts.push_back(it);
            } else if (contact_indicated) {
                // if contact is indicated and id is found in the estimated_contacts_, then correct using kinematics
                auto dimTheta = state_.getThetaDim();
                auto dimP = state_.getPDim();
                unsigned int startIdx;

                // fill H
                startIdx = H_.rows();
                H_.conservativeResize(startIdx + 3, dimP);
                H_.block(startIdx, 0, 3, dimP) = Eigen::MatrixXd::Zero(3, dimP);
                if (state_.getStateType() == StateType::WorldCentric) {
                    // -I
                    H_.block(startIdx, 6, 3, 3) = -Eigen::Matrix3d::Identity();
                    // I
                    H_.block(startIdx, 3 * it_estimated->second - dimTheta, 3, 3) = Eigen::Matrix3d::Identity();
                } else {
                    // I
                    H_.block(startIdx, 6, 3, 3) = Eigen::Matrix3d::Identity();
                    // -I
                    H_.block(startIdx, 3 * it_estimated->second - dimTheta, 3, 3) = -Eigen::Matrix3d::Identity();
                }

                // fill N
                startIdx = N_.rows();
                N_.conservativeResize(startIdx + 3, startIdx + 3);
                N_.block(startIdx, 0, 3, startIdx) = Eigen::MatrixXd::Zero(3, startIdx);
                N_.block(0, startIdx, startIdx, 3) = Eigen::MatrixXd::Zero(startIdx, 3);
                N_.block(startIdx, startIdx, 3, 3) = state_.getWorldRotation() * it.covariance.block<3, 3>(3, 3) * state_.getWorldRotation().transpose();

                // fill Z
                startIdx = Z_.rows();
                Z_.conservativeResize(startIdx + 3, 1);
                Eigen::Matrix3d R = state_.getRotation();
                Eigen::Vector3d p = state_.getPosition();
                Eigen::Vector3d d = state_.getVector(it_estimated->second);
                if (state_.getStateType() == StateType::WorldCentric) {
                    Z_.block(startIdx, 0, 3, 1) = R * it.pose.block<3, 1>(0, 3) - (d - p);
                } else {
                    Z_.block(startIdx, 0, 3, 1) = R.transpose() * (it.pose.block<3, 1>(0, 3) - (p - d));
                }
            } else {
                continue;
            }
        }

        // correct state estimates using stacked observations
        if (Z_.rows() > 0) {
            if (state_.getStateType() == StateType::WorldCentric) {
                this->correctRightInvariant(Z_, H_, N_);
            } else {
                this->correctLeftInvariant(Z_, H_, N_);
            }
        }

        // remove contacts from states
        if (!remove_contacts.empty()) {
            auto X_rem_Mat = state_.getXMatrix();
            auto P_rem = state_.getP();

            for (auto it = remove_contacts.begin(); it != remove_contacts.end(); ++it) {
                // remove corresponding row and column from X
                removeMatrixRowAndColumn(it->second, X_rem_Mat, 1);
                // remove 3 corresponding row and column from P
                unsigned int startIdx = 3 + 3 * (it->second - 3);
                removeMatrixRowAndColumn(startIdx, P_rem, 3);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                //removeMatrixRowAndColumn(startIdx, P_rem);
                // update all indices for estimated_landmarks_ and estimated_contact_positions_
                for (auto& estimated_landmark : estimated_landmarks_) {
                    if (estimated_landmark.second > it->second)
                        estimated_landmark.second -= 1;
                }
                for (auto& estimated_contact_position : estimated_contact_positions_) {
                    if (estimated_contact_position.second > it->second)
                        estimated_contact_position.second -= 1;
                }
                // update the indices of remove_contacts in the case where multiple contacts are being removed at once
                for (auto it2 = it; it2 != remove_contacts.end(); ++it2) {
                    if (it2->second > it->second)
                        it2->second -= 1;
                }
                // remove from estimated_contact_positions_
                estimated_contact_positions_.erase(it->first);
            }
            // update state estimate and covariance
            state_.setX(X_rem_Mat);
            state_.setP(P_rem);
        }

        // augment state with newly detected contacts
        if (!new_contacts.empty()) {
            auto X_aug_Mat = state_.getXMatrix();
            auto P_aug = state_.getP();
            for (const auto& new_contact : new_contacts) {
                // initialize new landmark estimated mean
                unsigned int startIdx = X_aug_Mat.rows();
                X_aug_Mat.conservativeResize(startIdx + 1, startIdx + 1);
                X_aug_Mat.block(startIdx, 0, 1, startIdx) = Eigen::MatrixXd::Zero(1, startIdx);
                X_aug_Mat.block(0, startIdx, startIdx, 1) = Eigen::MatrixXd::Zero(startIdx, 1);
                X_aug_Mat(startIdx, startIdx) = 1;
                if (state_.getStateType() == StateType::WorldCentric) {
                    X_aug_Mat.block(0, startIdx, 3, 1) = state_.getPosition() + state_.getRotation() * new_contact.pose.block<3, 1>(0, 3);
                } else {
                    X_aug_Mat.block(0, startIdx, 3, 1) = state_.getPosition() - new_contact.pose.block<3, 1>(0, 3);
                }

                // initialize new landmark covariance
                Eigen::MatrixXd F = Eigen::MatrixXd::Zero(state_.getPDim() + 3, state_.getPDim());
                // for old X
                F.block(0, 0, state_.getPDim() - state_.getThetaDim(), state_.getPDim() - state_.getThetaDim()) = Eigen::MatrixXd::Identity(state_.getPDim() - state_.getThetaDim(), state_.getPDim() - state_.getThetaDim());
                // for theta
                F.block(state_.getPDim() - state_.getThetaDim() + 3, state_.getPDim() - state_.getThetaDim(), state_.getThetaDim(), state_.getThetaDim()) = Eigen::MatrixXd::Identity(state_.getThetaDim(), state_.getThetaDim());
                Eigen::MatrixXd G = Eigen::MatrixXd::Zero(F.rows(), 3);
                // matrix blocks for new contacts
                if ((state_.getStateType() == StateType::WorldCentric && error_type_ == ErrorType::RightInvariant) ||
                    (state_.getStateType() == StateType::BodyCentric && error_type_ == ErrorType::LeftInvariant)) {
                    F.block(state_.getPDim() - state_.getThetaDim(), 6, 3, 3) = Eigen::Matrix3d::Identity();
                    G.block(G.rows() - state_.getThetaDim() - 3, 0, 3, 3) = state_.getWorldRotation();
                } else {
                    F.block(state_.getPDim() - state_.getThetaDim(), 6, 3, 3) = Eigen::Matrix3d::Identity();
                    F.block(state_.getPDim() - state_.getThetaDim(), 0, 3, 3) = skew(-new_contact.pose.block<3, 1>(0, 3));
                    G.block(G.rows() - state_.getThetaDim() - 3, 0, 3, 3) = Eigen::Matrix3d::Identity();
                }
                // augmented P
                P_aug = (F * P_aug * F.transpose() + G * new_contact.covariance.block<3, 3>(3, 3) * G.transpose()).eval();
                // update state and covariance
                state_.setX(X_aug_Mat);
                state_.setP(P_aug);
                // add to estimated_contact_positions_
                estimated_contact_positions_.insert(std::pair<unsigned int, unsigned int>(new_contact.id, startIdx));
            }
        }
    }

    void InEKF::correctLandmarks(const vectorLandmarks& measured_landmarks) {
        // reset and clear observation matrices: Z, H, N
        this->resetObservationMatrices();

        vectorLandmarks new_landmarks;
        std::vector<unsigned int> used_landmark_ids;

        for (const auto& measured_landmark : measured_landmarks) {
            // detect and skip if an ID is not unique
            if (std::find(used_landmark_ids.begin(), used_landmark_ids.end(), measured_landmark.id) != used_landmark_ids.end()) {
                std::cout << "Duplicate landmark ID detected! Skipping measurement... \n";
                continue;
            } else {
                used_landmark_ids.push_back(measured_landmark.id);
            }

            // if id in prior_landmarks_ or estimated_landmarks_
            auto it_prior = prior_landmarks_.find(measured_landmark.id);
            auto it_estimated = estimated_landmarks_.find(measured_landmark.id);
            if (it_prior != prior_landmarks_.end()) {
                // found in prior_landmarks_
                unsigned int dimP = state_.getPDim();
                unsigned int startIdx;

                // fill H
                startIdx = H_.rows();
                H_.conservativeResize(startIdx + 3, dimP);
                H_.block(startIdx, 0, 3, dimP) = Eigen::MatrixXd::Zero(3, dimP);
                if (state_.getStateType() == StateType::WorldCentric) {
                    // skew(p_wl)
                    H_.block(startIdx, 0, 3, 3) = skew(it_prior->second);
                    // -I
                    H_.block(startIdx, 6, 3, 3) = -Eigen::Matrix3d::Identity();
                } else {
                    // -skew(p_wl)
                    H_.block(startIdx, 0, 3, 3) = skew(-it_prior->second);
                    // I
                    H_.block(startIdx, 6, 3, 3) = Eigen::Matrix3d::Identity();
                }

                // fill N
                startIdx = N_.rows();
                N_.conservativeResize(startIdx + 3, startIdx + 3);
                N_.block(startIdx, 0, 3, startIdx) = Eigen::MatrixXd::Zero(3, startIdx);
                N_.block(0, startIdx, startIdx, 3) = Eigen::MatrixXd::Zero(startIdx, 3);
                N_.block(startIdx, startIdx, 3, 3) = state_.getWorldRotation() * measured_landmark.covariance * state_.getWorldRotation().transpose();

                // fill Z
                startIdx = Z_.rows();
                Z_.conservativeResize(startIdx + 3, 1);
                Eigen::Matrix3d R = state_.getRotation();
                Eigen::Vector3d p = state_.getPosition();
                Eigen::Vector3d l = state_.getVector(it_estimated->second);
                if (state_.getStateType() == StateType::WorldCentric) {
                    Z_.block(startIdx, 0, 3, 1) = R * measured_landmark.position - (l - it_prior->second);
                } else {
                    Z_.block(startIdx, 0, 3, 1) = R.transpose() * (measured_landmark.position - (p - it_prior->second));
                }
            } else if (it_estimated != estimated_landmarks_.end()) {
                // found in estimated_landmarks_
                unsigned int dimTheta = state_.getThetaDim();
                unsigned int dimP = state_.getPDim();
                unsigned int startIdx;

                // fill H
                startIdx = H_.rows();
                H_.conservativeResize(startIdx + 3, dimP);
                H_.block(startIdx, 0, 3, dimP) = Eigen::MatrixXd::Zero(3, dimP);
                if (state_.getStateType() == StateType::WorldCentric) {
                    // -I
                    H_.block(startIdx, 6, 3, 3) = -Eigen::Matrix3d::Identity();
                    // I
                    H_.block(startIdx, 3 * it_estimated->second - dimTheta, 3, 3) = Eigen::Matrix3d::Identity();
                } else {
                    // I
                    H_.block(startIdx, 6, 3, 3) = Eigen::Matrix3d::Identity();
                    // -I
                    H_.block(startIdx, 3 * it_estimated->second - dimTheta, 3, 3) = -Eigen::Matrix3d::Identity();
                }

                // fill N
                startIdx = N_.rows();
                N_.conservativeResize(startIdx + 3, startIdx + 3);
                N_.block(startIdx, 0, 3, startIdx) = Eigen::MatrixXd::Zero(3, startIdx);
                N_.block(0, startIdx, startIdx, 3) = Eigen::MatrixXd::Zero(startIdx, 3);
                N_.block(startIdx, startIdx, 3, 3) = state_.getWorldRotation() * measured_landmark.covariance * state_.getWorldRotation().transpose();

                // fill Z
                startIdx = Z_.rows();
                Z_.conservativeResize(startIdx + 3, 1);
                Eigen::Matrix3d R = state_.getRotation();
                Eigen::Vector3d p = state_.getPosition();
                Eigen::Vector3d l = state_.getVector(it_estimated->second);
                if (state_.getStateType() == StateType::WorldCentric) {
                    Z_.block(startIdx, 0, 3, 1) = R * measured_landmark.position - (l - p);
                } else {
                    Z_.block(startIdx, 0, 3, 1) = R.transpose() * (measured_landmark.position - (p - l));
                }
            } else {
                // first time landmark as been detected (add to list for later state augmentation)
                new_landmarks.push_back(measured_landmark);
            }
        }

        // correct state using stacked observations
        if (Z_.rows() > 0) {
            if (state_.getStateType() == StateType::WorldCentric) {
                this->correctRightInvariant(Z_, H_, N_);
            } else {
                this->correctLeftInvariant(Z_, H_, N_);
            }
        }

        // augment state with newly detected landmarks
        if (!new_landmarks.empty()) {
            auto X_aug_Mat = state_.getXMatrix();
            auto P_aug = state_.getP();
            for (const auto& new_landmark : new_landmarks) {
                // initialize new landmark estimated mean
                unsigned int startIdx = X_aug_Mat.rows();
                X_aug_Mat.conservativeResize(startIdx + 1, startIdx + 1);
                X_aug_Mat.block(startIdx, 0, 1, startIdx) = Eigen::MatrixXd::Zero(1, startIdx);
                X_aug_Mat.block(0, startIdx, startIdx, 1) = Eigen::MatrixXd::Zero(startIdx, 1);
                X_aug_Mat(startIdx, startIdx) = 1;
                X_aug_Mat.block(0, startIdx, 3, 1) = state_.getPosition() + state_.getRotation() * new_landmark.position;

                // initialize new landmark covariance
                Eigen::MatrixXd F = Eigen::MatrixXd::Zero(state_.getPDim() + 3, state_.getPDim());
                // for old X
                F.block(0, 0, state_.getPDim() - state_.getThetaDim(), state_.getPDim() - state_.getThetaDim()) = Eigen::MatrixXd::Identity(state_.getPDim() - state_.getThetaDim(), state_.getPDim() - state_.getThetaDim());
                // for theta
                F.block(state_.getPDim() - state_.getThetaDim() + 3, state_.getPDim() - state_.getThetaDim(), state_.getThetaDim(), state_.getThetaDim()) = Eigen::MatrixXd::Identity(state_.getThetaDim(), state_.getThetaDim());
                Eigen::MatrixXd G = Eigen::MatrixXd::Zero(F.rows(), 3);
                // matrix blocks for new landmarks
                if (error_type_ == ErrorType::RightInvariant) {
                    F.block(state_.getPDim() - state_.getThetaDim(), 6, 3, 3) = Eigen::Matrix3d::Identity();
                    G.block(G.rows() - state_.getThetaDim() - 3, 0, 3, 3) = state_.getRotation();
                } else {
                    F.block(state_.getPDim() - state_.getThetaDim(), 6, 3, 3) = Eigen::Matrix3d::Identity();
                    F.block(state_.getPDim() - state_.getThetaDim(), 0, 3, 3) = skew(-new_landmark.position);
                    G.block(G.rows() - state_.getThetaDim() - 3, 0, 3, 3) = Eigen::Matrix3d::Identity();
                }
                // augmented P
                P_aug = (F * P_aug * F.transpose() + G * new_landmark.covariance * G.transpose()).eval();
                // update state and covariance
                state_.setX(X_aug_Mat);
                state_.setP(P_aug);
                // add to estimated_landmarks_
                estimated_landmarks_.insert(std::pair<unsigned int, unsigned int>(new_landmark.id, startIdx));
            }
        }
    }

    void InEKF::correctContactPosition(const unsigned int id, const Eigen::Vector3d& measured_contact_position, const Eigen::Matrix3d& covariance, const Eigen::Vector3d& indices) {
        // reset and clear observation matrices: Z, H, N
        this->resetObservationMatrices();

        // mapping matrix for reducing matrix dimensions
        Eigen::MatrixXd PI;
        // matrices for reduced dimensions
        Eigen::MatrixXd Z, H, N;
        // if find id in the estimated_contact_positions_
        auto it_estimated = estimated_contact_positions_.find(id);
        if (it_estimated != estimated_contact_positions_.end()) {
            // fill PI
            unsigned int startIdx;
            for (int i = 0; i < 3; ++i) {
                if (indices(i) != 0) {
                    startIdx = PI.rows();
                    PI.conservativeResize(startIdx + 1, 3);
                    PI.block(startIdx, 0, 1, 3) = Eigen::MatrixXd::Zero(1, 3);
                    PI(startIdx, i) = 1;
                }
            }
            if (PI.rows() == 0) {
                return;
            }

            // fill observation matrices
            auto dimP = state_.getPDim();

            // get contact position
            Eigen::Vector3d d = state_.getVector(it_estimated->second);

            // fill H
            H_ = Eigen::MatrixXd::Zero(3, dimP);
            H_.block<3, 3>(0, 0) = -skew(d);
            H_.block<3, 3>(0, 3 * it_estimated->second - 6) = Eigen::Matrix3d::Identity();
            H = PI * H_;
            // fill N
            N_ = covariance;
            N = PI * N_ * PI.transpose();
            // fill Z
            Z_ = measured_contact_position - d;
            Z = PI * Z_;
            // correct
            this->correctRightInvariant(Z, H, N);
        }
    }

}
