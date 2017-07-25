#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;
    n_sigma_ = 2 * n_aug_ + 1;

    // initial state vector
    x_ = VectorXd::Zero(n_x_);

    // initial covariance matrix
    P_ = MatrixXd::Zero(n_x_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;

    // Process noise standard deviation yaw acceleration in rad/s^2
    // TODO adjust using NIS
    std_yawdd_ = 1;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    // init state covariance matrix along the diagonal
    // first two values are what we expect from x and y positions, other we leave as ones (like identity matrix)
    P_.diagonal() << std_laspx_ * std_laspx_, std_laspy_ * std_laspy_, 1, 1, 1;

    // initialize sigma points matrix
    Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_);

    // initialize weights
    weights_ = VectorXd::Zero(n_sigma_);
    weights_[0] = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < n_sigma_; i++) {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights_[i] = weight;
    }

    is_initialized_ = false;
    time_us_ = 0;
}

UKF::~UKF() = default;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    if (!is_initialized_) {
        Initialize(meas_package);
        is_initialized_ = true;
    }

    // compute time delta in seconds
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
        time_us_ = meas_package.timestamp_;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        UpdateLidar(meas_package);
        time_us_ = meas_package.timestamp_;
    }
}

void UKF::Initialize(MeasurementPackage &meas_package) {
    VectorXd raw = meas_package.raw_measurements_;

    // initialize time
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        x_ << raw[0], raw[1], 0, 0, 0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // map coordinates from polar to cartesian
        double ro = raw[0],
                theta = raw[1],
                ro_dot = raw[2];

        // convert polar coordinates from radar to cartesian
        double px = ro * cos(theta),
                py = ro * sin(theta),
                vx = ro_dot * cos(theta),
                vy = ro_dot * sin(theta),
                v = sqrt(vx * vx + vy * vy);

        // use directional speed v instead of vx & vy
        x_ << px, py, v, 0, 0;

    } else {
        cerr << "unknown measurement type" << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // generate sigma points (augmented due to noise)
    MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();

    // predict sigma points: just apply model function
    PredictSigmaPoints(Xsig_aug, delta_t);

    // predict mean & covariance, store into x_ and P_
    PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateAugmentedSigmaPoints() {
    VectorXd x_aug = VectorXd::Zero(n_aug_);
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sigma_);

    // create augmented mean state
    x_aug.head(n_x_) = x_;

    // create augmented covariance matrix
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug.bottomRightCorner(2, 2) << std_a_ * std_a_, 0,
            0, std_yawdd_ * std_yawdd_;

    // create square root matrix
    MatrixXd A_aug = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug.block(0, 1, n_aug_, n_aug_) = sqrt(3.0) * A_aug;
    Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = -sqrt(3.0) * A_aug;
    for (int col = 0; col < n_sigma_; col++) {
        Xsig_aug.col(col) += x_aug;
    }

    return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd &Xsig_aug, double delta_t) {
    Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_);

    // predict sigma points
    for (int i = 0; i < n_sigma_; i++) {
        // extract values for better readability
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        // add noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        // write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }
}

void UKF::PredictMeanAndCovariance() {
    VectorXd x = VectorXd::Zero(n_x_);
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_);

    // predict state
    for (int i = 0; i < n_sigma_; i++) {
        x += weights_[i] * Xsig_pred_.col(i);
    }
    x_ = x;

    // predict state covariance matrix
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        // normalize rho
        x_diff[3] = tools.normalizeAngle(x_diff[3]);
        P += weights_[i] * x_diff * x_diff.transpose();
    }
    P_ = P;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // px, py
    int n_z = 2;

    // transform Xsig_pred_ into measurement space - simply take the first 2 rows
    MatrixXd Zsig = Xsig_pred_.topRows(2);

    // calculate mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z);
    for (int i = 0; i < n_sigma_; i++) {
        z_pred += weights_[i] * Zsig.col(i);
    }

    // calculate measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);
    MatrixXd R = MatrixXd::Zero(n_z, n_z);
    R << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_[i] * z_diff * z_diff.transpose();
    }
    S += R;

    MatrixXd Tc = buildCrossCorrelationMatrix(n_z, Zsig, z_pred);

    // kalman gain
    MatrixXd K = Tc * S.inverse();

    // the actual measurement
    VectorXd z = meas_package.raw_measurements_;

    // residual
    VectorXd z_diff = z - z_pred;

    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // r, phi, r_dot
    int n_z = 3;

    // transform sigma points into measurement space
    MatrixXd Zsig = MatrixXd::Zero(n_z, n_sigma_);
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd x = Xsig_pred_.col(i);
        double px = x[0],
                py = x[1],
                v = x[2],
                phi = x[3],
                r = sqrt(px * px + py * py);

        // avoid division by zero
        if (fabs(r) < 0.001) {
            r = 0.001;
        }

        Zsig.col(i) << r,
                atan2(py, px),
                v * (px * cos(phi) + py * sin(phi)) / r;
    }


    // calculate mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z);
    for (int i = 0; i < n_sigma_; i++) {
        z_pred += weights_[i] * Zsig.col(i);
    }

    // calculate measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z, n_z);
    MatrixXd R = MatrixXd::Zero(n_z, n_z);
    R << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_[i] * z_diff * z_diff.transpose();
    }
    S += R;


    MatrixXd Tc = buildCrossCorrelationMatrix(n_z, Zsig, z_pred);

    // kalman gain
    MatrixXd K = Tc * S.inverse();

    // the actual measurement
    VectorXd z = meas_package.raw_measurements_;

    // residual
    VectorXd z_diff = z - z_pred;
    z_diff[1] = tools.normalizeAngle(z_diff[1]);

    // update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}

MatrixXd UKF::buildCrossCorrelationMatrix(int n_z, MatrixXd &Zsig, VectorXd &z_pred) {
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
    for (int i = 0; i < n_sigma_; i++) {

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff[1] = tools.normalizeAngle(z_diff[1]);

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff[3] = tools.normalizeAngle(x_diff[3]);

        Tc = Tc + weights_[i] * x_diff * z_diff.transpose();
    }
    return Tc;
}
