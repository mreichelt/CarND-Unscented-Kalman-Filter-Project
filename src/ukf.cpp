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
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;

    // Process noise standard deviation yaw acceleration in rad/s^2
    // TODO too high!
    std_yawdd_ = 30;

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
    weights_ = VectorXd(n_sigma_);
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
    time_us_ = meas_package.timestamp_;

    Prediction(delta_t);
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
}

void UKF::Initialize(MeasurementPackage &meas_package) {
    VectorXd &raw = meas_package.raw_measurements_;

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
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */


    // generate sigma points (augmented due to noise)
    MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();

    // TODO: predict sigma points


    // TODO: predict mean & covariance
}

MatrixXd UKF::GenerateAugmentedSigmaPoints() {
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

    // create augmented mean state
    x_aug.head(n_x_) = x_;

    //create augmented covariance matrix
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


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
}

