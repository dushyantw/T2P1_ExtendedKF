#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.0001 // A very small number
/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  R_radar_ = MatrixXd(3, 3);
  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  // Initializing the H Matrix
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  Hj_ = MatrixXd(3, 4);
  // Initializing the Jacobian Matrix
  Hj_ << 1,1,0,0,
         1,1,0,0,
         1,1,1,1;

  ekf_.F_ = MatrixXd(4,4);
  // Initializing the transition F Matrix
  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  // adding Q matrix
  // prediction uncertainty covariance
    ekf_.Q_ = MatrixXd(4,4);
    ekf_.Q_ << 0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    // ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        float rho = measurement_pack.raw_measurements_(0); // range
        float phi = measurement_pack.raw_measurements_(1); // bearing
        float rho_dot = measurement_pack.raw_measurements_(2); // velocity
        float x = rho * cos(phi);
        float y = rho * sin(phi);
        float vx = rho_dot * cos(phi);
        float vy = rho_dot * sin(phi);
        ekf_.x_ << x, y, vx, vy;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    ekf_.x_ << measurement_pack.raw_measurements_(0,0), measurement_pack.raw_measurements_(1,0), 0.0, 0.0;
    }

    // Deal with the special case initialisation problems
    if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
        ekf_.x_(0) = EPS;
        ekf_.x_(1) = EPS;
    }

    // Initial covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    // Print the initialization results
    cout << "EKF init: " << ekf_.x_ << endl;

    // Save the initiall timestamp for dt calculation
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    if(ekf_.x_[0]!=0 || ekf_.x_[1]!=0){
        is_initialized_ = true;
    }

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1.0e6; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // update the process noice covariance matrix
  float noise_ax = 9;
  float noise_ay = 9;
  float dt_square_by_2 = (dt * dt) / 2;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
           0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
           dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
           0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

   ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
          Hj_ = tools.CalculateJacobian(ekf_.x_);
          ekf_.H_ = Hj_;
          ekf_.R_ = R_radar_;
          ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
