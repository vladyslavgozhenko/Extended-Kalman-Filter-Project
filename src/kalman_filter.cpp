#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predicts the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * updates the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * updates the state by using Extended Kalman Filter equations
  */
  //recovers state parameters
	  float px = x_(0);
	  float py = x_(1);
	  float vx = x_(2);
	  float vy = x_(3);

	  // checks for 0 like values
	  if (fabs(px) < 0.00001 && fabs(py) < 0.00001) {
			px = 0.00001;
			py = 0.00001;
	  }
	  else if (fabs(px) < 0.00001) {
			px = 0.00001;
	  }

	  // accounts for polar data (change variable names)
	  float temp_1 = sqrt(px*px + py*py);
	  float temp_2 = atan2(py, px);
	  float temp_3 = (px*vx + py*vy) / temp_1;

	  // creates vector for z_pred
	  VectorXd hx(3, 1);

	  // assigns values to z_pred
	  hx << temp_1, temp_2, temp_3;

	  // computes kalman update as normal
	  VectorXd y = z - hx;

	  // normalizes y(1)
	  float pi = 3.14159;

	  while (y(1) < -pi) {
			y(1) += 2 * pi;
	  }

	  while (y(1) > pi) {
			y(1) -= 2 * pi;
	  }

	  MatrixXd Ht = H_.transpose();
	  MatrixXd PHt = P_ * Ht;
	  MatrixXd S = H_ * PHt + R_;
	  MatrixXd Si = S.inverse();
	  MatrixXd K = PHt * Si;

	  //new estimate
	  x_ = x_ + (K * y);
	  long x_size = x_.size();
	  MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
