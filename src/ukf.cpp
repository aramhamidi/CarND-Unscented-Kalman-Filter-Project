#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <cmath>

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

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    // P should remain symmetr
    //you can try setting the diagonal values by how much difference you expect between the true
    //state and the initialized x state vector, by considering the standard deviation of the lidar x and y
    //measurements.
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.2;//30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.2;//30;

    //--------measurement noise values are determined by the sensor manufacturer and should not be changed.
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

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;
    is_initialized_ = false;
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2 * n_aug_ + 1);

    time_us_ = 0.0;
    NIS_RADAR_ = 0.0;
    NIS_LIDAR_ = 0.0;
}
/*
UKF::UKF(const VectorXd &weights_){
	weights_(0) = lambda_ / (n_aug_ + lambda_);
	for (int i = 1; i < 2*n_aug_ + 1; i++){
		weights_(i) = 0.5 / (n_aug_ + lambda_);
	}
}
*/

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if(!is_initialized_) {
	cout << "Initializing the UKF: "<< endl;
        // At the start up you dont know  where the bycicle is, but after the first measurement
        // you can initialize px and py.
        // for the other variables in X try different values to see what works the best.
        // x = [px,py,v,yaw,yaw_dot]
	x_ << 1, 1, 1, 1, 1;
	P_ << 0.15,0,0,0,0,
	   0,0.15,0,0,0,
	   0,0,0.15,0,0,
	   0,0,0,0.15,0,
	   0,0,0,0,0.15;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            // TODO: check and make sure that phi is measured in radians?
            // VectorXd x(3);
	    x_[0] = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
	    x_[1] = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
            time_us_ = meas_package.timestamp_;
            is_initialized_ = true;

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
	      Initialize state.
             */
            //measurement matrix
            x_[0] = meas_package.raw_measurements_[0];
	    x_[1] = meas_package.raw_measurements_[1];
            time_us_ = meas_package.timestamp_;
            is_initialized_ = true;
        }
	

        // done initializing, no need to predict or update
        //is_initialized_ = true;
        return;
    }

    //float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    std::cout << "Prediction ..." << std::endl;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    std::cout << "delta_t: " << delta_t << endl;
    //std::cout << "Time_Stamp: "<< time_us_ << std::endl;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
	    std::cout << "Update Lidar ..." << endl;
	    UpdateLidar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
	    cout << "Update Radar ..." << endl;
	    UpdateRadar(meas_package);
    }

   // print NIS
   cout << "NIS_RADAR_ = " << NIS_RADAR_  << endl;
   cout << "NIS_LIDAR_ = " << NIS_LIDAR_  << endl;
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(float delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
    Xsig.fill(0.0);
    GenerateSigmaPoints(&Xsig);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0.0);
    AugmentedSigmaPoints(&Xsig_aug);
    Xsig_pred_.fill(0.0);
    SigmaPointPrediction(&Xsig_pred_, Xsig_aug, delta_t);
    VectorXd Xpred = VectorXd(n_x_);
    Xpred.fill(0.0);
    MatrixXd Ppred = MatrixXd(n_x_, n_x_);
    Ppred.fill(0.0);
    PredictMeanAndCovariance(&Xpred, &Ppred, Xsig_pred_);
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

  You'll also need to calculate the lidar NIS.(Normalized Innovation Squared)
  */
    int n_z = 2;
    /*VectorXd z = VectorXd(n_z);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
	    z = meas_package.raw_measurements_ ;
    }
     */
    //mean predicted measurement
    VectorXd Zpred = VectorXd(n_z);
    Zpred.fill(0.0);
    //measurement covariance matrix S
    MatrixXd Spred = MatrixXd(n_z, n_z);
    Spred.fill(0.0);
    //create matrix for sigma points in measured space
    MatrixXd Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig_.fill(0.0);
    PredictLidarMeasurement(&Zpred, &Spred, &Zsig_);
    //UpdateState(n_z, Zpred, Spred, z , Zsig_, meas_package);
    UpdateState(n_z, Zpred, Spred, Zsig_, meas_package);
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

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    /*VectorXd z = VectorXd(n_z);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
	    z = meas_package.raw_measurements_ ;
    }
     */
    //mean predicted measurement
    VectorXd Zpred = VectorXd(n_z);
    Zpred.fill(0.0);
    //measurement covariance matrix S
    MatrixXd Spred = MatrixXd(n_z, n_z);
    Spred.fill(0.0);
    //create matrix for sigma points in measurement space
    MatrixXd Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig_.fill(0.0);
    PredictRadarMeasurement(&Zpred,&Spred, &Zsig_);
    //UpdateState(n_z, Zpred, Spred, z, Zsig_, meas_package);
    UpdateState(n_z, Zpred, Spred, Zsig_, meas_package);

}

// Additional Functions from lecture
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
    //create sigma point matrix
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
    Xsig.fill(0.0);
    //calculate square root of P
    MatrixXd A = P_.llt().matrixL();
    //set first column of sigma point matrix
    Xsig.col(0)  = x_;
    //set remaining sigma points
    for (int i = 0; i < n_x_; i++)
    {
        Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
        Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
    }
    //print result
    //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

    //write result
    *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.fill(0.0);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0.0);

    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }

    //print result
    //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

    //write result
    *Xsig_out = Xsig_aug;
}
void UKF::SigmaPointPrediction(MatrixXd* Xsig_pred_, const MatrixXd Xsig_out, float delta_t) {

    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred.fill(0.0);
    //double delta_t = 0.1; //time diff in sec

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_out(0,i);
        double p_y = Xsig_out(1,i);
        double v = Xsig_out(2,i);
        double yaw = Xsig_out(3,i);
        double yawd = Xsig_out(4,i);
        double nu_a = Xsig_out(5,i);
        double nu_yawdd = Xsig_out(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+ yawd * delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }

    //print result
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

    //write result
    *Xsig_pred_ = Xsig_pred;

}
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, const MatrixXd Xsig_pred_) {
    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);

    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

    //predicted state mean
    x.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x = x+ weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P = P + weights_(i) * x_diff * x_diff.transpose() ;
    }

    //print result
    std::cout << "Predicted state" << std::endl;
    std::cout << x << std::endl;
    std::cout << "Predicted covariance matrix" << std::endl;
    std::cout << P << std::endl;

    //write result
    *x_out = x;
    *P_out = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_) {

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    *Zsig_ = Zsig;
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;

    //print result
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

    //write result
    *z_out = z_pred;
    *S_out = S;
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_) {

    //set measurement dimension, Lidar can measure px and py
    int n_z = 2;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);

        // measurement model
        Zsig(0,i) = p_x;                        //px
        Zsig(1,i) = p_y;                                 //py

    }
    *Zsig_ = Zsig;
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R << std_laspx_*std_laspy_, 0,
    0, std_laspy_*std_laspy_;

    S = S + R;

    //print result
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

    //write result
    *z_out = z_pred;
    *S_out = S;
}


//void UKF::UpdateState(int n_z, const VectorXd z_out, const MatrixXd S_out, const VectorXd z, const MatrixXd Zsig_, MeasurementPackage meas_package){
void UKF::UpdateState(int n_z, const VectorXd z_out, const MatrixXd S_out, const MatrixXd Zsig_, MeasurementPackage meas_package){
    
    // Get the new measurements
    VectorXd z = VectorXd(n_z);
    z = meas_package.raw_measurements_ ;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i =0; i < 2*n_aug_ +1 ; i++ ){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd z_diff = Zsig_.col(i) - z_out;
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
        //angle normalization
            while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
            while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
            while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
            while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        }
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    //calculate Kalman gain K;
    MatrixXd K = Tc * (S_out.inverse());

    //update state mean and covariance matrix
    //residual
    VectorXd z_diff = z - z_out;

    //angle normalization
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){

        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    }
    //Calculate NIS
    double NIS = z_diff.transpose() * S_out.inverse() * z_diff;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	    NIS_LIDAR_ = NIS;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
	    NIS_RADAR_ = NIS;
    }
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * (S_out) * K.transpose();

    //print result
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

    //write result
    //*x_out = x_;
    //*P_out = P_;
}
