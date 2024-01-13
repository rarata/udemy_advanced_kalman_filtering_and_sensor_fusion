// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Extended Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
}

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the Lidar Measurements in the 
        // section below.
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // HINT: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // HINT: The mapped-matched beacon position can be accessed by the variables
        // map_beacon.x and map_beacon.y
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1)
        {           
            // The map matched beacon positions can be accessed using: map_beacon.x AND map_beacon.y
            Vector2d z = Vector2d(); // z = measurement (in eigen matrix format)
            Vector2d z_hat = Vector2d(); // z_hat = expected measurement based on current state estimate
            Vector2d y = Vector2d(); // y = measurement innovation
            MatrixXd H = MatrixXd(2,4); // H = measurement jacobian matrix
            MatrixXd R = Matrix2d(); // R = measurment covariance matrix
            MatrixXd S = Matrix2d(); // S = innovation covariance matrix
            MatrixXd K = MatrixXd(4,2); // K = kalman filter gain matrix

            double dx_hat = map_beacon.x - state(0);
            double dy_hat = map_beacon.y - state(1);
            double z_hat_rng = sqrt(dx_hat*dx_hat + dy_hat*dy_hat);
            double z_hat_ang = wrapAngle(atan2(dy_hat,dx_hat)-state(2));
            z << meas.range,meas.theta;
            z_hat << z_hat_rng, z_hat_ang;
            y = z - z_hat;
            y(1) = wrapAngle(y(1));
            H << -dx_hat/z_hat_rng,-dy_hat/z_hat_rng,0,0, dy_hat/(z_hat_rng*z_hat_rng),-dx_hat/(z_hat_rng*z_hat_rng),-1,0;
            R << LIDAR_RANGE_STD*LIDAR_RANGE_STD,0, 0,LIDAR_THETA_STD*LIDAR_THETA_STD;
            S = H*cov*H.transpose() + R;
            K = cov*H.transpose()*S.inverse();

            state = state + K*y;
            state(2) = wrapAngle(state(2));
            cov = (Matrix4d::Identity() - K*H)*cov;

        }

        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // HINT: Assume the state vector has the form [PX, PY, PSI, V].
        // HINT: Use the Gyroscope measurement as an input into the prediction step.
        // HINT: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
        MatrixXd F = Matrix4d(); // F = Jacobian State matrix
        MatrixXd Q = Matrix4d(); // Q = process noise
        VectorXd delta_state = Vector4d();
        double v = state(3);
        double psi = state(2);

        delta_state << dt*v*cos(psi), dt*v*sin(psi), dt*gyro.psi_dot, 0;
        state = state + delta_state;
        state(2) = wrapAngle(state(2));

        F << 1,0,-dt*v*sin(psi),dt*cos(psi), 0,1,dt*v*cos(psi),dt*sin(psi), 0,0,1,0, 0,0,0,1;
        Q << 0,0,0,0, 0,0,0,0, 0,0,dt*dt*GYRO_STD*GYRO_STD,0, 0,0,0,dt*dt*ACCEL_STD*ACCEL_STD;

        cov = F*cov*F.transpose() + Q;
        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the EKF update state would just produce the same result.
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,4);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,1,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (Matrix4d::Identity() - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

        setState(state);
        setCovariance(cov);
    } 
             
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
