// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 0.05;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 2;
constexpr double INIT_PSI_STD = 5.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

// ----------------------------------------------------------------------- //
// USEFUL HELPER FUNCTIONS
double getLambda(int numStates)
{
    return (3.0 - numStates);
}

VectorXd normaliseState(VectorXd state)
{
    state(2) = wrapAngle(state(2));
    return state;
}
VectorXd normaliseLidarMeasurement(VectorXd meas)
{
    meas(1) = wrapAngle(meas(1));
    return meas;
}
std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov)
{
    std::vector<VectorXd> sigmaPoints;

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE
    int numStates = state.size();
    int numPoints = 2*numStates + 1;
    int lambda = getLambda(numStates);
    MatrixXd sqrtCov = cov.llt().matrixL();

    // x(0)
    sigmaPoints.push_back(state);

    // x(1->n)
    for (int i = 0; i < numStates; i++) {
        sigmaPoints.push_back(normaliseState(state + sqrt(numStates+lambda)*sqrtCov.col(i)));
    }

    // x(n+1->2n)
    for (int i = 0; i < numStates; i++) {
        sigmaPoints.push_back(normaliseState(state - sqrt(numStates+lambda)*sqrtCov.col(i)));
    }
    // ----------------------------------------------------------------------- //

    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates)
{
    std::vector<double> weights;

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE
    double lambda = getLambda(numStates);
    std::cout << "numStates: " << numStates << std::endl;
    std::cout << "lambda: " << lambda << std::endl;
    weights.push_back(lambda / (lambda + numStates));
    for (int i = 1; i <= 2*numStates; i++) {
        weights.push_back(1.0/(2.0*(lambda + numStates)));
    }
    // ----------------------------------------------------------------------- //

    return weights;
}

VectorXd lidarMeasurementModel(VectorXd augState, double beaconX, double beaconY)
{
    VectorXd z_hat = VectorXd::Zero(2);

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE

    // ----------------------------------------------------------------------- //

    return z_hat;
}

VectorXd vehicleProcessModel(VectorXd augState, double psi_dot, double dt)
{
    VectorXd new_state = VectorXd::Zero(4);

    // ----------------------------------------------------------------------- //
    // ENTER YOUR CODE HERE
    VectorXd state_delta = VectorXd(4);
    VectorXd state = augState.head(4);
    double omega_psi_dot = augState(4);
    double omega_a = augState(5);
    double v = state(3);
    double psi = state(2);
    state_delta << dt*v*cos(psi), dt*v*sin(psi), dt*(psi_dot+omega_psi_dot), dt*omega_a;

    new_state = normaliseState(state + state_delta);
    // ----------------------------------------------------------------------- //

    return new_state;
}
// ----------------------------------------------------------------------- //

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the Lidar Measurements in the 
        // section below.
        // HINT: Use the normaliseState() and normaliseLidarMeasurement() functions
        // to always keep angle values within correct range.
        // HINT: Do not normalise during sigma point calculation!
        // HINT: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // HINT: The mapped-matched beacon position can be accessed by the variables
        // map_beacon.x and map_beacon.y
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1) // Check that we have a valid beacon match
        {
           

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
        // HINT: Use the normaliseState() function to always keep angle values within correct range.
        // HINT: Do NOT normalise during sigma point calculation!
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
        int numStates = state.size();
        int numNoise = 2;
        int numAugStates = numStates + numNoise;

        // make augmented state, then generate sigma points and sigma weights
        VectorXd augState = VectorXd::Zero(numAugStates);
        MatrixXd augCov = MatrixXd::Zero(numAugStates,numAugStates);
        augState.head(numStates) = state;
        augCov.topLeftCorner(numStates,numStates) = cov;
        augCov(4,4) = GYRO_STD*GYRO_STD;
        augCov(5,5) = ACCEL_STD*ACCEL_STD;
        std::vector<double> weights = generateSigmaWeights(numAugStates);
        std::vector<VectorXd> sigmaPoints = generateSigmaPoints(augState,augCov);
        std::cout << "sigma weights" << std::endl;
        for (const auto& element : weights) {
            std::cout << element << " ";
        }
        std::cout << std::endl;

        // translate sigma points
        std::vector<VectorXd> translatedSigmaPoints;
        for (int i = 0; i < sigmaPoints.size(); i++) {
            VectorXd thisTranslatedPoint = vehicleProcessModel(sigmaPoints[i],gyro.psi_dot,dt);
            translatedSigmaPoints.push_back(thisTranslatedPoint);
            if (i == 0) {
                std::cout << "sigma point" << std::endl;
                std::cout << sigmaPoints[i] << std::endl;
                std::cout << "Translated sigma point" << std::endl;
                std::cout << thisTranslatedPoint << std::endl;
            }
            
        }

        // generate a priori state estimate
        state = VectorXd::Zero(numStates);
        for (int i = 0; i < weights.size(); i++) {
            state = normaliseState(state + weights[i]*translatedSigmaPoints[i]);
        }

        std::cout << "new state: " << std::endl;
        std::cout << state << std::endl;

        // generate a priori covariance matrix
        cov = MatrixXd::Zero(numStates,numStates);
        for (int i = 0; i < weights.size(); i++) {
            VectorXd stateDiff = normaliseState(translatedSigmaPoints[i]-state);
            cov = cov + weights[i]*stateDiff*stateDiff.transpose();
        }


        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the UKF update state would just produce the same result.
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
        cov = (MatrixXd::Identity(4,4) - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // You may modify this initialisation routine if you can think of a more
        // robust and accuracy way of initialising the filter.
        // ----------------------------------------------------------------------- //
        // YOU ARE FREE TO MODIFY THE FOLLOWING CODE HERE

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

        // ----------------------------------------------------------------------- //
    }             
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
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
