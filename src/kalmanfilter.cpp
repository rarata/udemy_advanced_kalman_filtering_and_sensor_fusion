// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Extended Kalman Filter
// Capstone project
// Ryan Arata
//
// Problem:
// Build a kalman filter to handle the following with the best possible performance:
// Vehicle travels on a 2D plane.  Turning control is provided by a gyro.  Acceleration is not
// a known control input or measured, so it must be assuemd to be a random variable.
// Position measurements are provided by GPS, but there may be erroneous readings.
// Lidar provides range and angle measurements from beacons of known location.
// 
// Problem elements:
// Gyroscope input may have bias and bias drift
// Faulty GPS sensor measurements
// Areas of the simulation are GPS-denied
// Lidar beacon data association
// Initialization of the kalman filter

// Task list
// 1. lidar data association
// 2. initialization
// 3. handle faultly GPS data
// 4. handle gyro bias
// 5. if necessary, handle gps-denied locations
// 6. mess with consants to improve performance

#include "kalmanfilter.h"
#include "utils.h"
#include <limits>
#include <algorithm>

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double INIT_GYRO_BIAS_STD = 0.01/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 4.0;
constexpr double LIDAR_THETA_STD = 0.02;
constexpr double LIDAR_MAX_RANGE = 90.0;
constexpr int NUM_STATE = 5; // x,y,psi,v,gyro_bias
constexpr int MIN_INIT_SAMPLES = 5;
// -------------------------------------------------- //
// --------------Helper functions-------------------- //
double angleDelta(double theta1, double theta2) {
    // returns the delta between the angles in the range [-pi,pi)
    double delta = theta1 - theta2;
    if (delta > M_PI) {
        delta = delta - 2*M_PI;
    }
    else if (delta <= -M_PI) {
        delta = delta + 2*M_PI;
    }
    return delta;
}

double calculateAverageAngle(const std::vector<double>& angles) {
    if (angles.empty()) {
        // Handle the case where the vector is empty
        std::cerr << "Error: Vector is empty." << std::endl;
        return 0.0; // Return a default angle
    }

    // Convert angles to Cartesian coordinates
    std::vector<double> xCoordinates;
    std::vector<double> yCoordinates;

    for (double angle : angles) {
        xCoordinates.push_back(cos(angle));
        yCoordinates.push_back(sin(angle));
    }

    // Calculate the average Cartesian coordinates
    double avgX = std::accumulate(xCoordinates.begin(), xCoordinates.end(), 0.0) / xCoordinates.size();
    double avgY = std::accumulate(yCoordinates.begin(), yCoordinates.end(), 0.0) / yCoordinates.size();

    // Calculate the average angle from the average Cartesian coordinates
    double avgAngle = atan2(avgY, avgX);

    // Normalize the result to the range [0, 2*pi)
    avgAngle = wrapAngle(avgAngle);

    return avgAngle;
}

VectorXd normaliseState(VectorXd state) {
    state(2) = wrapAngle(state(2));
    return state;
}

VectorXd normaliseLidarMeasurement(VectorXd meas) {
    meas(1) = wrapAngle(meas(1));
    return meas;
}

double getLidarDistanceError(VectorXd state, LidarMeasurement meas, BeaconData beacon) {
    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double x_hat = x + meas.range*cos(meas.theta + psi);
    double y_hat = y + meas.range*sin(meas.theta + psi);
    double distance = sqrt((beacon.x - x_hat)*(beacon.x - x_hat) + (beacon.y - y_hat)*(beacon.y - y_hat));
    return distance;
}

int associateBeaconID(VectorXd state, LidarMeasurement meas, const BeaconMap& map) {
    double x = state[0];
    double y = state[1];
    
    // get all beacons possibly within range
    double range = LIDAR_MAX_RANGE + 6*GPS_POS_STD; 
    std::vector<BeaconData> beacons = map.getBeaconsWithinRange(x, y, range);

    // find the beacon closest to the location the measurement predicts
    double minDistance = std::numeric_limits<double>::max();
    int closestBeaconID = -1;

    for (const BeaconData& beacon : beacons) {
        // calculate distance between expected location and this beacon
        double distance = getLidarDistanceError(state, meas, beacon);

        // Update closest beacon if the current distance is smaller
        if (distance < minDistance) {
            minDistance = distance;
            closestBeaconID = beacon.id;
        }
    }

    return closestBeaconID;
}

double estimatePsiFromLidar(double x, double y, std::vector<LidarMeasurement> dataset, const BeaconMap& map) {
    // provides a best estimate value for psi, the angle of the vehicle, given a list of lidar measurements and estimated vehicle location
    //
    // Method:
    // Using the first lidar measurement, find all beacons it could map to
    // Test the psi angles that the measurement would suggest for each of those beacons by assuming that
    // orientation and mapping all the lidar measurements against beacons.  Select the psi value that has
    // the lowest RMSE across the whole dataset mapping.
    double firstMeasRange = dataset[0].range;
    double firstMeasTheta = dataset[0].theta;
    double potentialBeaconMinRange = firstMeasRange - 3*LIDAR_RANGE_STD;
    double potentialBeaconMaxRange = firstMeasRange + 3*LIDAR_RANGE_STD;
    std::vector<BeaconData> beacons = map.getBeaconsWithinRange(x, y, potentialBeaconMaxRange);
    std::vector<double> psiOptions;

    for (const BeaconData& beacon: beacons) {
        double beaconRange = sqrt((beacon.x - x)*(beacon.x - x) + (beacon.y - y)*(beacon.y - y));
        if (beaconRange > potentialBeaconMinRange) { // max range already captured by "getBeaconsWithinRange" above
            double beaconTheta = atan2((beacon.y - y), (beacon.x -x));
            psiOptions.push_back(wrapAngle(beaconTheta - firstMeasTheta));
            std::cout << "beaconX: " << beacon.x << ", beaconY: " << beacon.y << std::endl;
            std::cout << "beaconTheta: " << beaconTheta << ", measTheta: " << firstMeasTheta << ", psi: " << wrapAngle(beaconTheta - firstMeasTheta) << std::endl; 
        }
    }

    double bestPsi = 0;
    double minRMSE = std::numeric_limits<double>::max();
    for (const double& psi : psiOptions) {
        double RMSE = 0;
        VectorXd state = VectorXd(3);
        state << x, y, psi;
        for (const LidarMeasurement& meas : dataset) {
            int beaconID = associateBeaconID(state, meas, map);
            BeaconData beacon = map.getBeaconWithId(beaconID);
            RMSE += pow(getLidarDistanceError(state, meas, beacon),2);
        }
        if (RMSE < minRMSE) {
            bestPsi = psi;
            minRMSE = RMSE;
        }
        std::cout << "psi = " << psi << ", RMSE = " << RMSE << std::endl;
    }

    std::cout << "best psi: " << bestPsi << std::endl;
    return bestPsi;
}

bool withinChiSq1PctThreshold(double normalizedInnovationSquared, int numMeasStates) {
    double chiSquareThreshold;
    if (numMeasStates == 1) {
        chiSquareThreshold = 6.63;
    } else if (numMeasStates == 2) {
        chiSquareThreshold = 9.21;
    } else if (numMeasStates == 3) {
        chiSquareThreshold = 11.34;
    } else if (numMeasStates == 4) {
        chiSquareThreshold = 13.28;
    } else if (numMeasStates == 5) {
        chiSquareThreshold = 15.09;
    } else {
        std::cerr << "chiSquareTest1Pct only supports 1-5 measurement variables" << std::endl;
        return false;
    }
    if (normalizedInnovationSquared > chiSquareThreshold) {
        return false;
    } else {
        return true;
    }
}

bool withinChiSq5PctThreshold(double normalizedInnovationSquared, int numMeasStates) {
    double chiSquareThreshold;
    if (numMeasStates == 1) {
        chiSquareThreshold = 3.84;
    } else if (numMeasStates == 2) {
        chiSquareThreshold = 5.99;
    } else if (numMeasStates == 3) {
        chiSquareThreshold = 7.81;
    } else if (numMeasStates == 4) {
        chiSquareThreshold = 9.49;
    } else if (numMeasStates == 5) {
        chiSquareThreshold = 11.07;
    } else {
        std::cerr << "chiSquareTest5Pct only supports 1-5 measurement variables" << std::endl;
        return false;
    }
    if (normalizedInnovationSquared > chiSquareThreshold) {
        return false;
    } else {
        return true;
    }
}
// ----------------End helper functions ----------------------- //

void KalmanFilter::initializeKalmanFilter(const BeaconMap& map) {
    // attempt to initialize the filter using this set of lidar measurements and all GPS measurements up to this point

    // Method:
    // Ensure there are a minimum number of each measurement
    // Get the initial x and y estimates from the 1st gps measurement
    // Get the initial psi estimate by determining the best fit for the first lidar measurement set
    // Run the kalman filter for the remaining points across a sweep of possible speeds
        // determine the range of speeds by getting the range of gps-position-derivative speeds across the five samples with some margin
        // for each run, keep track of the normalised innovation error
        // select the starting velocity that resulted in the minimum innovation error at measurement 5

    std::vector<GPSMeasurement> gpsMeasurements = getGPSInitializationMeasurements();
    std::vector<std::vector<LidarMeasurement>> lidarDatasets = getLidarInitilizationDatasets();
    std::vector<GyroMeasurement> gyroMeasurements = getGyroInitializationMeasurements();
    int nLidar = lidarDatasets.size();
    int nGPS = gpsMeasurements.size();
    int nGyro = gyroMeasurements.size();
    if (nLidar < MIN_INIT_SAMPLES || nGPS < MIN_INIT_SAMPLES || nGyro < MIN_INIT_SAMPLES) {
        // not enough data to initilize yet
        return;
    }

    // x0 and y0 from 1st gps measurement
    double x0 = gpsMeasurements[0].x;
    double y0 = gpsMeasurements[0].y;

    // psi0 from 1st lidar measurement set
    double psi0 = estimatePsiFromLidar(x0, y0, lidarDatasets[0], map);
    
    // determine range of velocities to explore
    double minVelocity = 0;
    double maxVelocity = 0;
    double dt = gpsMeasurements[1].time - gpsMeasurements[0].time;
    for (int i = 0; i < gpsMeasurements.size()-1; i++) {
        double distance = sqrt(pow(gpsMeasurements[i+1].x-gpsMeasurements[i].x,2) + pow(gpsMeasurements[i+1].y-gpsMeasurements[i].y,2));
        double velocity = distance / dt;
        maxVelocity = (velocity > maxVelocity) ? velocity : maxVelocity;
        minVelocity = (velocity < minVelocity) ? velocity : minVelocity;
    }
    double velocityRange = maxVelocity - minVelocity;
    minVelocity = minVelocity - velocityRange/2.0;
    maxVelocity = maxVelocity + velocityRange/2.0;
    int nStep = 20; // explore 20 possible options for starting velocity within the possible range
    double vStep = (maxVelocity-minVelocity) / (nStep - 1); // this vStep will also be used as the initial velocity std

    // run the kalman filter through the initialization measurements for each starting velocity option
    double minInnovationSq = std::numeric_limits<double>::max();
    double bestv0;
    VectorXd bestState = Vector4d();
    MatrixXd bestCov = Matrix4d();
    for (int i = 0; i < nStep; i++) {
        double v0 = minVelocity + i*vStep;
        double thisInnovationSq = 0;
        
        VectorXd state = Vector4d();
        MatrixXd cov = Matrix4d::Zero();
        state << x0, y0, psi0, v0;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = cov(0,0);
        cov(2,2) = LIDAR_THETA_STD*LIDAR_THETA_STD;
        cov(3,3) = vStep*vStep;

        setState(state);
        setCovariance(cov);

        int iGyro = 1;
        int iGPS = 1;
        int iLidar = 1;
        int nGyro = gyroMeasurements.size();
        double nextGyroTime = gyroMeasurements[iGyro].time;
        double nextGPSTime = gpsMeasurements[iGPS].time;
        double nextLidarTime = lidarDatasets[iLidar][0].time;

        while (true) {
            if (iLidar == nLidar && iGPS == nGPS && iGyro == nGyro) {
                break;
            }
            if (nextLidarTime < nextGyroTime && nextLidarTime < nextGPSTime) {
                thisInnovationSq += handleLidarMeasurements(lidarDatasets[iLidar], map);
                iLidar++;
                nextLidarTime = (iLidar < nLidar) ? lidarDatasets[iLidar][0].time : std::numeric_limits<double>::max();
                // std::cout << "lidar update" << std::endl;
            }
            else if (nextGPSTime < nextGyroTime) {
                thisInnovationSq += handleGPSMeasurement(gpsMeasurements[iGPS]);
                iGPS++;
                nextGPSTime = (iGPS < nGPS) ? gpsMeasurements[iGPS].time : std::numeric_limits<double>::max();
                // std::cout << "GPS update" << std::endl;
            }
            else {
                double dt = gyroMeasurements[iGyro].time - gyroMeasurements[iGyro-1].time;
                predictionStep(gyroMeasurements[iGyro], dt);
                iGyro++;
                nextGyroTime = (iGyro < nGyro) ? gyroMeasurements[iGyro].time : std::numeric_limits<double>::max();
            }

        }

        if (thisInnovationSq < minInnovationSq) {
            std::cout << thisInnovationSq << std::endl;
            minInnovationSq = thisInnovationSq;
            bestv0 = v0;
            bestState = getState();
            bestCov = getCovariance();
        }
    }

    std::cout << "Best Initial State: x0_" << x0 << " y0_" << y0 << "psi0_" << psi0 << "v0_" << bestv0 << std::endl;

    // --- DEBUG SECTION --- //
    // re-run the filter starting at the best initial state and cout the steps
    // VectorXd state = Vector4d();
    // MatrixXd cov = Matrix4d::Zero();
    // state << x0, y0, psi0, bestv0;
    // cov(0,0) = GPS_POS_STD*GPS_POS_STD;
    // cov(1,1) = cov(0,0);
    // cov(2,2) = LIDAR_THETA_STD*LIDAR_THETA_STD;
    // cov(3,3) = vStep*vStep;

    // setState(state);
    // setCovariance(cov);

    // int iGyro = 1;
    // int iGPS = 1;
    // int iLidar = 1;
    // nGyro = gyroMeasurements.size();
    // double nextGyroTime = gyroMeasurements[iGyro].time;
    // double nextGPSTime = gpsMeasurements[iGPS].time;
    // double nextLidarTime = lidarDatasets[iLidar][0].time;
    // std::cout << "Initial state: " << std::endl << state << std::endl;

    // while (true) {
    //     if (iLidar == nLidar && iGPS == nGPS && iGyro == nGyro) {
    //         break;
    //     }
    //     if (nextLidarTime < nextGyroTime && nextLidarTime < nextGPSTime) {
    //         handleLidarMeasurements(lidarDatasets[iLidar], map);
    //         iLidar++;
    //         nextLidarTime = (iLidar < nLidar) ? lidarDatasets[iLidar][0].time : std::numeric_limits<double>::max();
    //         std::cout << "lidar update" << std::endl;
    //     }
    //     else if (nextGPSTime < nextGyroTime) {
    //         handleGPSMeasurement(gpsMeasurements[iGPS]);
    //         iGPS++;
    //         nextGPSTime = (iGPS < nGPS) ? gpsMeasurements[iGPS].time : std::numeric_limits<double>::max();
    //         std::cout << "GPS update" << std::endl;
    //     }
    //     else {
    //         double dt = gyroMeasurements[iGyro].time - gyroMeasurements[iGyro-1].time;
    //         predictionStep(gyroMeasurements[iGyro], dt);
    //         iGyro++;
    //         nextGyroTime = (iGyro < nGyro) ? gyroMeasurements[iGyro].time : std::numeric_limits<double>::max();
    //         std::cout << "Prediction step" << std::endl;
    //     }

    //     std::cout << "new state: " << std::endl << getState() << std::endl;

    // }
    

    // --- END DEBUG SECTION --- //

    setState(bestState);
    setCovariance(bestCov);
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt) {
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
    else 
    {
        addGyroInitializationMeasurement(gyro);
    }
}

double KalmanFilter::handleGPSMeasurement(GPSMeasurement meas) {
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
        VectorXd v = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*v;
        cov = (Matrix4d::Identity() - K*H) * cov;

        double nis = v.transpose()*S.inverse()*v;

        if (withinChiSq1PctThreshold(nis, v.size())) {
            // valid measurement
            setState(state);
            setCovariance(cov);
            return nis;
        } else {
            // invalid measurement; don't update state or covariance
            std::cout << "Invalid gps measurement discarded" << std::endl;
        }
        return nis;
    }
    else
    {
        if (isLidarEnabled()) {
            addGPSInitializationMeasurement(meas);
        }
        else
        {
            // No complex initialization for gps-only case; just set the position, assume zero velocity and psi, and start the filter
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
        return 0;
    } 
             
}

double KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map) {
    if (isInitialised()) {
        // Assume No Correlation between the Measurements and Update Sequentially
        double nisSum = 0;
        for(const auto& meas : dataset) {nisSum += handleLidarMeasurement(meas, map);}
        return nisSum;
    }
    else
    {
        std::cout << "not initialized" << std::endl;
        addLidarInitializationDatasets(dataset);
        initializeKalmanFilter(map);
        return 0;
    }
    
}

double KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map) {
    // returns normalized innovation squared
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();
        meas.id = associateBeaconID(state, meas,map);

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
        double nis = 1000; // arbitrarily high innovation squared (in case measurement doesn't have id)
        if (meas.id != -1 && map_beacon.id != -1)
        {           
            // The map matched beacon positions can be accessed using: map_beacon.x AND map_beacon.y
            Vector2d z = Vector2d(); // z = measurement (in eigen matrix format)
            Vector2d z_hat = Vector2d(); // z_hat = expected measurement based on current state estimate
            Vector2d v = Vector2d(); // v = measurement innovation
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
            v = z - z_hat;
            v(1) = wrapAngle(v(1));
            H << -dx_hat/z_hat_rng,-dy_hat/z_hat_rng,0,0, dy_hat/(z_hat_rng*z_hat_rng),-dx_hat/(z_hat_rng*z_hat_rng),-1,0;
            R << LIDAR_RANGE_STD*LIDAR_RANGE_STD,0, 0,LIDAR_THETA_STD*LIDAR_THETA_STD;
            S = H*cov*H.transpose() + R;
            K = cov*H.transpose()*S.inverse();

            state = state + K*v;
            state(2) = wrapAngle(state(2));
            cov = (Matrix4d::Identity() - K*H)*cov;

            double nis = v.transpose()*S.inverse()*v; // normalised innovation squared

            if (withinChiSq1PctThreshold(nis, v.size())) {
                // valid measurement
                setState(state);
                setCovariance(cov);
                std::cout << "GOOD LIDAR; NIS = " << nis << std::endl;
            } else {
                // invalid measurement; don't update state or covariance
                std::cout << "Invalid lidar measurement discarded" << std::endl;
                std::cout << "NIS: " << nis << std::endl;
                std::cout << "State: x_" << state[0] << ", y_" << state[1] << ", psi_" << state[2] << std::endl;
                std::cout << "Meas: range_" << meas.range << ", theta_" << meas.theta << std::endl;
                std::cout << "Expected meas: range_" << z_hat_rng << ", theta_" << z_hat_ang << std::endl;
                std::cout << "Beacon: x_" << map_beacon.x << ", y_" << map_beacon.y << std::endl;
                std::cout << "Expected beacon: x_" << state[0]+dx_hat << ", y_" << state[1]+dy_hat << std::endl;
            }
        }

        // ----------------------------------------------------------------------- //
        return nis;
    }
    return 0;
}

VectorXd KalmanFilter::lidarMeasurementModel(VectorXd state, double beaconX, double beaconY){
    VectorXd z_hat = VectorXd::Zero(2);

    double x = state(0);
    double y = state(1);
    double psi = state(2);

    double dx = beaconX - x;
    double dy = beaconY - y;
    double r = sqrt(dx*dx + dy*dy);
    double theta = wrapAngle(atan2(dy,dx) - psi);

    z_hat << r,theta;

    return z_hat;
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance() {
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState() {
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
