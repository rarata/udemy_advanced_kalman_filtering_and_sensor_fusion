#ifndef INCLUDE_AKFSFSIM_KALMANFILTER_H
#define INCLUDE_AKFSFSIM_KALMANFILTER_H

#include <vector>
#include <Eigen/Dense>

#include "car.h"
#include "sensors.h"
#include "beacons.h"

using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector4d;

using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix4d;

class KalmanFilterBase
{
    public:

        KalmanFilterBase():m_initialised(false){}
        virtual ~KalmanFilterBase(){}
        void reset(bool lidar_enabled){m_initialised = false; m_gps_initialization_measurements.clear(); m_lidar_initialization_datasets.clear(); m_gyro_initialization_measurements.clear(); m_lidar_enabled = lidar_enabled;}
        bool isInitialised() const {return m_initialised;}
        bool isLidarEnabled() const {return m_lidar_enabled;}

    protected:
        VectorXd getState() const {return m_state;}
        MatrixXd getCovariance()const {return m_covariance;}
        void setState(const VectorXd& state ) {m_state = state; m_initialised = true;}
        void setCovariance(const MatrixXd& cov ){m_covariance = cov;}
        void addGPSInitializationMeasurement(const GPSMeasurement meas) {m_gps_initialization_measurements.push_back(meas);}
        void addLidarInitializationDatasets(const std::vector<LidarMeasurement> dataset) {m_lidar_initialization_datasets.push_back(dataset);}
        void addGyroInitializationMeasurement(const GyroMeasurement meas) {m_gyro_initialization_measurements.push_back(meas);}
        std::vector<GPSMeasurement> getGPSInitializationMeasurements() const {return m_gps_initialization_measurements;}
        std::vector<std::vector<LidarMeasurement>> getLidarInitilizationDatasets() const {return m_lidar_initialization_datasets;}
        std::vector<GyroMeasurement> getGyroInitializationMeasurements() const {return m_gyro_initialization_measurements;}

    private:
        bool m_lidar_enabled;
        bool m_initialised;
        VectorXd m_state;
        MatrixXd m_covariance;
        std::vector<GPSMeasurement> m_gps_initialization_measurements;
        std::vector<std::vector<LidarMeasurement>> m_lidar_initialization_datasets;
        std::vector<GyroMeasurement> m_gyro_initialization_measurements;
};

class KalmanFilter : public KalmanFilterBase
{
    public:

        VehicleState getVehicleState();
        Matrix2d getVehicleStatePositionCovariance();

        void attemptInitialization(const BeaconMap& map);
        void predictionStep(double dt);
        void predictionStep(GyroMeasurement gyro, double dt);
        double handleLidarMeasurements(const std::vector<LidarMeasurement>& meas, const BeaconMap& map);
        double handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map);
        double handleGPSMeasurement(GPSMeasurement meas);

};

#endif  // INCLUDE_AKFSFSIM_KALMANFILTER_H