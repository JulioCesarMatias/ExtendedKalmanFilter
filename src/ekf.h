#pragma once

#include "ekfMath.h"

// acceleration due to gravity in m/s/s
#define GRAVITY_MSS 9.80665f

extern float rotX, rotY, rotZ;
extern float gForceX, gForceY, gForceZ;
extern int32_t roll_sensor;
extern int32_t pitch_sensor;
extern int32_t yaw_sensor;
void ekf_update(void);

typedef union
{
    struct
    {
        uint16_t attitude : 1;           // 0 - true if attitude estimate is valid
        uint16_t horiz_vel : 1;          // 1 - true if horizontal velocity estimate is valid
        uint16_t vert_vel : 1;           // 2 - true if the vertical velocity estimate is valid
        uint16_t horiz_pos_rel : 1;      // 3 - true if the relative horizontal position estimate is valid
        uint16_t horiz_pos_abs : 1;      // 4 - true if the absolute horizontal position estimate is valid
        uint16_t vert_pos : 1;           // 5 - true if the vertical position estimate is valid
        uint16_t terrain_alt : 1;        // 6 - true if the terrain height estimate is valid
        uint16_t const_pos_mode : 1;     // 7 - true if we are in const position mode
        uint16_t pred_horiz_pos_rel : 1; // 8 - true if filter expects it can produce a good relative horizontal position estimate - used before takeoff
        uint16_t pred_horiz_pos_abs : 1; // 9 - true if filter expects it can produce a good absolute horizontal position estimate - used before takeoff
        uint16_t takeoff_detected : 1;   // 10 - true if optical flow takeoff has been detected
        uint16_t takeoff : 1;            // 11 - true if filter is compensating for baro errors during takeoff
        uint16_t touchdown : 1;          // 12 - true if filter is compensating for baro errors during touchdown
        uint16_t using_gps : 1;          // 13 - true if we are using GPS position
    } flags;
    uint16_t value;
} nav_filter_status;

typedef struct
{
    fpQuaternion_t quat;        // 0..3
    fpVector3_t velocity;       // 4..6
    fpVector3_t position;       // 7..9
    fpVector3_t gyro_bias;      // 10..12
    float accel_zbias1;         // 13
    fpVector3_t wind_vel;       // 14..15
    fpVector3_t earth_magfield; // 16..18
    fpVector3_t body_magfield;  // 19..21
    float accel_zbias2;         // 22
    fpVector3_t vel1;           // 23 .. 25
    float posD1;                // 26
    fpVector3_t vel2;           // 27 .. 29
    float posD2;                // 30
    fpVector3_t omega;          // 31 .. 33
} state_elements;

// states held by optical flow fusion across time steps
// optical flow X,Y motion compensated rate measurements are fused across two time steps
// to level computational load as this can be an expensive operation
typedef struct
{
    uint8_t obsIndex;
    float SH_LOS[4];
    float SK_LOS[10];
    float q0;
    float q1;
    float q2;
    float q3;
    float vn;
    float ve;
    float vd;
    float pd;
    float losPred[2];
} flow_state_s;

typedef struct
{
    bool bad_xmag : 1;
    bool bad_ymag : 1;
    bool bad_zmag : 1;
    bool bad_airspeed : 1;
    bool bad_sideslip : 1;
} faultStatus_s;

// states held by magnetomter fusion across time steps
// magnetometer X,Y,Z measurements are fused across three time steps
// to level computational load as this is an expensive operation
typedef struct
{
    float q0;
    float q1;
    float q2;
    float q3;
    float magN;
    float magE;
    float magD;
    float magXbias;
    float magYbias;
    float magZbias;
    uint8_t obsIndex;
    float DCM[3][3];
    float MagPred[3];
    float R_MAG;
    float SH_MAG[9];
} mag_state_s;

// This function is used to initialise the filter whilst moving, using the AHRS DCM solution
// It should NOT be used to re-initialise after a timeout as DCM will also be corrupted
bool ekf_InitialiseFilterDynamic(void);

// Initialise the states from accelerometer and magnetometer data (if present)
// This method can only be used when the vehicle is static
bool ekf_InitialiseFilterBootstrap(void);

// Update Filter States - this should be called whenever new IMU data is available
void ekf_UpdateFilter(void);

// Check basic filter health metrics and return a consolidated health status
bool ekf_healthy(void);

// Return the last calculated NED position relative to the reference point (m).
// If a calculated solution is not available, use the best available data and return false
// If false returned, do not use for flight control
bool ekf_getPosNED(fpVector3_t *pos);

// return NED velocity in m/s
void ekf_getVelNED(fpVector3_t *vel);

// This returns the specific forces in the NED frame
void ekf_getAccelNED(fpVector3_t *accelNED);

// return body axis gyro bias estimates in rad/sec
void ekf_getGyroBias(fpVector3_t *gyroBias);

// reset body axis gyro bias estimates
void ekf_resetGyroBias(void);

// Resets the baro so that it reads zero at the current height
// Resets the EKF height to zero
// Adjusts the EKf origin height so that the EKF height + origin height is the same as before
// Returns true if the height datum reset has been performed
// If using a range finder for height no reset is performed and it returns false
bool ekf_resetHeightDatum(void);

// Commands the EKF to not use GPS.
// This command must be sent prior to arming as it will only be actioned when the filter is in static mode
// This command is forgotten by the EKF each time it goes back into static mode (eg the vehicle disarms)
// Returns 0 if command rejected
// Returns 1 if attitude, vertical velocity and vertical position will be provided
// Returns 2 if attitude, 3D-velocity, vertical position and relative horizontal position will be provided
uint8_t ekf_setInhibitGPS(void);

// return the horizontal speed limit in m/s set by optical flow sensor limits
// return the scale factor to be applied to navigation velocity gains to compensate for increase in velocity noise with height when using optical flow
void ekf_getEkfControlLimits(float *ekfGndSpdLimit, float *ekfNavVelGainScaler);

// return weighting of first IMU in blending function
void ekf_getIMU1Weighting(float *ret);

// return the individual Z-accel bias estimates in m/s^2
void ekf_getAccelZBias(float *zbias1, float *zbias2);

// return the NED wind speed estimates in m/s (positive is air moving in the direction of the axis)
void ekf_getWind(fpVector3_t *wind);

// return earth magnetic field estimates in measurement units / 1000
void ekf_getMagNED(fpVector3_t *magNED);

// return body magnetic field estimates in measurement units / 1000
void ekf_getMagXYZ(fpVector3_t *magXYZ);

// Return estimated magnetometer offsets
// Return true if magnetometer offsets are valid
bool ekf_getMagOffsets(fpVector3_t *magOffsets);

// return the latitude and longitude and height used to set the NED origin
// All NED positions calculated by the filter are relative to this location
// Returns false if the origin has not been set
// bool getOriginLLH(struct Location *loc);

// set the latitude and longitude and height used to set the NED origin
// All NED positions calcualted by the filter will be relative to this location
// The origin cannot be set if the filter is in a flight mode (eg vehicle armed)
// Returns false if the filter has rejected the attempt to set the origin
// bool setOriginLLH(struct Location *loc);

// return estimated height above ground level
// return false if ground height is not being estimated.
bool ekf_getHAGL(float *HAGL);

// return the Euler roll, pitch and yaw angle in radians
void ekf_getEulerAngles(fpVector3_t *eulers);

// return the transformation matrix from XYZ (body) to NED axes
void ekf_getRotationBodyToNED(fpMatrix3_t *mat);

// return the quaternions defining the rotation from NED to XYZ (body) axes
void ekf_getQuaternion(fpQuaternion_t *quat);

// return the innovations for the NED Pos, NED Vel, XYZ Mag and Vtas measurements
void ekf_getInnovations(fpVector3_t *velInnov, fpVector3_t *posInnov, fpVector3_t *magInnov, float *tasInnov);

// return the innovation consistency test ratios for the velocity, position, magnetometer and true airspeed measurements
void ekf_getVariances(float *velVar, float *posVar, float *hgtVar, fpVector3_t *magVar, float *tasVar, fpVector3_t *offset);

// should we use the compass? This is public so it can be used for
// reporting via ahrs.ekf_use_compass()
bool ekf_use_compass(void);

// write the raw optical flow measurements
// rawFlowQuality is a measured of quality between 0 and 255, with 255 being the best quality
// rawFlowRates are the optical flow rates in rad/sec about the X and Y sensor axes.
// rawGyroRates are the sensor rotation rates in rad/sec measured by the sensors internal gyro
// The sign convention is that a RH physical rotation of the sensor about an axis produces both a positive flow and gyro rate
// msecFlowMeas is the scheduler time in msec when the optical flow data was received from the sensor.
void ekf_writeOptFlowMeas(uint8_t *rawFlowQuality, fpVector3_t *rawFlowRates, fpVector3_t *rawGyroRates, uint32_t *msecFlowMeas);

// called by vehicle code to specify that a takeoff is happening
// causes the EKF to compensate for expected barometer errors due to ground effect
void ekf_setTakeoffExpected(bool val);

// called by vehicle code to specify that a touchdown is expected to happen
// causes the EKF to compensate for expected barometer errors due to ground effect
void ekf_setTouchdownExpected(bool val);

/*
return the filter fault status as a bitmasked integer
 0 = quaternions are NaN
 1 = velocities are NaN
 2 = badly conditioned X magnetometer fusion
 3 = badly conditioned Y magnetometer fusion
 5 = badly conditioned Z magnetometer fusion
 6 = badly conditioned airspeed fusion
 7 = badly conditioned synthetic sideslip fusion
 7 = filter is not initialised
*/
void ekf_getFilterFaults(uint8_t *faults);

/*
return filter timeout status as a bitmasked integer
 0 = position measurement timeout
 1 = velocity measurement timeout
 2 = height measurement timeout
 3 = magnetometer measurement timeout
 5 = unassigned
 6 = unassigned
 7 = unassigned
 7 = unassigned
*/
void ekf_getFilterTimeouts(uint8_t *timeouts);

/*
return filter status flags
*/
void ekf_getFilterStatus(nav_filter_status *status);

// provides the height limit to be observed by the control loops
// returns false if no height limiting is required
// this is needed to ensure the vehicle does not fly too high when using optical flow navigation
bool ekf_getHeightControlLimit(float *height);

// provides the quaternion that was used by the INS calculation to rotate from the previous orientation to the orientaion at the current time step
// returns a zero rotation quaternion if the INS calculation was not performed on that time step.
fpQuaternion_t ekf_getDeltaQuaternion(void);

// return the amount of yaw angle change due to the last yaw angle reset in radians
// returns true if a reset yaw angle has been updated and not queried
// this function should not have more than one client
bool ekf_getLastYawResetAngle(float *yawAng);

// update the quaternion, velocity and position states using IMU measurements
void ekf_UpdateStrapdownEquationsNED(void);

// calculate the predicted state covariance matrix
void ekf_CovariancePrediction(void);

// force symmetry on the state covariance matrix
void ekf_ForceSymmetry(void);

// copy covariances across from covariance prediction calculation and fix numerical errors
void ekf_CopyAndFixCovariances(void);

// constrain variances (diagonal terms) in the state covariance matrix
void ekf_ConstrainVariances(void);

// constrain states
void ekf_ConstrainStates(void);

// fuse selected position, velocity and height measurements
void ekf_FuseVelPosNED(void);

// fuse magnetometer measurements
void ekf_FuseMagnetometer(void);

// fuse true airspeed measurements
void ekf_Fekf_useAirspeed(void);

// fuse sythetic sideslip measurement of zero
void ekf_FuseSideslip(void);

// zero specified range of rows in the state covariance matrix
void ekf_zeroRows(float covMat[22][22], uint8_t first, uint8_t last);

// zero specified range of columns in the state covariance matrix
void ekf_zeroCols(float covMat[22][22], uint8_t first, uint8_t last);

// store states along with system time stamp in msces
void ekf_StoreStates(void);

// Reset the stored state history and store the current state
void ekf_StoreStatesReset(void);

// recall state vector stored at closest time to the one specified by msec
void ekf_RecallStates(state_elements *statesForFusion, uint32_t msec);

// calculate the NED earth spin vector in rad/sec
void ekf_calcEarthRateNED(fpVector3_t *omega, int32_t latitude);

// calculate whether the flight vehicle is on the ground or flying from height, airspeed and GPS speed
void ekf_SetFlightAndFusionModes(void);

// initialise the covariance matrix
void ekf_CovarianceInit(void);

// update IMU delta angle and delta velocity measurements
void ekf_readIMUData(void);

// check for new valid GPS data and update stored measurement if available
void ekf_readGpsData(void);

// check for new altitude measurement data and update stored measurement if available
void ekf_readHgtData(void);

// check for new magnetometer data and update store measurements if available
void ekf_readMagData(void);

// check for new airspeed data and update stored measurements if available
void ekf_readAirSpdData(void);

// determine when to perform fusion of GPS position and  velocity measurements
void ekf_SelectVelPosFusion(void);

// determine when to perform fusion of true airspeed measurements
void ekf_SelectTasFusion(void);

// determine when to perform fusion of synthetic sideslp measurements
void ekf_SelectBetaFusion(void);

// determine when to perform fusion of magnetometer measurements
void ekf_SelectMagFusion(void);

// force alignment of the yaw angle using GPS velocity data
void ekf_alignYawGPS(void);

// Forced alignment of the wind velocity states so that they are set to the reciprocal of
// the ground speed and scaled to 6 m/s. This is used when launching a fly-forward vehicle without an airspeed sensor
// on the assumption that launch will be into wind and 6 is representative global average at height
// http://maps.google.com/gallery/details?id=zJuaSgXp_WLc.kTBytKPmNODY*hl=en
void ekf_setWindVelStates();

// initialise the earth magnetic field states using declination and current attitude and magnetometer meaasurements
// and return attitude quaternion
fpQuaternion_t ekf_calcQuatAndFieldStates(float roll, float pitch);

// zero stored variables
void ekf_InitialiseVariables(void);

// reset the horizontal position states uing the last GPS measurement
void ekf_ResetPosition(void);

// reset velocity states using the last GPS measurement
void ekf_ResetVelocity(void);

// reset the vertical position state using the last height measurement
void ekf_ResetHeight(void);

// return true if we should use the airspeed sensor
bool ekf_useAirspeed(void);

// decay GPS horizontal position offset to close to zero at a rate of 1 m/s
// this allows large GPS position jumps to be accomodated gradually
void ekf_decayGpsOffset(void);

// return true if optical flow data is available
bool ekf_optFlowDataPresent(void);

// determine when to perform fusion of optical flow measurements
void ekf_SelectFlowFusion();

// recall omega (angular rate vector) average from time specified by msec to current time
// this is useful for motion compensation of optical flow measurements
void ekf_RecallOmega(fpVector3_t *omegaAvg, uint32_t msecStart, uint32_t msecEnd);

// Estimate terrain offset using a single state EKF
void ekf_EstimateTerrainOffset(void);

// fuse optical flow measurements into the main filter
void ekf_FuseOptFlow(void);

// Check arm status and perform required checks and mode changes
void ekf_performArmingChecks(void);

// Set the NED origin to be used until the next filter reset
void ekf_setOrigin(void);

// determine if a takeoff is expected so that we can compensate for expected barometer errors due to ground effect
bool ekf_getTakeoffExpected(void);

// determine if a touchdown is expected so that we can compensate for expected barometer errors due to ground effect
bool ekf_getTouchdownExpected();

// Assess GPS data quality and return true if good enough to align the EKF
bool ekf_calcGpsGoodToAlign(void);

// Read the range finder and take new measurements if available
// Apply a median filter to range finder data
void ekf_readRangeFinder(void);

// check if the vehicle has taken off during optical flow navigation by looking at inertial and range finder data
void ekf_detectOptFlowTakeoff(void);

// align the NE earth magnetic field states with the published declination
void ekf_alignMagStateDeclination(void);
