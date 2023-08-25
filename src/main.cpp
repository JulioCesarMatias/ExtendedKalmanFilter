#include "ekf.h"
#include <Wire.h>

TwoWire HWire(2, I2C_FAST_MODE); // Initiate I2C port 2 at 400kHz

enum gyro_fsr_e
{
  INV_FSR_250DPS = 0,
  INV_FSR_500DPS = 1,
  INV_FSR_1000DPS = 2,
  INV_FSR_2000DPS = 3,
  NUM_GYRO_FSR
};

enum accel_fsr_e
{
  INV_FSR_2G = 0,
  INV_FSR_4G = 1,
  INV_FSR_8G = 2,
  INV_FSR_16G = 4,
  NUM_ACCEL_FSR
};

enum lpf_e
{
  INV_FILTER_256HZ_NOLPF2 = 0,
  INV_FILTER_188HZ = 1,
  INV_FILTER_98HZ = 2,
  INV_FILTER_42HZ = 3,
  INV_FILTER_20HZ = 4,
  INV_FILTER_10HZ = 5,
  INV_FILTER_5HZ = 6,
  INV_FILTER_2100HZ_NOLPF = 7,
  NUM_FILTER
};

enum task_e
{
  IMU_TASK = 0,
  EKF_TASK,

  IMU_DELTA_TIME,
  EKF_DELTA_TIME,

  SIZE_OF_TASK_ENUM
};

// Used to ignore results for functions marked as warn unused
#define IGNORE_RETURN(x) \
  do                     \
  {                      \
    if (x)               \
    {                    \
    }                    \
  } while (0)

// Gyroscope scale (uncertain where the 0.01745 value comes from)
#define MPU6050_GYRO_SCALE_2000 (0.0174532f / 16.4f)
#define MPU6050_GYRO_SCALE_1000 (0.0174532f / 32.8f)
#define MPU6050_GYRO_SCALE_500 (0.0174532f / 65.5f)
#define MPU6050_GYRO_SCALE_250 (0.0174532f / 131f)

// Accelerometer scale adjustment
#define MPU6050_ACCEL_SCALE_16G (GRAVITY_MSS / 2048.0f)
#define MPU6050_ACCEL_SCALE_8G (GRAVITY_MSS / 4096.0f)
#define MPU6050_ACCEL_SCALE_4G (GRAVITY_MSS / 8192.0f)
#define MPU6050_ACCEL_SCALE_2G (GRAVITY_MSS / 16384.0f)

#define MPU6050_ADDR 0x68
#define IMU_RATE_HZ 1000 // Value in Hz. Range:4Hz to 1Khz

uint32_t task_started_usec[SIZE_OF_TASK_ENUM];
uint32_t task_previous_usec[SIZE_OF_TASK_ENUM];
uint32_t task_dt[SIZE_OF_TASK_ENUM];

#define RUN_TASK(name, func, hz)                                          \
  task_started_usec[name] = micros();                                     \
  task_dt[name] = task_started_usec[name] - task_previous_usec[name];     \
  if (task_dt[name] >= (1000000 / hz))                                    \
  {                                                                       \
    if (task_dt[name] <= (1000000 / hz) + (((1000000 / hz)) * 0.20f))     \
    {                                                                     \
      func;                                                               \
      task_previous_usec[name + (SIZE_OF_TASK_ENUM / 2)] = task_dt[name]; \
    }                                                                     \
    task_previous_usec[name] = task_started_usec[name];                   \
  }

#define GET_TASK_DELTA_TIME(name, getDeltaTime) getDeltaTime = (float)task_previous_usec[name] * 1.0e-6f;

// #define PRINTLN_IMU_DATA

bool IMU_Healthy;
bool gyroCalibrationOk;
fpVector3_t gyroDataOffSet;

void readIMURegisters(void)
{
  if (!IMU_Healthy)
  {
    return;
  }

  HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
  HWire.write(0x3B);                     // Start reading
  HWire.endTransmission();               // End the transmission
  HWire.requestFrom(MPU6050_ADDR, 14);   // Request 14 bytes from the MPU-6050

  accData.x = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_ACCEL_SCALE_2G;
  accData.y = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_ACCEL_SCALE_2G;
  accData.z = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_ACCEL_SCALE_2G;

  IGNORE_RETURN((HWire.read() << 8) | HWire.read()); // Temperature

  gyroData.x = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_GYRO_SCALE_2000;
  gyroData.y = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_GYRO_SCALE_2000;
  gyroData.z = ((int16_t)((HWire.read() << 8) | HWire.read())) * MPU6050_GYRO_SCALE_2000;

  if (gyroCalibrationOk)
  {
    gyroData.x -= gyroDataOffSet.x;
    gyroData.y -= gyroDataOffSet.y;
    gyroData.z -= gyroDataOffSet.z;
  }
}

void initGyroCalibtation(void)
{
  fpVector3_t last_average;
  fpVector3_t best_avg;
  fpVector3_t new_gyro_offset;
  float best_diff;
  bool converged;

  Serial.print("Init Gyro Calibration");

  // Remove existing gyro offsets
  gyroDataOffSet.x = 0.0f;
  gyroDataOffSet.y = 0.0f;
  gyroDataOffSet.z = 0.0f;
  best_avg.x = 0.0f;
  best_avg.y = 0.0f;
  best_avg.z = 0.0f;
  new_gyro_offset.x = 0.0f;
  new_gyro_offset.y = 0.0f;
  new_gyro_offset.z = 0.0f;
  best_diff = 0;
  last_average.x = 0.0f;
  last_average.y = 0.0f;
  last_average.z = 0.0f;
  converged = false;

  for (uint8_t i = 0; i < 5; i++)
  {
    delay(5);
    readIMURegisters();
  }

  // The strategy is to average 50 points over 0.5 seconds, then do it again and see if the 2nd average is within a small margin of the first

  uint8_t num_converged = 0;

  // We try to get a good calibration estimate for up to 30 seconds if the gyro are stable, we should get it in 1 second
  for (int16_t j = 0; j <= 120 && num_converged < 1; j++)
  {
    fpVector3_t gyro_sum;
    fpVector3_t gyro_avg;
    fpVector3_t gyro_diff;
    fpVector3_t accel_start;
    float diff_norm;
    uint8_t i;

    Serial.print(".");

    gyro_sum.x = 0.0f;
    gyro_sum.y = 0.0f;
    gyro_sum.z = 0.0f;

    accel_start.x = accData.x;
    accel_start.y = accData.y;
    accel_start.z = accData.z;

    for (i = 0; i < 50; i++)
    {
      readIMURegisters();
      gyro_sum.x += gyroData.x;
      gyro_sum.y += gyroData.y;
      gyro_sum.z += gyroData.z;
      delay(5);
    }

    fpVector3_t accel_diff = {.v = {accData.x - accel_start.x, accData.y - accel_start.y, accData.z - accel_start.z}};

    if (calc_length_pythagorean_3D(accel_diff.x, accel_diff.y, accel_diff.z) > 0.2f)
    {
      // The accelerometers changed during the gyro sum. Skip this sample.
      // This copes with doing gyro cal on a  steadily moving platform. The value 0.2 corresponds with around 5 degrees/second of rotation.
      // continue;
    }

    gyro_avg.x = gyro_sum.x / i;
    gyro_avg.y = gyro_sum.y / i;
    gyro_avg.z = gyro_sum.z / i;
    gyro_diff.x = last_average.x - gyro_avg.x;
    gyro_diff.y = last_average.y - gyro_avg.y;
    gyro_diff.z = last_average.z - gyro_avg.z;
    diff_norm = calc_length_pythagorean_3D(gyro_diff.x, gyro_diff.y, gyro_diff.z);

    if (j == 0)
    {
      best_diff = diff_norm;
      best_avg.x = gyro_avg.x;
      best_avg.y = gyro_avg.y;
      best_avg.z = gyro_avg.z;
    }
    else if (calc_length_pythagorean_3D(gyro_diff.x, gyro_diff.y, gyro_diff.z) < (0.1f * (PI / 180.0f)))
    {
      // We want the average to be within 0.1 bit, which is 0.04 degrees/s
      last_average.x = (gyro_avg.x * 0.5f) + (last_average.x * 0.5f);
      last_average.y = (gyro_avg.y * 0.5f) + (last_average.y * 0.5f);
      last_average.z = (gyro_avg.z * 0.5f) + (last_average.z * 0.5f);
      if (!converged || calc_length_pythagorean_3D(last_average.x, last_average.y, last_average.z) < calc_length_pythagorean_3D(new_gyro_offset.x, new_gyro_offset.y, new_gyro_offset.z))
      {
        new_gyro_offset.x = last_average.x;
        new_gyro_offset.y = last_average.y;
        new_gyro_offset.z = last_average.z;
      }
      if (!converged)
      {
        converged = true;
        num_converged++;
      }
    }
    else if (diff_norm < best_diff)
    {
      best_diff = diff_norm;
      best_avg.x = (gyro_avg.x * 0.5f) + (last_average.x * 0.5f);
      best_avg.y = (gyro_avg.y * 0.5f) + (last_average.y * 0.5f);
      best_avg.z = (gyro_avg.z * 0.5f) + (last_average.z * 0.5f);
    }
    last_average.x = gyro_avg.x;
    last_average.y = gyro_avg.y;
    last_average.z = gyro_avg.z;
  }

  if (!converged)
  {
    Serial.println("");
    Serial.print("Gyro calibration fail diff: ");
    Serial.print(best_diff * (180.0f / PI));
    Serial.println(" dps");
    gyroDataOffSet.x = best_avg.x;
    gyroDataOffSet.y = best_avg.y;
    gyroDataOffSet.z = best_avg.z;
    gyroCalibrationOk = false;
  }
  else
  {
    Serial.println("");
    Serial.print("Gyro calibration X: ");
    Serial.print(new_gyro_offset.x * (180.0f / PI));
    Serial.print(" Gyro calibration Y: ");
    Serial.print(new_gyro_offset.y * (180.0f / PI));
    Serial.print(" Gyro calibration Z: ");
    Serial.println(new_gyro_offset.z * (180.0f / PI));
    gyroCalibrationOk = true;
    gyroDataOffSet.x = new_gyro_offset.x;
    gyroDataOffSet.y = new_gyro_offset.y;
    gyroDataOffSet.z = new_gyro_offset.z;
  }
}

void setup()
{
  delay(8000);
  Serial.begin(115200);

  HWire.begin();
  delay(150);

  HWire.beginTransmission(MPU6050_ADDR);    // Start communication with the MPU-6050
  uint8_t status = HWire.endTransmission(); // End the transmission and register the exit status

  if (status == 0)
  {
    Serial.println("MPU6050 Healthy");
    IMU_Healthy = true;
  }
  else
  {
    Serial.println("MPU6050 bad");
  }

  if (IMU_Healthy)
  {
    HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
    HWire.write(0x6B);                     // We want to write to the PWR_MGMT_1 register
    HWire.write(0x00);                     // Set the register bits as 0x00 to activate the gyro
    HWire.endTransmission();               // End the transmission with the gyro

    HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
    HWire.write(0x1B);                     // We want to write to the GYRO_CONFIG register
    HWire.write(INV_FSR_2000DPS << 3);     // 2000dps full scale
    HWire.endTransmission();               // End the transmission with the gyro

    HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
    HWire.write(0x1C);                     // We want to write to the ACCEL_CONFIG register
    HWire.write(INV_FSR_2G << 3);          // +/- 2g full scale range
    HWire.endTransmission();               // End the transmission with the gyro

    HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
    HWire.write(0x19);                     // We want to write to the RATE_DIV register
    HWire.write((1000 / IMU_RATE_HZ) - 1); // Set sampling rate in 1KHz
    HWire.endTransmission();               // End the transmission with the gyro

    HWire.beginTransmission(MPU6050_ADDR); // Start communication with the MPU-6050
    HWire.write(0x1A);                     // We want to write to the CONFIG register
    HWire.write(INV_FILTER_20HZ);          // Set Digital Low Pass Filter
    HWire.endTransmission();               // End the transmission with the gyro

    initGyroCalibtation();
  }
}

void loop()
{
  RUN_TASK(IMU_TASK, readIMURegisters(), 1000); // IMU Readings running at 1KHz
  GET_TASK_DELTA_TIME(IMU_DELTA_TIME, accUpdate);

  RUN_TASK(EKF_TASK, ekf_update(), 400); // EKF calculations running at 400Hz

#ifdef PRINTLN_IMU_DATA
  Serial.print(accData.x);
  Serial.print(" ");
  Serial.print(accData.y);
  Serial.print(" ");
  Serial.print(accData.z);

  Serial.print("      ");

  Serial.print(gyroData.x);
  Serial.print(" ");
  Serial.print(gyroData.y);
  Serial.print(" ");
  Serial.println(gyroData.z);
#endif
}