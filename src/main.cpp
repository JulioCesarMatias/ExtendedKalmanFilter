#include "ekf.h"
#include <Wire.h>

int16_t accelX, accelY, accelZ;
int16_t gyroX, gyroY, gyroZ;

void setup()
{
  delay(8000);
  Serial.begin(115200);

  Wire.begin();

  Wire.beginTransmission(0x68); // This is the I2C address of the MPU (b1101000/b1101001 for AC0 low/high datasheet sec. 9.2)
  Wire.write(0x6B);             // Accessing the register 6B - Power Management (Sec. 4.28)
  Wire.write(0x03);             // Setting SLEEP register to 0. (Required; see Note on p. 9)
  Wire.endTransmission();

  Wire.beginTransmission(0x68); // I2C address of the MPU
  Wire.write(0x1B);             // Accessing the register 1B - Gyroscope Configuration (Sec. 4.4)
  Wire.write(0x18);             // Setting the gyro to full scale (2000 deg/sec)
  Wire.endTransmission();

  Wire.beginTransmission(0x68); // I2C address of the MPU
  Wire.write(0x1C);             // Accessing the register 1C - Acccelerometer Configuration (Sec. 4.5)
  Wire.write(0x10);             // Setting the accel to +/- 8g
  Wire.endTransmission();
}

void readAccelRegisters(void)
{
  Wire.beginTransmission(0x68); // I2C address of the MPU
  Wire.write(0x3B);             // Starting register for Accel Readings
  Wire.endTransmission();

  Wire.requestFrom(0x68, 6); // Request Accel Registers (3B - 40)

  while (Wire.available() < 6)
  {
  }

  accelX = (Wire.read() << 8) | Wire.read(); // Store first two bytes into accelX
  accelY = (Wire.read() << 8) | Wire.read(); // Store middle two bytes into accelY
  accelZ = (Wire.read() << 8) | Wire.read(); // Store last two bytes into accelZ

  gForceX = accelX / 16384.0f;
  gForceY = accelY / 16384.0f;
  gForceZ = accelZ / 16384.0f;
}

void readGyroRegisters(void)
{
  Wire.beginTransmission(0x68); // I2C address of the MPU
  Wire.write(0x43);             // Starting register for Gyro Readings
  Wire.endTransmission();

  Wire.requestFrom(0x68, 6); // Request Gyro Registers (43 - 48)

  while (Wire.available() < 6)
  {
  }

  gyroX = (Wire.read() << 8) | Wire.read(); // Store first two bytes into accelX
  gyroY = (Wire.read() << 8) | Wire.read(); // Store middle two bytes into accelY
  gyroZ = (Wire.read() << 8) | Wire.read(); // Store last two bytes into accelZ

  rotX = gyroX / 131.0f;
  rotY = gyroY / 131.0f;
  rotZ = gyroZ / 131.0f;
}

void loop()
{
  readAccelRegisters();
  readGyroRegisters();

  ekf_update();

  delay(1); // 1KHz loop
}