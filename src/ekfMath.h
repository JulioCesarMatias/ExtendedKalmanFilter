#pragma once

#include "Arduino.h"

#include <stdbool.h>
#include <stdint.h>

typedef struct {
    float q0, q1, q2, q3;
} fpQuaternion_t;

typedef union {
    float v[3];
    struct {
       float x,y,z;
    };
} fpVector3_t;

typedef struct {
    float m[3][3];
} fpMatrix3_t;

void quaternion_initialise(fpQuaternion_t *q);
void quaternion_from_euler(fpQuaternion_t *q, float roll, float pitch, float yaw);
void quaternion_from_axis_angle(fpVector3_t v, fpQuaternion_t *q);
void quaternion_multiply(fpQuaternion_t qv, fpQuaternion_t *q);
void quaternion_normalize(fpQuaternion_t *q);
void quaternion_rotation_matrix(fpQuaternion_t q, fpMatrix3_t *m);
void quaternion_to_euler(fpQuaternion_t q, float *roll, float *pitch, float *yaw);
void matrix_from_euler(fpMatrix3_t *m, float roll, float pitch, float yaw);
fpVector3_t multiply_matrix_by_vector(fpMatrix3_t m, fpVector3_t v);
void matrix_transpose3x3(fpMatrix3_t mIn, fpMatrix3_t *mOut);