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
void matrixIdentity(fpMatrix3_t *m);
void matrix_from_euler(fpMatrix3_t *m, float roll, float pitch, float yaw);
fpVector3_t multiply_matrix_by_vector(fpMatrix3_t m, fpVector3_t v);
void matrix_transpose3x3(fpMatrix3_t mIn, fpMatrix3_t *mOut);

// copied from https://code.google.com/p/cxutil/source/browse/include/cxutil/utility.h#70
#define _CHOOSE2(binoper, lexpr, lvar, rexpr, rvar) \
    (__extension__({                                \
        __typeof__(lexpr) lvar = (lexpr);           \
        __typeof__(rexpr) rvar = (rexpr);           \
        lvar binoper rvar ? lvar : rvar;            \
    }))
#define _CHOOSE_VAR2(prefix, unique) prefix##unique
#define _CHOOSE_VAR(prefix, unique) _CHOOSE_VAR2(prefix, unique)
#define _CHOOSE(binoper, lexpr, rexpr)          \
    _CHOOSE2(                                   \
        binoper,                                \
        lexpr, _CHOOSE_VAR(_left, __COUNTER__), \
        rexpr, _CHOOSE_VAR(_right, __COUNTER__))
#define MIN(a, b) _CHOOSE(<, a, b)
#define MAX(a, b) _CHOOSE(>, a, b)

// function to calculate the normalization (pythagoras) of a 2-dimensional vector
static inline float calc_length_pythagorean_2D(const float firstElement, const float secondElement)
{
    return sqrtf(sq(firstElement) + sq(secondElement));
}

// function to calculate the normalization (pythagoras) of a 3-dimensional vector
static inline float calc_length_pythagorean_3D(const float firstElement, const float secondElement, const float thirdElement)
{
    return sqrtf(sq(firstElement) + sq(secondElement) + sq(thirdElement));
}

static inline float constrainf(float amt, float low, float high)
{
    if (amt < low)
        return low;
    else if (amt > high)
        return high;
    else
        return amt;
}