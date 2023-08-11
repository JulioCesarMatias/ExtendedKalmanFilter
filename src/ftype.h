#pragma once

#include <Arduino.h>


#define acosF(x) acosf(x)
#define asinF(x) asinf(x)
#define sinF(x) sinf(x)
#define cosF(x) cosf(x)
#define tanF(x) tanf(x)
#define atanF(x) atanf(x)
#define atan2F(x,y) atan2f(x,y)
#define sqrtF(x) sqrtf(x)
#define fmaxF(x,y) fmaxf(x,y)
#define powF(x,y) powf(x,y)
#define logF(x) logf(x)
#define fabsF(x) fabsf(x)
#define ceilF(x) ceilf(x)
#define fminF(x,y) fminf(x,y)
#define fmodF(x,y) fmodf(x,y)
#define fabsF(x) fabsf(x)

#if MATH_CHECK_INDEXES
#define ZERO_FARRAY(a) a.zero()
#else
#define ZERO_FARRAY(a) memset(a, 0, sizeof(a))
#endif

/*
 * @brief: Check whether a float is zero
 */
inline bool is_zero(const float x) {
    return fabsf(x) < FLT_EPSILON;
}

/*
 * @brief: Check whether a double is zero
 */
inline bool is_zero(const double x) {
#ifdef ALLOW_DOUBLE_MATH_FUNCTIONS
  return fabs(x) < FLT_EPSILON;
#else
  return fabsf(static_cast<float>(x)) < FLT_EPSILON;
#endif
}
