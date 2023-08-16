#include "ekfMath.h"

static bool is_zero(const float x)
{
  return fabsf(x) < 1.19209290e-7F;
}

// initialise the quaternion to no rotation
void quaternion_initialise(fpQuaternion_t *q)
{
  q->q0 = 1.0f;
  q->q1 = q->q2 = q->q3 = 0.0f;
}

// create a quaternion from Euler angles
void quaternion_from_euler(fpQuaternion_t *q, float roll, float pitch, float yaw)
{
  const float cr2 = cos(roll * 0.5f);
  const float cp2 = cos(pitch * 0.5f);
  const float cy2 = cos(yaw * 0.5f);
  const float sr2 = sin(roll * 0.5f);
  const float sp2 = sin(pitch * 0.5f);
  const float sy2 = sin(yaw * 0.5f);

  q->q0 = cr2 * cp2 * cy2 + sr2 * sp2 * sy2;
  q->q1 = sr2 * cp2 * cy2 - cr2 * sp2 * sy2;
  q->q2 = cr2 * sp2 * cy2 + sr2 * cp2 * sy2;
  q->q3 = cr2 * cp2 * sy2 - sr2 * sp2 * cy2;
}

// create a quaternion from its axis-angle representation
void quaternion_from_axis_angle(fpVector3_t v, fpQuaternion_t *q)
{
  const float theta = sqrtf(sq(v.x) + sq(v.y) + sq(v.z));

  if (is_zero(theta))
  {
    q->q0 = 1.0f;
    q->q1 = q->q2 = q->q3 = 0.0f;
    return;
  }

  v.x /= theta;
  v.y /= theta;
  v.z /= theta;

  const float st2 = sin(0.5 * theta);

  q->q0 = cos(0.5 * theta);
  q->q1 = v.x * st2;
  q->q2 = v.y * st2;
  q->q3 = v.z * st2;
}

void quaternion_multiply(fpQuaternion_t qv, fpQuaternion_t *q)
{
  const float w1 = q->q0;
  const float x1 = q->q1;
  const float y1 = q->q2;
  const float z1 = q->q3;

  const float w2 = qv.q0;
  const float x2 = qv.q1;
  const float y2 = qv.q2;
  const float z2 = qv.q3;

  q->q0 = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;
  q->q1 = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2;
  q->q2 = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2;
  q->q3 = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2;
}

void quaternion_normalize(fpQuaternion_t *q)
{
  const float quatMag = sqrt(sq(q->q0) + sq(q->q1) + sq(q->q2) + sq(q->q3));

  if (!is_zero(quatMag))
  {
    const float quatMagInv = 1.0f / quatMag;
    q->q0 *= quatMagInv;
    q->q1 *= quatMagInv;
    q->q2 *= quatMagInv;
    q->q3 *= quatMagInv;
  }
}

// populate the supplied rotation matrix equivalent from this quaternion
void quaternion_rotation_matrix(fpQuaternion_t q, fpMatrix3_t *m)
{
  const float q3q3 = q.q2 * q.q2;
  const float q3q4 = q.q2 * q.q3;
  const float q2q2 = q.q1 * q.q1;
  const float q2q3 = q.q1 * q.q2;
  const float q2q4 = q.q1 * q.q3;
  const float q1q2 = q.q0 * q.q1;
  const float q1q3 = q.q0 * q.q2;
  const float q1q4 = q.q0 * q.q3;
  const float q4q4 = q.q3 * q.q3;

  m->m[0][0] = 1.0f - 2.0f * (q3q3 + q4q4);
  m->m[0][1] = 2.0f * (q2q3 - q1q4);
  m->m[0][2] = 2.0f * (q2q4 + q1q3);
  m->m[1][0] = 2.0f * (q2q3 + q1q4);
  m->m[1][1] = 1.0f - 2.0f * (q2q2 + q4q4);
  m->m[1][2] = 2.0f * (q3q4 - q1q2);
  m->m[2][0] = 2.0f * (q2q4 - q1q3);
  m->m[2][1] = 2.0f * (q3q4 + q1q2);
  m->m[2][2] = 1.0f - 2.0f * (q2q2 + q3q3);
}

void quaternion_to_euler(fpQuaternion_t q, float *roll, float *pitch, float *yaw)
{
  *roll = atan2f(2.0f * (q.q0 * q.q1 + q.q2 * q.q3), 1.0f - 2.0f * (q.q1 * q.q1 + q.q2 * q.q2));
  *pitch = asin(2.0f * (q.q0 * q.q2 - q.q3 * q.q1));
  *yaw = atan2f(2.0f * (q.q0 * q.q3 + q.q1 * q.q2), 1.0f - 2.0f * (q.q2 * q.q2 + q.q3 * q.q3));
}

void matrix_from_euler(fpMatrix3_t *m, float roll, float pitch, float yaw)
{
  const float cp = cos(pitch);
  const float sp = sin(pitch);
  const float sr = sin(roll);
  const float cr = cos(roll);
  const float sy = sin(yaw);
  const float cy = cos(yaw);

  m->m[0][0] = cp * cy;
  m->m[0][1] = (sr * sp * cy) - (cr * sy);
  m->m[0][2] = (cr * sp * cy) + (sr * sy);
  m->m[1][0] = cp * sy;
  m->m[1][1] = (sr * sp * sy) + (cr * cy);
  m->m[1][2] = (cr * sp * sy) - (sr * cy);
  m->m[2][0] = -sp;
  m->m[2][1] = sr * cp;
  m->m[2][2] = cr * cp;
}

// multiplication by a vector
fpVector3_t multiply_matrix_by_vector(fpMatrix3_t m, fpVector3_t v)
{
  fpVector3_t vRet;

  vRet.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
  vRet.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
  vRet.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;

  return vRet;
}

void matrix_transpose3x3(fpMatrix3_t mIn, fpMatrix3_t *mOut)
{
  fpMatrix3_t mCopy;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      mCopy.m[j][i] = mIn.m[i][j];
    }
  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      mOut->m[i][j] = mCopy.m[i][j];
    }
  }
}