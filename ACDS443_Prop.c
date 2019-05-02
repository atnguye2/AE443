/*
 * ACDS443_Prop.c
 *
 * Classroom License -- for classroom instructional use only.  Not for
 * government, commercial, academic research, or other organizational use.
 *
 * Code generation for model "ACDS443_Prop".
 *
 * Model version              : 1.25
 * Simulink Coder version : 9.0 (R2018b) 24-May-2018
 * C source code generated on : Wed May  1 19:40:39 2019
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "ACDS443_Prop.h"
#include "ACDS443_Prop_private.h"

/* Block signals (default storage) */
B_ACDS443_Prop_T ACDS443_Prop_B;

/* Continuous states */
X_ACDS443_Prop_T ACDS443_Prop_X;

/* Periodic continuous states */
PeriodicIndX_ACDS443_Prop_T ACDS443_Prop_PeriodicIndX;
PeriodicRngX_ACDS443_Prop_T ACDS443_Prop_PeriodicRngX;

/* Block states (default storage) */
DW_ACDS443_Prop_T ACDS443_Prop_DW;

/* Real-time model */
RT_MODEL_ACDS443_Prop_T ACDS443_Prop_M_;
RT_MODEL_ACDS443_Prop_T *const ACDS443_Prop_M = &ACDS443_Prop_M_;
static void rate_scheduler(void);
int32_T plook_s32dd_bincp(real_T u, const real_T bp[], uint32_T maxIndex, real_T
  *fraction, int32_T *prevIndex)
{
  int32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_s32d_prevIdx(u, bp, (uint32_T)*prevIndex, maxIndex);
    *fraction = (u - bp[(uint32_T)bpIndex]) / (bp[bpIndex + 1U] - bp[(uint32_T)
      bpIndex]);
  } else {
    bpIndex = (int32_T)(maxIndex - 1U);
    *fraction = 1.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

int32_T binsearch_s32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;
  uint32_T found;

  /* Binary Search using Previous Index */
  bpIdx = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx - 1U;
      bpIdx = (iRght + iLeft) >> 1U;
    } else if (u < bp[bpIdx + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIdx + 1U;
      bpIdx = (iRght + iLeft) >> 1U;
    }
  }

  return (int32_T)bpIdx;
}

/*
 *   This function updates active task flag for each subrate.
 * The function is called at model base rate, hence the
 * generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (ACDS443_Prop_M->Timing.TaskCounters.TID[2])++;
  if ((ACDS443_Prop_M->Timing.TaskCounters.TID[2]) > 79) {/* Sample time: [8.0s, 0.0s] */
    ACDS443_Prop_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/* State reduction function */
void local_stateReduction(real_T* x, int_T* p, int_T n, real_T* r)
{
  int_T i, j;
  for (i = 0, j = 0; i < n; ++i, ++j) {
    int_T k = p[i];
    real_T lb = r[j++];
    real_T xk = x[k]-lb;
    real_T rk = r[j]-lb;
    int_T q = (int_T) floor(xk/rk);
    if (q) {
      x[k] = xk-q*rk+lb;
    }
  }
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 28;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  ACDS443_Prop_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  ACDS443_Prop_step();
  ACDS443_Prop_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  ACDS443_Prop_step();
  ACDS443_Prop_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  local_stateReduction(x, rtsiGetPeriodicContStateIndices(si), 3,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T lo;
  uint32_T hi;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T y;
  real_T sr;
  real_T si;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

void rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(const real_T u0[3], const real_T u1[9],
  real_T y[3])
{
  real_T u1_0[9];
  int32_T ONE;
  int32_T THREE;
  int32_T r1;
  int32_T r2;
  real_T maxval;
  real_T a21;
  real_T x;
  real_T u0_idx_0;
  real_T u0_idx_1;
  real_T u0_idx_2;
  u0_idx_0 = u0[0];
  u0_idx_1 = u0[1];
  u0_idx_2 = u0[2];
  memcpy(&u1_0[0], &u1[0], 9U * sizeof(real_T));
  THREE = 2;
  r1 = 0;
  r2 = 1;
  x = u1_0[0];
  x = fabs(x);
  maxval = x;
  x = u1_0[1];
  x = fabs(x);
  a21 = x;
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  x = u1_0[2];
  x = fabs(x);
  a21 = x;
  if (a21 > maxval) {
    r1 = 2;
    r2 = 1;
    THREE = 0;
  }

  u1_0[r2] /= u1_0[r1];
  u1_0[THREE] /= u1_0[r1];
  u1_0[3 + r2] -= u1_0[3 + r1] * u1_0[r2];
  u1_0[3 + THREE] -= u1_0[3 + r1] * u1_0[THREE];
  u1_0[6 + r2] -= u1_0[6 + r1] * u1_0[r2];
  u1_0[6 + THREE] -= u1_0[6 + r1] * u1_0[THREE];
  x = u1_0[3 + THREE];
  x = fabs(x);
  a21 = x;
  x = u1_0[3 + r2];
  x = fabs(x);
  maxval = x;
  if (a21 > maxval) {
    ONE = r2 + 1;
    r2 = THREE;
    THREE = ONE - 1;
  }

  u1_0[3 + THREE] /= u1_0[3 + r2];
  u1_0[6 + THREE] -= u1_0[3 + THREE] * u1_0[6 + r2];
  y[r1] = u0_idx_0 / u1_0[r1];
  y[r2] = u0_idx_1 - u1_0[3 + r1] * y[r1];
  y[THREE] = u0_idx_2 - u1_0[6 + r1] * y[r1];
  y[r2] /= u1_0[3 + r2];
  y[THREE] -= u1_0[6 + r2] * y[r2];
  y[THREE] /= u1_0[6 + THREE];
  y[r2] -= u1_0[3 + THREE] * y[THREE];
  y[r1] -= y[THREE] * u1_0[THREE];
  y[r1] -= y[r2] * u1_0[r2];
}

/* Model step function */
void ACDS443_Prop_step(void)
{
  if (rtmIsMajorTimeStep(ACDS443_Prop_M)) {
    /* set solver stop time */
    if (!(ACDS443_Prop_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&ACDS443_Prop_M->solverInfo,
                            ((ACDS443_Prop_M->Timing.clockTickH0 + 1) *
        ACDS443_Prop_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&ACDS443_Prop_M->solverInfo,
                            ((ACDS443_Prop_M->Timing.clockTick0 + 1) *
        ACDS443_Prop_M->Timing.stepSize0 + ACDS443_Prop_M->Timing.clockTickH0 *
        ACDS443_Prop_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(ACDS443_Prop_M)) {
    ACDS443_Prop_M->Timing.t[0] = rtsiGetT(&ACDS443_Prop_M->solverInfo);
  }

  {
    real_T (*lastU)[3];
    real_T lastTime;
    real_T y[9];
    real_T u;
    real_T u_0;
    int32_T i;
    real_T u_idx_0;
    real_T u_idx_1;
    real_T tmp;
    real_T tmp_0;
    real_T tmp_1;
    real_T u0;
    int32_T u0_0;

    /* Integrator: '<S7>/phi theta psi' */
    ACDS443_Prop_B.phithetapsi[0] = ACDS443_Prop_X.phi[0];
    ACDS443_Prop_B.phithetapsi[1] = ACDS443_Prop_X.phi[1];
    ACDS443_Prop_B.phithetapsi[2] = ACDS443_Prop_X.phi[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* Integrator: '<S3>/p,q,r ' */
    ACDS443_Prop_B.pqr[0] = ACDS443_Prop_X.p[0];
    ACDS443_Prop_B.pqr[1] = ACDS443_Prop_X.p[1];
    ACDS443_Prop_B.pqr[2] = ACDS443_Prop_X.p[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* Integrator: '<S4>/Integrator' */
    ACDS443_Prop_B.Integrator = ACDS443_Prop_X.Integrator_CSTATE;

    /* Integrator: '<S5>/Integrator' */
    ACDS443_Prop_B.Integrator_a = ACDS443_Prop_X.Integrator_CSTATE_i;

    /* Integrator: '<S6>/Integrator' */
    ACDS443_Prop_B.Integrator_a1 = ACDS443_Prop_X.Integrator_CSTATE_k;
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_B.Integrator1[0] = ACDS443_Prop_X.Integrator1_CSTATE[0];

    /* Saturate: '<S6>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1[0];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentumSaturation[0] = u0;

    /* Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_B.Integrator1[1] = ACDS443_Prop_X.Integrator1_CSTATE[1];

    /* Saturate: '<S6>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1[1];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentumSaturation[1] = u0;

    /* Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_B.Integrator1[2] = ACDS443_Prop_X.Integrator1_CSTATE[2];

    /* Saturate: '<S6>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1[2];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentumSaturation[2] = u0;

    /* Derivative: '<S6>/Derivative' */
    if ((ACDS443_Prop_DW.TimeStampA >= ACDS443_Prop_M->Timing.t[0]) &&
        (ACDS443_Prop_DW.TimeStampB >= ACDS443_Prop_M->Timing.t[0])) {
      ACDS443_Prop_B.Derivative[0] = 0.0;
      ACDS443_Prop_B.Derivative[1] = 0.0;
      ACDS443_Prop_B.Derivative[2] = 0.0;
    } else {
      lastTime = ACDS443_Prop_DW.TimeStampA;
      lastU = &ACDS443_Prop_DW.LastUAtTimeA;
      if (ACDS443_Prop_DW.TimeStampA < ACDS443_Prop_DW.TimeStampB) {
        if (ACDS443_Prop_DW.TimeStampB < ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB;
        }
      } else {
        if (ACDS443_Prop_DW.TimeStampA >= ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB;
        }
      }

      lastTime = ACDS443_Prop_M->Timing.t[0] - lastTime;
      ACDS443_Prop_B.Derivative[0] = (ACDS443_Prop_B.MomentumSaturation[0] -
        (*lastU)[0]) / lastTime;
      ACDS443_Prop_B.Derivative[1] = (ACDS443_Prop_B.MomentumSaturation[1] -
        (*lastU)[1]) / lastTime;
      ACDS443_Prop_B.Derivative[2] = (ACDS443_Prop_B.MomentumSaturation[2] -
        (*lastU)[2]) / lastTime;
    }

    /* End of Derivative: '<S6>/Derivative' */

    /* Saturate: '<S6>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative[0];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation[0] = u0;

    /* Product: '<S6>/Multiply' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply[0] = ACDS443_Prop_P.MomentumWheel3_NormalVector[0] *
      ACDS443_Prop_B.MomentInputSaturation[0];

    /* Saturate: '<S6>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative[1];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation[1] = u0;

    /* Product: '<S6>/Multiply' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply[1] = ACDS443_Prop_P.MomentumWheel3_NormalVector[1] *
      ACDS443_Prop_B.MomentInputSaturation[1];

    /* Saturate: '<S6>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative[2];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSat) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSat;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSat) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSat;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation[2] = u0;

    /* Product: '<S6>/Multiply' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply[2] = ACDS443_Prop_P.MomentumWheel3_NormalVector[2] *
      ACDS443_Prop_B.MomentInputSaturation[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_B.Integrator1_b[0] = ACDS443_Prop_X.Integrator1_CSTATE_f[0];

    /* Saturate: '<S5>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_b[0];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_b) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_b;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_l) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_l;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_c[0] = u0;

    /* Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_B.Integrator1_b[1] = ACDS443_Prop_X.Integrator1_CSTATE_f[1];

    /* Saturate: '<S5>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_b[1];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_b) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_b;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_l) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_l;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_c[1] = u0;

    /* Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_B.Integrator1_b[2] = ACDS443_Prop_X.Integrator1_CSTATE_f[2];

    /* Saturate: '<S5>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_b[2];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_b) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_b;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_l) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_l;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_c[2] = u0;

    /* Derivative: '<S5>/Derivative' */
    if ((ACDS443_Prop_DW.TimeStampA_i >= ACDS443_Prop_M->Timing.t[0]) &&
        (ACDS443_Prop_DW.TimeStampB_j >= ACDS443_Prop_M->Timing.t[0])) {
      ACDS443_Prop_B.Derivative_p[0] = 0.0;
      ACDS443_Prop_B.Derivative_p[1] = 0.0;
      ACDS443_Prop_B.Derivative_p[2] = 0.0;
    } else {
      lastTime = ACDS443_Prop_DW.TimeStampA_i;
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_f;
      if (ACDS443_Prop_DW.TimeStampA_i < ACDS443_Prop_DW.TimeStampB_j) {
        if (ACDS443_Prop_DW.TimeStampB_j < ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB_j;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB_a;
        }
      } else {
        if (ACDS443_Prop_DW.TimeStampA_i >= ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB_j;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB_a;
        }
      }

      lastTime = ACDS443_Prop_M->Timing.t[0] - lastTime;
      ACDS443_Prop_B.Derivative_p[0] = (ACDS443_Prop_B.MomentumSaturation_c[0] -
        (*lastU)[0]) / lastTime;
      ACDS443_Prop_B.Derivative_p[1] = (ACDS443_Prop_B.MomentumSaturation_c[1] -
        (*lastU)[1]) / lastTime;
      ACDS443_Prop_B.Derivative_p[2] = (ACDS443_Prop_B.MomentumSaturation_c[2] -
        (*lastU)[2]) / lastTime;
    }

    /* End of Derivative: '<S5>/Derivative' */

    /* Saturate: '<S5>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_p[0];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_e) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_e;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_o) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_o;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_o[0] = u0;

    /* Product: '<S5>/Multiply' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b[0] = ACDS443_Prop_P.MomentumWheel2_NormalVector[0]
      * ACDS443_Prop_B.MomentInputSaturation_o[0];

    /* Saturate: '<S5>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_p[1];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_e) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_e;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_o) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_o;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_o[1] = u0;

    /* Product: '<S5>/Multiply' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b[1] = ACDS443_Prop_P.MomentumWheel2_NormalVector[1]
      * ACDS443_Prop_B.MomentInputSaturation_o[1];

    /* Saturate: '<S5>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_p[2];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_e) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_e;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_o) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_o;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_o[2] = u0;

    /* Product: '<S5>/Multiply' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b[2] = ACDS443_Prop_P.MomentumWheel2_NormalVector[2]
      * ACDS443_Prop_B.MomentInputSaturation_o[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_B.Integrator1_f[0] = ACDS443_Prop_X.Integrator1_CSTATE_a[0];

    /* Saturate: '<S4>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_f[0];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_g) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_g;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_m) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_m;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_h[0] = u0;

    /* Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_B.Integrator1_f[1] = ACDS443_Prop_X.Integrator1_CSTATE_a[1];

    /* Saturate: '<S4>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_f[1];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_g) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_g;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_m) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_m;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_h[1] = u0;

    /* Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_B.Integrator1_f[2] = ACDS443_Prop_X.Integrator1_CSTATE_a[2];

    /* Saturate: '<S4>/Momentum Saturation' */
    u0 = ACDS443_Prop_B.Integrator1_f[2];
    if (u0 > ACDS443_Prop_P.MomentumSaturation_UpperSat_g) {
      u0 = ACDS443_Prop_P.MomentumSaturation_UpperSat_g;
    } else {
      if (u0 < ACDS443_Prop_P.MomentumSaturation_LowerSat_m) {
        u0 = ACDS443_Prop_P.MomentumSaturation_LowerSat_m;
      }
    }

    ACDS443_Prop_B.MomentumSaturation_h[2] = u0;

    /* Derivative: '<S4>/Derivative' */
    if ((ACDS443_Prop_DW.TimeStampA_j >= ACDS443_Prop_M->Timing.t[0]) &&
        (ACDS443_Prop_DW.TimeStampB_m >= ACDS443_Prop_M->Timing.t[0])) {
      ACDS443_Prop_B.Derivative_m[0] = 0.0;
      ACDS443_Prop_B.Derivative_m[1] = 0.0;
      ACDS443_Prop_B.Derivative_m[2] = 0.0;
    } else {
      lastTime = ACDS443_Prop_DW.TimeStampA_j;
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_b;
      if (ACDS443_Prop_DW.TimeStampA_j < ACDS443_Prop_DW.TimeStampB_m) {
        if (ACDS443_Prop_DW.TimeStampB_m < ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB_m;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB_c;
        }
      } else {
        if (ACDS443_Prop_DW.TimeStampA_j >= ACDS443_Prop_M->Timing.t[0]) {
          lastTime = ACDS443_Prop_DW.TimeStampB_m;
          lastU = &ACDS443_Prop_DW.LastUAtTimeB_c;
        }
      }

      lastTime = ACDS443_Prop_M->Timing.t[0] - lastTime;
      ACDS443_Prop_B.Derivative_m[0] = (ACDS443_Prop_B.MomentumSaturation_h[0] -
        (*lastU)[0]) / lastTime;
      ACDS443_Prop_B.Derivative_m[1] = (ACDS443_Prop_B.MomentumSaturation_h[1] -
        (*lastU)[1]) / lastTime;
      ACDS443_Prop_B.Derivative_m[2] = (ACDS443_Prop_B.MomentumSaturation_h[2] -
        (*lastU)[2]) / lastTime;
    }

    /* End of Derivative: '<S4>/Derivative' */

    /* Saturate: '<S4>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_m[0];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_h) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_h;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_a) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_a;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_g[0] = u0;

    /* Product: '<S4>/Multiply' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b4[0] = ACDS443_Prop_P.MomentumWheel1_NormalVector[0]
      * ACDS443_Prop_B.MomentInputSaturation_g[0];

    /* Saturate: '<S4>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_m[1];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_h) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_h;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_a) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_a;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_g[1] = u0;

    /* Product: '<S4>/Multiply' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b4[1] = ACDS443_Prop_P.MomentumWheel1_NormalVector[1]
      * ACDS443_Prop_B.MomentInputSaturation_g[1];

    /* Saturate: '<S4>/Moment Input Saturation' */
    u0 = ACDS443_Prop_B.Derivative_m[2];
    if (u0 > ACDS443_Prop_P.MomentInputSaturation_UpperSa_h) {
      u0 = ACDS443_Prop_P.MomentInputSaturation_UpperSa_h;
    } else {
      if (u0 < ACDS443_Prop_P.MomentInputSaturation_LowerSa_a) {
        u0 = ACDS443_Prop_P.MomentInputSaturation_LowerSa_a;
      }
    }

    ACDS443_Prop_B.MomentInputSaturation_g[2] = u0;

    /* Product: '<S4>/Multiply' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    ACDS443_Prop_B.Multiply_b4[2] = ACDS443_Prop_P.MomentumWheel1_NormalVector[2]
      * ACDS443_Prop_B.MomentInputSaturation_g[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
    }

    /* RelationalOperator: '<S1>/Less Than' incorporates:
     *  Constant: '<S1>/Saturation1'
     */
    ACDS443_Prop_B.LessThan[0] = (ACDS443_Prop_B.Integrator <
      ACDS443_Prop_P.Saturation1_Value);
    ACDS443_Prop_B.LessThan[1] = (ACDS443_Prop_B.Integrator_a <
      ACDS443_Prop_P.Saturation1_Value);
    ACDS443_Prop_B.LessThan[2] = (ACDS443_Prop_B.Integrator_a1 <
      ACDS443_Prop_P.Saturation1_Value);

    /* RelationalOperator: '<S1>/GreaterThan' incorporates:
     *  Constant: '<S1>/Saturation'
     */
    ACDS443_Prop_B.GreaterThan[0] = (ACDS443_Prop_B.Integrator >
      ACDS443_Prop_P.Saturation_Value);
    ACDS443_Prop_B.GreaterThan[1] = (ACDS443_Prop_B.Integrator_a >
      ACDS443_Prop_P.Saturation_Value);
    ACDS443_Prop_B.GreaterThan[2] = (ACDS443_Prop_B.Integrator_a1 >
      ACDS443_Prop_P.Saturation_Value);

    /* ZeroOrderHold: '<S1>/Zero-Order Hold1' incorporates:
     *  ZeroOrderHold: '<S1>/Zero-Order Hold'
     */
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[2] == 0) {
      ACDS443_Prop_B.ZeroOrderHold1[0] = ACDS443_Prop_B.LessThan[0];

      /* Product: '<S1>/Multiply1' incorporates:
       *  Constant: '<S1>/Constant3'
       */
      ACDS443_Prop_B.Multiply1[0] = (real_T)ACDS443_Prop_B.ZeroOrderHold1[0] *
        ACDS443_Prop_P.Constant3_Value[0];
      ACDS443_Prop_B.ZeroOrderHold1[1] = ACDS443_Prop_B.LessThan[1];

      /* Product: '<S1>/Multiply1' incorporates:
       *  Constant: '<S1>/Constant3'
       */
      ACDS443_Prop_B.Multiply1[1] = (real_T)ACDS443_Prop_B.ZeroOrderHold1[1] *
        ACDS443_Prop_P.Constant3_Value[1];
      ACDS443_Prop_B.ZeroOrderHold1[2] = ACDS443_Prop_B.LessThan[2];

      /* Product: '<S1>/Multiply1' incorporates:
       *  Constant: '<S1>/Constant3'
       */
      ACDS443_Prop_B.Multiply1[2] = (real_T)ACDS443_Prop_B.ZeroOrderHold1[2] *
        ACDS443_Prop_P.Constant3_Value[2];
      ACDS443_Prop_B.ZeroOrderHold[0] = ACDS443_Prop_B.GreaterThan[0];

      /* Product: '<S1>/Multiply' incorporates:
       *  Constant: '<S1>/Constant2'
       */
      ACDS443_Prop_B.Multiply_n[0] = (real_T)ACDS443_Prop_B.ZeroOrderHold[0] *
        ACDS443_Prop_P.Constant2_Value[0];

      /* Sum: '<Root>/Add1' */
      ACDS443_Prop_B.Add1[0] = ACDS443_Prop_B.Multiply1[0] +
        ACDS443_Prop_B.Multiply_n[0];
      ACDS443_Prop_B.ZeroOrderHold[1] = ACDS443_Prop_B.GreaterThan[1];

      /* Product: '<S1>/Multiply' incorporates:
       *  Constant: '<S1>/Constant2'
       */
      ACDS443_Prop_B.Multiply_n[1] = (real_T)ACDS443_Prop_B.ZeroOrderHold[1] *
        ACDS443_Prop_P.Constant2_Value[1];

      /* Sum: '<Root>/Add1' */
      ACDS443_Prop_B.Add1[1] = ACDS443_Prop_B.Multiply1[1] +
        ACDS443_Prop_B.Multiply_n[1];
      ACDS443_Prop_B.ZeroOrderHold[2] = ACDS443_Prop_B.GreaterThan[2];

      /* Product: '<S1>/Multiply' incorporates:
       *  Constant: '<S1>/Constant2'
       */
      ACDS443_Prop_B.Multiply_n[2] = (real_T)ACDS443_Prop_B.ZeroOrderHold[2] *
        ACDS443_Prop_P.Constant2_Value[2];

      /* Sum: '<Root>/Add1' */
      ACDS443_Prop_B.Add1[2] = ACDS443_Prop_B.Multiply1[2] +
        ACDS443_Prop_B.Multiply_n[2];
    }

    /* End of ZeroOrderHold: '<S1>/Zero-Order Hold1' */

    /* Sum: '<Root>/Add' */
    ACDS443_Prop_B.Add[0] = ((ACDS443_Prop_B.Multiply[0] - ACDS443_Prop_B.Add1[0])
      + ACDS443_Prop_B.Multiply_b[0]) + ACDS443_Prop_B.Multiply_b4[0];
    ACDS443_Prop_B.Add[1] = ((ACDS443_Prop_B.Multiply[1] - ACDS443_Prop_B.Add1[1])
      + ACDS443_Prop_B.Multiply_b[1]) + ACDS443_Prop_B.Multiply_b4[1];
    ACDS443_Prop_B.Add[2] = ((ACDS443_Prop_B.Multiply[2] - ACDS443_Prop_B.Add1[2])
      + ACDS443_Prop_B.Multiply_b[2]) + ACDS443_Prop_B.Multiply_b4[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
      /* RandomNumber: '<Root>/Aerodynamic Drag' */
      u0 = ACDS443_Prop_DW.NextOutput[0];
      ACDS443_Prop_B.AerodynamicDrag[0] = u0;
      u0 = ACDS443_Prop_DW.NextOutput[1];
      ACDS443_Prop_B.AerodynamicDrag[1] = u0;
      u0 = ACDS443_Prop_DW.NextOutput[2];
      ACDS443_Prop_B.AerodynamicDrag[2] = u0;
    }

    /* Gain: '<Root>/Gain' */
    ACDS443_Prop_B.Gain[0] = ACDS443_Prop_P.Gain_Gain *
      ACDS443_Prop_B.phithetapsi[0];

    /* Gain: '<Root>/Gain1' */
    ACDS443_Prop_B.Gain1[0] = ACDS443_Prop_P.Gain1_Gain * ACDS443_Prop_B.pqr[0];

    /* Integrator: '<Root>/Integrator' */
    ACDS443_Prop_B.Integrator_g[0] = ACDS443_Prop_X.Integrator_CSTATE_j[0];

    /* Gain: '<Root>/Gain2' */
    ACDS443_Prop_B.Gain2[0] = ACDS443_Prop_P.Gain2_Gain *
      ACDS443_Prop_B.Integrator_g[0];

    /* Gain: '<Root>/Gain' */
    ACDS443_Prop_B.Gain[1] = ACDS443_Prop_P.Gain_Gain *
      ACDS443_Prop_B.phithetapsi[1];

    /* Gain: '<Root>/Gain1' */
    ACDS443_Prop_B.Gain1[1] = ACDS443_Prop_P.Gain1_Gain * ACDS443_Prop_B.pqr[1];

    /* Integrator: '<Root>/Integrator' */
    ACDS443_Prop_B.Integrator_g[1] = ACDS443_Prop_X.Integrator_CSTATE_j[1];

    /* Gain: '<Root>/Gain2' */
    ACDS443_Prop_B.Gain2[1] = ACDS443_Prop_P.Gain2_Gain *
      ACDS443_Prop_B.Integrator_g[1];

    /* Gain: '<Root>/Gain' */
    ACDS443_Prop_B.Gain[2] = ACDS443_Prop_P.Gain_Gain *
      ACDS443_Prop_B.phithetapsi[2];

    /* Gain: '<Root>/Gain1' */
    ACDS443_Prop_B.Gain1[2] = ACDS443_Prop_P.Gain1_Gain * ACDS443_Prop_B.pqr[2];

    /* Integrator: '<Root>/Integrator' */
    ACDS443_Prop_B.Integrator_g[2] = ACDS443_Prop_X.Integrator_CSTATE_j[2];

    /* Gain: '<Root>/Gain2' */
    ACDS443_Prop_B.Gain2[2] = ACDS443_Prop_P.Gain2_Gain *
      ACDS443_Prop_B.Integrator_g[2];
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[2] == 0) {
      /* DataTypeConversion: '<S1>/Data Type Conversion' */
      ACDS443_Prop_B.DataTypeConversion[0] = ACDS443_Prop_B.ZeroOrderHold[0];

      /* DataTypeConversion: '<S1>/Data Type Conversion1' */
      ACDS443_Prop_B.DataTypeConversion1[0] = ACDS443_Prop_B.ZeroOrderHold1[0];

      /* Sum: '<S1>/Sum2' */
      ACDS443_Prop_B.Sum2[0] = ACDS443_Prop_B.DataTypeConversion[0] +
        ACDS443_Prop_B.DataTypeConversion1[0];

      /* DataTypeConversion: '<S1>/Data Type Conversion' */
      ACDS443_Prop_B.DataTypeConversion[1] = ACDS443_Prop_B.ZeroOrderHold[1];

      /* DataTypeConversion: '<S1>/Data Type Conversion1' */
      ACDS443_Prop_B.DataTypeConversion1[1] = ACDS443_Prop_B.ZeroOrderHold1[1];

      /* Sum: '<S1>/Sum2' */
      ACDS443_Prop_B.Sum2[1] = ACDS443_Prop_B.DataTypeConversion[1] +
        ACDS443_Prop_B.DataTypeConversion1[1];

      /* DataTypeConversion: '<S1>/Data Type Conversion' */
      ACDS443_Prop_B.DataTypeConversion[2] = ACDS443_Prop_B.ZeroOrderHold[2];

      /* DataTypeConversion: '<S1>/Data Type Conversion1' */
      ACDS443_Prop_B.DataTypeConversion1[2] = ACDS443_Prop_B.ZeroOrderHold1[2];

      /* Sum: '<S1>/Sum2' */
      ACDS443_Prop_B.Sum2[2] = ACDS443_Prop_B.DataTypeConversion[2] +
        ACDS443_Prop_B.DataTypeConversion1[2];
    }

    /* DotProduct: '<S4>/Dot Product' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    u0 = ACDS443_Prop_P.MomentumWheel1_NormalVector[0];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_g[0];
    u = u0 * tmp_1;

    /* DotProduct: '<S5>/Dot Product' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    u0 = ACDS443_Prop_P.MomentumWheel2_NormalVector[0];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_o[0];

    /* DotProduct: '<S4>/Dot Product' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    u_idx_0 = u0;
    tmp = tmp_1;
    u0 = ACDS443_Prop_P.MomentumWheel1_NormalVector[1];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_g[1];
    u += u0 * tmp_1;

    /* DotProduct: '<S5>/Dot Product' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    u0 = ACDS443_Prop_P.MomentumWheel2_NormalVector[1];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_o[1];

    /* DotProduct: '<S4>/Dot Product' incorporates:
     *  Constant: '<S4>/Wheel Orientation'
     */
    u_idx_1 = u0;
    tmp_0 = tmp_1;
    u0 = ACDS443_Prop_P.MomentumWheel1_NormalVector[2];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_g[2];
    u += u0 * tmp_1;

    /* DotProduct: '<S5>/Dot Product' incorporates:
     *  Constant: '<S5>/Wheel Orientation'
     */
    u0 = ACDS443_Prop_P.MomentumWheel2_NormalVector[2];
    tmp_1 = ACDS443_Prop_B.MomentInputSaturation_o[2];

    /* DotProduct: '<S4>/Dot Product' */
    ACDS443_Prop_B.DotProduct = u;

    /* DotProduct: '<S5>/Dot Product' */
    lastTime = tmp;
    u = u_idx_0;
    u_0 = u * lastTime;

    /* DotProduct: '<S6>/Dot Product' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    u = ACDS443_Prop_P.MomentumWheel3_NormalVector[0];
    lastTime = ACDS443_Prop_B.MomentInputSaturation[0];
    u_idx_0 = u;
    tmp = lastTime;

    /* DotProduct: '<S5>/Dot Product' */
    lastTime = tmp_0;
    u = u_idx_1;
    u_0 += u * lastTime;

    /* DotProduct: '<S6>/Dot Product' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    u = ACDS443_Prop_P.MomentumWheel3_NormalVector[1];
    lastTime = ACDS443_Prop_B.MomentInputSaturation[1];
    u_idx_1 = u;
    tmp_0 = lastTime;

    /* DotProduct: '<S5>/Dot Product' */
    lastTime = tmp_1;
    u = u0;
    u_0 += u * lastTime;

    /* DotProduct: '<S6>/Dot Product' incorporates:
     *  Constant: '<S6>/Wheel Orientation'
     */
    u = ACDS443_Prop_P.MomentumWheel3_NormalVector[2];
    lastTime = ACDS443_Prop_B.MomentInputSaturation[2];
    u0 = u;
    tmp_1 = lastTime;

    /* DotProduct: '<S5>/Dot Product' */
    ACDS443_Prop_B.DotProduct_j = u_0;

    /* DotProduct: '<S6>/Dot Product' */
    u = u_idx_0 * tmp;
    u += u_idx_1 * tmp_0;
    u += u0 * tmp_1;
    ACDS443_Prop_B.DotProduct_f = u;
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
      /* RandomNumber: '<Root>/SRP and Gravity ' */
      u0 = ACDS443_Prop_DW.NextOutput_k[0];
      ACDS443_Prop_B.SRPandGravity[0] = u0;
      u0 = ACDS443_Prop_DW.NextOutput_k[1];
      ACDS443_Prop_B.SRPandGravity[1] = u0;
      u0 = ACDS443_Prop_DW.NextOutput_k[2];
      ACDS443_Prop_B.SRPandGravity[2] = u0;
    }

    /* SignalConversion: '<S15>/TmpSignal ConversionAtsincosInport1' */
    ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[0] =
      ACDS443_Prop_B.phithetapsi[2];
    ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[1] =
      ACDS443_Prop_B.phithetapsi[1];
    ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[2] =
      ACDS443_Prop_B.phithetapsi[0];

    /* Trigonometry: '<S15>/sincos' */
    u0 = ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[0];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1[0] = lastTime;
    ACDS443_Prop_B.sincos_o2[0] = u0;
    u0 = ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[1];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1[1] = lastTime;
    ACDS443_Prop_B.sincos_o2[1] = u0;
    u0 = ACDS443_Prop_B.TmpSignalConversionAtsincosInpo[2];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1[2] = lastTime;
    ACDS443_Prop_B.sincos_o2[2] = u0;

    /* Fcn: '<S15>/Fcn11' */
    ACDS443_Prop_B.VectorConcatenate[0] = ACDS443_Prop_B.sincos_o2[1] *
      ACDS443_Prop_B.sincos_o2[0];

    /* Fcn: '<S15>/Fcn21' */
    ACDS443_Prop_B.VectorConcatenate[1] = ACDS443_Prop_B.sincos_o1[2] *
      ACDS443_Prop_B.sincos_o1[1] * ACDS443_Prop_B.sincos_o2[0] -
      ACDS443_Prop_B.sincos_o2[2] * ACDS443_Prop_B.sincos_o1[0];

    /* Fcn: '<S15>/Fcn31' */
    ACDS443_Prop_B.VectorConcatenate[2] = ACDS443_Prop_B.sincos_o2[2] *
      ACDS443_Prop_B.sincos_o1[1] * ACDS443_Prop_B.sincos_o2[0] +
      ACDS443_Prop_B.sincos_o1[2] * ACDS443_Prop_B.sincos_o1[0];

    /* Fcn: '<S15>/Fcn12' */
    ACDS443_Prop_B.VectorConcatenate[3] = ACDS443_Prop_B.sincos_o2[1] *
      ACDS443_Prop_B.sincos_o1[0];

    /* Fcn: '<S15>/Fcn22' */
    ACDS443_Prop_B.VectorConcatenate[4] = ACDS443_Prop_B.sincos_o1[2] *
      ACDS443_Prop_B.sincos_o1[1] * ACDS443_Prop_B.sincos_o1[0] +
      ACDS443_Prop_B.sincos_o2[2] * ACDS443_Prop_B.sincos_o2[0];

    /* Fcn: '<S15>/Fcn32' */
    ACDS443_Prop_B.VectorConcatenate[5] = ACDS443_Prop_B.sincos_o2[2] *
      ACDS443_Prop_B.sincos_o1[1] * ACDS443_Prop_B.sincos_o1[0] -
      ACDS443_Prop_B.sincos_o1[2] * ACDS443_Prop_B.sincos_o2[0];

    /* Fcn: '<S15>/Fcn13' */
    ACDS443_Prop_B.VectorConcatenate[6] = -ACDS443_Prop_B.sincos_o1[1];

    /* Fcn: '<S15>/Fcn23' */
    ACDS443_Prop_B.VectorConcatenate[7] = ACDS443_Prop_B.sincos_o1[2] *
      ACDS443_Prop_B.sincos_o2[1];

    /* Fcn: '<S15>/Fcn33' */
    ACDS443_Prop_B.VectorConcatenate[8] = ACDS443_Prop_B.sincos_o2[2] *
      ACDS443_Prop_B.sincos_o2[1];

    /* Trigonometry: '<S16>/sincos' */
    u0 = ACDS443_Prop_B.phithetapsi[0];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1_m[0] = lastTime;
    ACDS443_Prop_B.sincos_o2_i[0] = u0;
    u0 = ACDS443_Prop_B.phithetapsi[1];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1_m[1] = lastTime;
    ACDS443_Prop_B.sincos_o2_i[1] = u0;
    u0 = ACDS443_Prop_B.phithetapsi[2];
    lastTime = sin(u0);
    u0 = cos(u0);
    ACDS443_Prop_B.sincos_o1_m[2] = lastTime;
    ACDS443_Prop_B.sincos_o2_i[2] = u0;

    /* Fcn: '<S16>/phidot' */
    ACDS443_Prop_B.phidot = (ACDS443_Prop_B.pqr[1] * ACDS443_Prop_B.sincos_o1_m
      [0] + ACDS443_Prop_B.pqr[2] * ACDS443_Prop_B.sincos_o2_i[0]) *
      (ACDS443_Prop_B.sincos_o1_m[1] / ACDS443_Prop_B.sincos_o2_i[1]) +
      ACDS443_Prop_B.pqr[0];

    /* Fcn: '<S16>/thetadot' */
    ACDS443_Prop_B.thetadot = ACDS443_Prop_B.pqr[1] *
      ACDS443_Prop_B.sincos_o2_i[0] - ACDS443_Prop_B.pqr[2] *
      ACDS443_Prop_B.sincos_o1_m[0];

    /* Fcn: '<S16>/psidot' */
    ACDS443_Prop_B.psidot = (ACDS443_Prop_B.pqr[1] * ACDS443_Prop_B.sincos_o1_m
      [0] + ACDS443_Prop_B.pqr[2] * ACDS443_Prop_B.sincos_o2_i[0]) /
      ACDS443_Prop_B.sincos_o2_i[1];

    /* SignalConversion: '<S7>/TmpSignal ConversionAtphi theta psiInport1' */
    ACDS443_Prop_B.TmpSignalConversionAtphithetaps[0] = ACDS443_Prop_B.phidot;
    ACDS443_Prop_B.TmpSignalConversionAtphithetaps[1] = ACDS443_Prop_B.thetadot;
    ACDS443_Prop_B.TmpSignalConversionAtphithetaps[2] = ACDS443_Prop_B.psidot;

    /* Integrator: '<S9>/mass  ' */
    /* Limited  Integrator  (w/ Saturation Port) */
    if (ACDS443_Prop_X.mass >= ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_i) {
      ACDS443_Prop_X.mass = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_i;
      ACDS443_Prop_B.mass_o2 = 1.0;
    } else if (ACDS443_Prop_X.mass <=
               ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_l) {
      ACDS443_Prop_X.mass = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_l;
      ACDS443_Prop_B.mass_o2 = -1.0;
    } else {
      ACDS443_Prop_B.mass_o2 = 0.0;
    }

    ACDS443_Prop_B.mass_o1 = ACDS443_Prop_X.mass;

    /* End of Integrator: '<S9>/mass  ' */

    /* SignalConversion: '<S23>/TmpSignal ConversionAtPreLook-Up Index SearchInport2' incorporates:
     *  Constant: '<S23>/Constant'
     *  Constant: '<S23>/Constant1'
     */
    ACDS443_Prop_B.TmpSignalConversionAtPreLookUpI[0] =
      ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_l;
    ACDS443_Prop_B.TmpSignalConversionAtPreLookUpI[1] =
      ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_i;

    /* PreLookup: '<S23>/PreLook-Up Index Search' */
    ACDS443_Prop_B.PreLookUpIndexSearch_o1 = plook_s32dd_bincp
      (ACDS443_Prop_B.mass_o1, ACDS443_Prop_B.TmpSignalConversionAtPreLookUpI,
       1U, &ACDS443_Prop_B.PreLookUpIndexSearch_o2,
       &ACDS443_Prop_DW.PreLookUpIndexSearch_DWORK1);

    /* LookupNDDirect: '<S25>/[0]'
     *
     * About '<S25>/[0]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix,
     *  which is contiguous for column-major array
     */
    u0_0 = ACDS443_Prop_B.PreLookUpIndexSearch_o1;
    if (u0_0 > 1) {
      u0_0 = 1;
    } else {
      if (u0_0 < 0) {
        u0_0 = 0;
      }
    }

    u0_0 *= 9;

    /* Fcn: '<S25>/Fcn' */
    ACDS443_Prop_B.Fcn = 1.0 - ACDS443_Prop_B.PreLookUpIndexSearch_o2;
    for (i = 0; i < 9; i++) {
      /* LookupNDDirect: '<S25>/[0]'
       *
       * About '<S25>/[0]':
       *  3-dimensional Direct Look-Up returning a 2-D Matrix,
       *  which is contiguous for column-major array
       */
      ACDS443_Prop_B.u[i] = ACDS443_Prop_P.InterpolateInertia_matrix[i + u0_0];

      /* Product: '<S25>/1-lambda_x' */
      ACDS443_Prop_B.ulambda_x[i] = ACDS443_Prop_B.u[i] * ACDS443_Prop_B.Fcn;
    }

    /* Sum: '<S25>/Sum' incorporates:
     *  Constant: '<S25>/offset for upper index 1'
     */
    u0_0 = ACDS443_Prop_B.PreLookUpIndexSearch_o1;
    i = ACDS443_Prop_P.offsetforupperindex1_Value;
    if ((u0_0 < 0) && (i < MIN_int32_T - u0_0)) {
      u0_0 = MIN_int32_T;
    } else if ((u0_0 > 0) && (i > MAX_int32_T - u0_0)) {
      u0_0 = MAX_int32_T;
    } else {
      u0_0 += i;
    }

    ACDS443_Prop_B.Sum_dq = u0_0;

    /* End of Sum: '<S25>/Sum' */

    /* LookupNDDirect: '<S25>/[1]'
     *
     * About '<S25>/[1]':
     *  3-dimensional Direct Look-Up returning a 2-D Matrix,
     *  which is contiguous for column-major array
     */
    u0_0 = ACDS443_Prop_B.Sum_dq;
    if (u0_0 > 1) {
      u0_0 = 1;
    } else {
      if (u0_0 < 0) {
        u0_0 = 0;
      }
    }

    u0_0 *= 9;

    /* Sum: '<Root>/Sum3' */
    ACDS443_Prop_B.Sum3_p = (ACDS443_Prop_B.Gain2[0] + ACDS443_Prop_B.Gain2[1])
      + ACDS443_Prop_B.Gain2[2];

    /* Product: '<S23>/Product' incorporates:
     *  Constant: '<S23>/slope'
     */
    lastTime = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_i -
      ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_l;
    for (i = 0; i < 9; i++) {
      /* LookupNDDirect: '<S25>/[1]'
       *
       * About '<S25>/[1]':
       *  3-dimensional Direct Look-Up returning a 2-D Matrix,
       *  which is contiguous for column-major array
       */
      ACDS443_Prop_B.u_o[i] = ACDS443_Prop_P.InterpolateInertia_matrix[i + u0_0];

      /* Product: '<S25>/lambda_x.' */
      ACDS443_Prop_B.lambda_x[i] = ACDS443_Prop_B.u_o[i] *
        ACDS443_Prop_B.PreLookUpIndexSearch_o2;

      /* Sum: '<S25>/Sum3' */
      ACDS443_Prop_B.Sum3[i] = ACDS443_Prop_B.ulambda_x[i] +
        ACDS443_Prop_B.lambda_x[i];

      /* Product: '<S23>/Product' incorporates:
       *  Constant: '<S23>/slope'
       */
      u0 = (ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_j[i] -
            ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_d[i]) / lastTime;
      ACDS443_Prop_B.Product[i] = u0 * ACDS443_Prop_B.Sum3_p;
    }

    for (u0_0 = 0; u0_0 < 3; u0_0++) {
      /* Concatenate: '<S9>/Matrix Concatenation' */
      ACDS443_Prop_B.MatrixConcatenation[6 * u0_0] = ACDS443_Prop_B.Sum3[3 *
        u0_0];
      ACDS443_Prop_B.MatrixConcatenation[3 + 6 * u0_0] = ACDS443_Prop_B.Product
        [3 * u0_0];

      /* Selector: '<S8>/Selector' */
      ACDS443_Prop_B.Selector[3 * u0_0] = ACDS443_Prop_B.MatrixConcatenation[6 *
        u0_0];

      /* Concatenate: '<S9>/Matrix Concatenation' */
      ACDS443_Prop_B.MatrixConcatenation[1 + 6 * u0_0] = ACDS443_Prop_B.Sum3[3 *
        u0_0 + 1];
      ACDS443_Prop_B.MatrixConcatenation[4 + 6 * u0_0] = ACDS443_Prop_B.Product
        [3 * u0_0 + 1];

      /* Selector: '<S8>/Selector' */
      ACDS443_Prop_B.Selector[1 + 3 * u0_0] =
        ACDS443_Prop_B.MatrixConcatenation[6 * u0_0 + 1];

      /* Concatenate: '<S9>/Matrix Concatenation' */
      ACDS443_Prop_B.MatrixConcatenation[2 + 6 * u0_0] = ACDS443_Prop_B.Sum3[3 *
        u0_0 + 2];
      ACDS443_Prop_B.MatrixConcatenation[5 + 6 * u0_0] = ACDS443_Prop_B.Product
        [3 * u0_0 + 2];

      /* Selector: '<S8>/Selector' */
      ACDS443_Prop_B.Selector[2 + 3 * u0_0] =
        ACDS443_Prop_B.MatrixConcatenation[6 * u0_0 + 2];
    }

    /* Product: '<S19>/Product' */
    memcpy(&y[0], &ACDS443_Prop_B.Selector[0], 9U * sizeof(real_T));
    u_idx_0 = ACDS443_Prop_B.pqr[0];
    u_idx_1 = ACDS443_Prop_B.pqr[1];
    u0 = ACDS443_Prop_B.pqr[2];
    for (u0_0 = 0; u0_0 < 3; u0_0++) {
      ACDS443_Prop_B.Product_h[u0_0] = 0.0;
      ACDS443_Prop_B.Product_h[u0_0] += y[u0_0] * u_idx_0;

      /* Selector: '<S8>/Selector1' */
      ACDS443_Prop_B.Selector1[3 * u0_0] = ACDS443_Prop_B.MatrixConcatenation[6 *
        u0_0 + 3];
      ACDS443_Prop_B.Product_h[u0_0] += y[u0_0 + 3] * u_idx_1;

      /* Selector: '<S8>/Selector1' */
      ACDS443_Prop_B.Selector1[1 + 3 * u0_0] =
        ACDS443_Prop_B.MatrixConcatenation[6 * u0_0 + 4];
      ACDS443_Prop_B.Product_h[u0_0] += y[u0_0 + 6] * u0;

      /* Selector: '<S8>/Selector1' */
      ACDS443_Prop_B.Selector1[2 + 3 * u0_0] =
        ACDS443_Prop_B.MatrixConcatenation[6 * u0_0 + 5];
    }

    /* End of Product: '<S19>/Product' */

    /* Product: '<S21>/i x j' */
    ACDS443_Prop_B.ixj = ACDS443_Prop_B.pqr[0] * ACDS443_Prop_B.Product_h[1];

    /* Product: '<S21>/j x k' */
    ACDS443_Prop_B.jxk = ACDS443_Prop_B.pqr[1] * ACDS443_Prop_B.Product_h[2];

    /* Product: '<S21>/k x i' */
    ACDS443_Prop_B.kxi = ACDS443_Prop_B.pqr[2] * ACDS443_Prop_B.Product_h[0];

    /* Product: '<S22>/i x k' */
    ACDS443_Prop_B.ixk = ACDS443_Prop_B.pqr[0] * ACDS443_Prop_B.Product_h[2];

    /* Product: '<S22>/j x i' */
    ACDS443_Prop_B.jxi = ACDS443_Prop_B.pqr[1] * ACDS443_Prop_B.Product_h[0];

    /* Product: '<S22>/k x j' */
    ACDS443_Prop_B.kxj = ACDS443_Prop_B.pqr[2] * ACDS443_Prop_B.Product_h[1];

    /* Sum: '<S18>/Sum' */
    ACDS443_Prop_B.Sum[0] = ACDS443_Prop_B.jxk - ACDS443_Prop_B.kxj;
    ACDS443_Prop_B.Sum[1] = ACDS443_Prop_B.kxi - ACDS443_Prop_B.ixk;
    ACDS443_Prop_B.Sum[2] = ACDS443_Prop_B.ixj - ACDS443_Prop_B.jxi;

    /* Product: '<S20>/Product' */
    memcpy(&y[0], &ACDS443_Prop_B.Selector1[0], 9U * sizeof(real_T));
    u_idx_0 = ACDS443_Prop_B.pqr[0];
    u_idx_1 = ACDS443_Prop_B.pqr[1];
    u0 = ACDS443_Prop_B.pqr[2];
    for (u0_0 = 0; u0_0 < 3; u0_0++) {
      ACDS443_Prop_B.Product_g[u0_0] = 0.0;
      ACDS443_Prop_B.Product_g[u0_0] += y[u0_0] * u_idx_0;
      ACDS443_Prop_B.Product_g[u0_0] += y[u0_0 + 3] * u_idx_1;
      ACDS443_Prop_B.Product_g[u0_0] += y[u0_0 + 6] * u0;
    }

    /* End of Product: '<S20>/Product' */
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
      /* Sum: '<Root>/Sum' */
      ACDS443_Prop_B.Sum_e[0] = ACDS443_Prop_B.AerodynamicDrag[0] +
        ACDS443_Prop_B.SRPandGravity[0];
      ACDS443_Prop_B.Sum_e[1] = ACDS443_Prop_B.AerodynamicDrag[1] +
        ACDS443_Prop_B.SRPandGravity[1];
      ACDS443_Prop_B.Sum_e[2] = ACDS443_Prop_B.AerodynamicDrag[2] +
        ACDS443_Prop_B.SRPandGravity[2];
    }

    for (i = 0; i < 3; i++) {
      /* Sum: '<Root>/Subtract' */
      ACDS443_Prop_B.Subtract[i] = ACDS443_Prop_B.Sum_e[i] -
        ACDS443_Prop_B.Add[i];

      /* Sum: '<S8>/Sum2' */
      ACDS443_Prop_B.Sum2_k[i] = (ACDS443_Prop_B.Subtract[i] -
        ACDS443_Prop_B.Product_g[i]) - ACDS443_Prop_B.Sum[i];

      /* Selector: '<S8>/Selector2' */
      ACDS443_Prop_B.Selector2[3 * i] = ACDS443_Prop_B.MatrixConcatenation[6 * i];
      ACDS443_Prop_B.Selector2[1 + 3 * i] = ACDS443_Prop_B.MatrixConcatenation[6
        * i + 1];
      ACDS443_Prop_B.Selector2[2 + 3 * i] = ACDS443_Prop_B.MatrixConcatenation[6
        * i + 2];
    }

    /* Product: '<S8>/Product2' */
    rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(ACDS443_Prop_B.Sum2_k,
      ACDS443_Prop_B.Selector2, ACDS443_Prop_B.Product2);

    /* Switch: '<S9>/Switch' */
    if (ACDS443_Prop_B.mass_o2 != 0.0) {
      ACDS443_Prop_B.Switch = 0.0;
    } else {
      /* Product: '<S9>/Product' */
      ACDS443_Prop_B.Product_b = 0.0 * ACDS443_Prop_B.Sum3_p;
      ACDS443_Prop_B.Switch = ACDS443_Prop_B.Product_b;
    }

    /* End of Switch: '<S9>/Switch' */

    /* Sum: '<S9>/Sum' incorporates:
     *  Constant: '<Root>/Constant'
     */
    ACDS443_Prop_B.Sum_c[0] = ACDS443_Prop_P.Constant_Value[0] +
      ACDS443_Prop_B.Switch;

    /* Product: '<S3>/Product' */
    ACDS443_Prop_B.Product_d[0] = ACDS443_Prop_B.Sum_c[0] /
      ACDS443_Prop_B.mass_o1;

    /* Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_B.ubvbwb[0] = ACDS443_Prop_X.U[0];

    /* Sum: '<S9>/Sum' incorporates:
     *  Constant: '<Root>/Constant'
     */
    ACDS443_Prop_B.Sum_c[1] = ACDS443_Prop_P.Constant_Value[1] +
      ACDS443_Prop_B.Switch;

    /* Product: '<S3>/Product' */
    ACDS443_Prop_B.Product_d[1] = ACDS443_Prop_B.Sum_c[1] /
      ACDS443_Prop_B.mass_o1;

    /* Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_B.ubvbwb[1] = ACDS443_Prop_X.U[1];

    /* Sum: '<S9>/Sum' incorporates:
     *  Constant: '<Root>/Constant'
     */
    ACDS443_Prop_B.Sum_c[2] = ACDS443_Prop_P.Constant_Value[2] +
      ACDS443_Prop_B.Switch;

    /* Product: '<S3>/Product' */
    ACDS443_Prop_B.Product_d[2] = ACDS443_Prop_B.Sum_c[2] /
      ACDS443_Prop_B.mass_o1;

    /* Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_B.ubvbwb[2] = ACDS443_Prop_X.U[2];

    /* Product: '<S26>/j x k' */
    ACDS443_Prop_B.jxk_e = ACDS443_Prop_B.ubvbwb[1] * ACDS443_Prop_B.pqr[2];

    /* Product: '<S26>/k x i' */
    ACDS443_Prop_B.kxi_k = ACDS443_Prop_B.ubvbwb[2] * ACDS443_Prop_B.pqr[0];

    /* Product: '<S26>/i x j' */
    ACDS443_Prop_B.ixj_f = ACDS443_Prop_B.ubvbwb[0] * ACDS443_Prop_B.pqr[1];

    /* Product: '<S27>/k x j' */
    ACDS443_Prop_B.kxj_e = ACDS443_Prop_B.ubvbwb[2] * ACDS443_Prop_B.pqr[1];

    /* Product: '<S27>/i x k' */
    ACDS443_Prop_B.ixk_a = ACDS443_Prop_B.ubvbwb[0] * ACDS443_Prop_B.pqr[2];

    /* Product: '<S27>/j x i' */
    ACDS443_Prop_B.jxi_m = ACDS443_Prop_B.ubvbwb[1] * ACDS443_Prop_B.pqr[0];

    /* Sum: '<S10>/Sum' */
    ACDS443_Prop_B.Sum_d[0] = ACDS443_Prop_B.jxk_e - ACDS443_Prop_B.kxj_e;
    ACDS443_Prop_B.Sum_d[1] = ACDS443_Prop_B.kxi_k - ACDS443_Prop_B.ixk_a;
    ACDS443_Prop_B.Sum_d[2] = ACDS443_Prop_B.ixj_f - ACDS443_Prop_B.jxi_m;
    for (i = 0; i < 3; i++) {
      /* Sum: '<S3>/Sum' */
      ACDS443_Prop_B.Sum_i[i] = ACDS443_Prop_B.Product_d[i] +
        ACDS443_Prop_B.Sum_d[i];

      /* Math: '<S3>/Transpose' */
      ACDS443_Prop_B.Transpose[3 * i] = ACDS443_Prop_B.VectorConcatenate[i];
      ACDS443_Prop_B.Transpose[1 + 3 * i] = ACDS443_Prop_B.VectorConcatenate[i +
        3];
      ACDS443_Prop_B.Transpose[2 + 3 * i] = ACDS443_Prop_B.VectorConcatenate[i +
        6];
    }

    /* Product: '<S14>/Product' */
    memcpy(&y[0], &ACDS443_Prop_B.Transpose[0], 9U * sizeof(real_T));
    u_idx_0 = ACDS443_Prop_B.ubvbwb[0];
    u_idx_1 = ACDS443_Prop_B.ubvbwb[1];
    u0 = ACDS443_Prop_B.ubvbwb[2];
    for (i = 0; i < 3; i++) {
      ACDS443_Prop_B.Product_o[i] = 0.0;
      ACDS443_Prop_B.Product_o[i] += y[i] * u_idx_0;
      ACDS443_Prop_B.Product_o[i] += y[i + 3] * u_idx_1;
      ACDS443_Prop_B.Product_o[i] += y[i + 6] * u0;

      /* Integrator: '<S3>/xe,ye,ze' */
      ACDS443_Prop_B.xeyeze[i] = ACDS443_Prop_X.Xe[i];

      /* Sum: '<Root>/Sum1' */
      ACDS443_Prop_B.Sum1[i] = (ACDS443_Prop_B.Gain1[i] + ACDS443_Prop_B.Gain[i])
        - ACDS443_Prop_B.Add1[i];
    }

    /* End of Product: '<S14>/Product' */
  }

  if (rtmIsMajorTimeStep(ACDS443_Prop_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(ACDS443_Prop_M->rtwLogInfo, (ACDS443_Prop_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(ACDS443_Prop_M)) {
    real_T (*lastU)[3];
    real_T u1;
    real_T u2;

    /* Update for Derivative: '<S6>/Derivative' */
    if (ACDS443_Prop_DW.TimeStampA == (rtInf)) {
      ACDS443_Prop_DW.TimeStampA = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA;
    } else if (ACDS443_Prop_DW.TimeStampB == (rtInf)) {
      ACDS443_Prop_DW.TimeStampB = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB;
    } else if (ACDS443_Prop_DW.TimeStampA < ACDS443_Prop_DW.TimeStampB) {
      ACDS443_Prop_DW.TimeStampA = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA;
    } else {
      ACDS443_Prop_DW.TimeStampB = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB;
    }

    (*lastU)[0] = ACDS443_Prop_B.MomentumSaturation[0];
    (*lastU)[1] = ACDS443_Prop_B.MomentumSaturation[1];
    (*lastU)[2] = ACDS443_Prop_B.MomentumSaturation[2];

    /* End of Update for Derivative: '<S6>/Derivative' */

    /* Update for Derivative: '<S5>/Derivative' */
    if (ACDS443_Prop_DW.TimeStampA_i == (rtInf)) {
      ACDS443_Prop_DW.TimeStampA_i = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_f;
    } else if (ACDS443_Prop_DW.TimeStampB_j == (rtInf)) {
      ACDS443_Prop_DW.TimeStampB_j = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB_a;
    } else if (ACDS443_Prop_DW.TimeStampA_i < ACDS443_Prop_DW.TimeStampB_j) {
      ACDS443_Prop_DW.TimeStampA_i = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_f;
    } else {
      ACDS443_Prop_DW.TimeStampB_j = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB_a;
    }

    (*lastU)[0] = ACDS443_Prop_B.MomentumSaturation_c[0];
    (*lastU)[1] = ACDS443_Prop_B.MomentumSaturation_c[1];
    (*lastU)[2] = ACDS443_Prop_B.MomentumSaturation_c[2];

    /* End of Update for Derivative: '<S5>/Derivative' */

    /* Update for Derivative: '<S4>/Derivative' */
    if (ACDS443_Prop_DW.TimeStampA_j == (rtInf)) {
      ACDS443_Prop_DW.TimeStampA_j = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_b;
    } else if (ACDS443_Prop_DW.TimeStampB_m == (rtInf)) {
      ACDS443_Prop_DW.TimeStampB_m = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB_c;
    } else if (ACDS443_Prop_DW.TimeStampA_j < ACDS443_Prop_DW.TimeStampB_m) {
      ACDS443_Prop_DW.TimeStampA_j = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeA_b;
    } else {
      ACDS443_Prop_DW.TimeStampB_m = ACDS443_Prop_M->Timing.t[0];
      lastU = &ACDS443_Prop_DW.LastUAtTimeB_c;
    }

    (*lastU)[0] = ACDS443_Prop_B.MomentumSaturation_h[0];
    (*lastU)[1] = ACDS443_Prop_B.MomentumSaturation_h[1];
    (*lastU)[2] = ACDS443_Prop_B.MomentumSaturation_h[2];

    /* End of Update for Derivative: '<S4>/Derivative' */
    if (rtmIsMajorTimeStep(ACDS443_Prop_M) &&
        ACDS443_Prop_M->Timing.TaskCounters.TID[1] == 0) {
      /* Update for RandomNumber: '<Root>/Aerodynamic Drag' */
      u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[0];
      u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[0];
      ACDS443_Prop_DW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed[0]) * u1 + u2;
      u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[1];
      u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[1];
      ACDS443_Prop_DW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed[1]) * u1 + u2;
      u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[2];
      u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[2];
      ACDS443_Prop_DW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed[2]) * u1 + u2;

      /* Update for RandomNumber: '<Root>/SRP and Gravity ' */
      u1 = ACDS443_Prop_P.SRPandGravity_StdDev[0];
      u2 = ACDS443_Prop_P.SRPandGravity_Mean[0];
      ACDS443_Prop_DW.NextOutput_k[0] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed_c[0]) * u1 + u2;
      u1 = ACDS443_Prop_P.SRPandGravity_StdDev[1];
      u2 = ACDS443_Prop_P.SRPandGravity_Mean[1];
      ACDS443_Prop_DW.NextOutput_k[1] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed_c[1]) * u1 + u2;
      u1 = ACDS443_Prop_P.SRPandGravity_StdDev[2];
      u2 = ACDS443_Prop_P.SRPandGravity_Mean[2];
      ACDS443_Prop_DW.NextOutput_k[2] = rt_nrand_Upu32_Yd_f_pw_snf
        (&ACDS443_Prop_DW.RandSeed_c[2]) * u1 + u2;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(ACDS443_Prop_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(ACDS443_Prop_M)!=-1) &&
          !((rtmGetTFinal(ACDS443_Prop_M)-(((ACDS443_Prop_M->Timing.clockTick1+
               ACDS443_Prop_M->Timing.clockTickH1* 4294967296.0)) * 0.1)) >
            (((ACDS443_Prop_M->Timing.clockTick1+
               ACDS443_Prop_M->Timing.clockTickH1* 4294967296.0)) * 0.1) *
            (DBL_EPSILON))) {
        rtmSetErrorStatus(ACDS443_Prop_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&ACDS443_Prop_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++ACDS443_Prop_M->Timing.clockTick0)) {
      ++ACDS443_Prop_M->Timing.clockTickH0;
    }

    ACDS443_Prop_M->Timing.t[0] = rtsiGetSolverStopTime
      (&ACDS443_Prop_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.1s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.1, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      ACDS443_Prop_M->Timing.clockTick1++;
      if (!ACDS443_Prop_M->Timing.clockTick1) {
        ACDS443_Prop_M->Timing.clockTickH1++;
      }
    }

    rate_scheduler();
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void ACDS443_Prop_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_ACDS443_Prop_T *_rtXdot;
  _rtXdot = ((XDot_ACDS443_Prop_T *) ACDS443_Prop_M->derivs);

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE = ACDS443_Prop_B.DotProduct;

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_i = ACDS443_Prop_B.DotProduct_j;

  /* Derivatives for Integrator: '<S6>/Integrator' */
  _rtXdot->Integrator_CSTATE_k = ACDS443_Prop_B.DotProduct_f;

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[0] = ACDS443_Prop_B.TmpSignalConversionAtphithetaps[0];

  /* Derivatives for Integrator: '<S3>/p,q,r ' */
  _rtXdot->p[0] = ACDS443_Prop_B.Product2[0];

  /* Derivatives for Integrator: '<S6>/Integrator1' */
  _rtXdot->Integrator1_CSTATE[0] = ACDS443_Prop_B.Sum1[0];

  /* Derivatives for Integrator: '<S5>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_f[0] = ACDS443_Prop_B.Sum1[0];

  /* Derivatives for Integrator: '<S4>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_a[0] = ACDS443_Prop_B.Sum1[0];

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[0] = ACDS443_Prop_B.Sum2[0];

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[1] = ACDS443_Prop_B.TmpSignalConversionAtphithetaps[1];

  /* Derivatives for Integrator: '<S3>/p,q,r ' */
  _rtXdot->p[1] = ACDS443_Prop_B.Product2[1];

  /* Derivatives for Integrator: '<S6>/Integrator1' */
  _rtXdot->Integrator1_CSTATE[1] = ACDS443_Prop_B.Sum1[1];

  /* Derivatives for Integrator: '<S5>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_f[1] = ACDS443_Prop_B.Sum1[1];

  /* Derivatives for Integrator: '<S4>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_a[1] = ACDS443_Prop_B.Sum1[1];

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[1] = ACDS443_Prop_B.Sum2[1];

  /* Derivatives for Integrator: '<S7>/phi theta psi' */
  _rtXdot->phi[2] = ACDS443_Prop_B.TmpSignalConversionAtphithetaps[2];

  /* Derivatives for Integrator: '<S3>/p,q,r ' */
  _rtXdot->p[2] = ACDS443_Prop_B.Product2[2];

  /* Derivatives for Integrator: '<S6>/Integrator1' */
  _rtXdot->Integrator1_CSTATE[2] = ACDS443_Prop_B.Sum1[2];

  /* Derivatives for Integrator: '<S5>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_f[2] = ACDS443_Prop_B.Sum1[2];

  /* Derivatives for Integrator: '<S4>/Integrator1' */
  _rtXdot->Integrator1_CSTATE_a[2] = ACDS443_Prop_B.Sum1[2];

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[2] = ACDS443_Prop_B.Sum2[2];

  /* Derivatives for Integrator: '<S9>/mass  ' */
  lsat = (ACDS443_Prop_X.mass <= ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_l);
  usat = (ACDS443_Prop_X.mass >= ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_i);
  if (((!lsat) && (!usat)) || (lsat && (ACDS443_Prop_B.Sum3_p > 0.0)) || (usat &&
       (ACDS443_Prop_B.Sum3_p < 0.0))) {
    _rtXdot->mass = ACDS443_Prop_B.Sum3_p;
  } else {
    /* in saturation */
    _rtXdot->mass = 0.0;
  }

  /* End of Derivatives for Integrator: '<S9>/mass  ' */

  /* Derivatives for Integrator: '<S3>/ub,vb,wb' */
  _rtXdot->U[0] = ACDS443_Prop_B.Sum_i[0];

  /* Derivatives for Integrator: '<S3>/xe,ye,ze' */
  _rtXdot->Xe[0] = ACDS443_Prop_B.Product_o[0];

  /* Derivatives for Integrator: '<S3>/ub,vb,wb' */
  _rtXdot->U[1] = ACDS443_Prop_B.Sum_i[1];

  /* Derivatives for Integrator: '<S3>/xe,ye,ze' */
  _rtXdot->Xe[1] = ACDS443_Prop_B.Product_o[1];

  /* Derivatives for Integrator: '<S3>/ub,vb,wb' */
  _rtXdot->U[2] = ACDS443_Prop_B.Sum_i[2];

  /* Derivatives for Integrator: '<S3>/xe,ye,ze' */
  _rtXdot->Xe[2] = ACDS443_Prop_B.Product_o[2];
}

/* Model initialize function */
void ACDS443_Prop_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)ACDS443_Prop_M, 0,
                sizeof(RT_MODEL_ACDS443_Prop_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&ACDS443_Prop_M->solverInfo,
                          &ACDS443_Prop_M->Timing.simTimeStep);
    rtsiSetTPtr(&ACDS443_Prop_M->solverInfo, &rtmGetTPtr(ACDS443_Prop_M));
    rtsiSetStepSizePtr(&ACDS443_Prop_M->solverInfo,
                       &ACDS443_Prop_M->Timing.stepSize0);
    rtsiSetdXPtr(&ACDS443_Prop_M->solverInfo, &ACDS443_Prop_M->derivs);
    rtsiSetContStatesPtr(&ACDS443_Prop_M->solverInfo, (real_T **)
                         &ACDS443_Prop_M->contStates);
    rtsiSetNumContStatesPtr(&ACDS443_Prop_M->solverInfo,
      &ACDS443_Prop_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&ACDS443_Prop_M->solverInfo,
      &ACDS443_Prop_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&ACDS443_Prop_M->solverInfo,
      &ACDS443_Prop_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&ACDS443_Prop_M->solverInfo,
      &ACDS443_Prop_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&ACDS443_Prop_M->solverInfo, (&rtmGetErrorStatus
      (ACDS443_Prop_M)));
    rtsiSetRTModelPtr(&ACDS443_Prop_M->solverInfo, ACDS443_Prop_M);
  }

  rtsiSetSimTimeStep(&ACDS443_Prop_M->solverInfo, MAJOR_TIME_STEP);
  ACDS443_Prop_M->intgData.y = ACDS443_Prop_M->odeY;
  ACDS443_Prop_M->intgData.f[0] = ACDS443_Prop_M->odeF[0];
  ACDS443_Prop_M->intgData.f[1] = ACDS443_Prop_M->odeF[1];
  ACDS443_Prop_M->intgData.f[2] = ACDS443_Prop_M->odeF[2];
  ACDS443_Prop_M->contStates = ((X_ACDS443_Prop_T *) &ACDS443_Prop_X);
  ACDS443_Prop_M->periodicContStateIndices = ((int_T*) ACDS443_Prop_PeriodicIndX);
  ACDS443_Prop_M->periodicContStateRanges = ((real_T*) ACDS443_Prop_PeriodicRngX);
  rtsiSetSolverData(&ACDS443_Prop_M->solverInfo, (void *)
                    &ACDS443_Prop_M->intgData);
  rtsiSetSolverName(&ACDS443_Prop_M->solverInfo,"ode3");
  rtmSetTPtr(ACDS443_Prop_M, &ACDS443_Prop_M->Timing.tArray[0]);
  rtmSetTFinal(ACDS443_Prop_M, 6000.0);
  ACDS443_Prop_M->Timing.stepSize0 = 0.1;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    ACDS443_Prop_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(ACDS443_Prop_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(ACDS443_Prop_M->rtwLogInfo, (NULL));
    rtliSetLogT(ACDS443_Prop_M->rtwLogInfo, "tout");
    rtliSetLogX(ACDS443_Prop_M->rtwLogInfo, "");
    rtliSetLogXFinal(ACDS443_Prop_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(ACDS443_Prop_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(ACDS443_Prop_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(ACDS443_Prop_M->rtwLogInfo, 0);
    rtliSetLogDecimation(ACDS443_Prop_M->rtwLogInfo, 1);
    rtliSetLogY(ACDS443_Prop_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(ACDS443_Prop_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(ACDS443_Prop_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &ACDS443_Prop_B), 0,
                sizeof(B_ACDS443_Prop_T));

  /* states (continuous) */
  {
    (void) memset((void *)&ACDS443_Prop_X, 0,
                  sizeof(X_ACDS443_Prop_T));
  }

  /* Periodic continuous states */
  {
    (void) memset((void*) ACDS443_Prop_PeriodicIndX, 0,
                  3*sizeof(int_T));
    (void) memset((void*) ACDS443_Prop_PeriodicRngX, 0,
                  6*sizeof(real_T));
  }

  /* states (dwork) */
  (void) memset((void *)&ACDS443_Prop_DW, 0,
                sizeof(DW_ACDS443_Prop_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(ACDS443_Prop_M->rtwLogInfo, 0.0, rtmGetTFinal
    (ACDS443_Prop_M), ACDS443_Prop_M->Timing.stepSize0, (&rtmGetErrorStatus
    (ACDS443_Prop_M)));

  {
    uint32_T tseed;
    int32_T r;
    int32_T t;
    real_T u0;
    real_T u1;
    real_T u2;

    /* InitializeConditions for Integrator: '<S4>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE = ACDS443_Prop_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE_i = ACDS443_Prop_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S6>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE_k = ACDS443_Prop_P.Integrator_IC_l;

    /* InitializeConditions for Derivative: '<S6>/Derivative' */
    ACDS443_Prop_DW.TimeStampA = (rtInf);
    ACDS443_Prop_DW.TimeStampB = (rtInf);

    /* InitializeConditions for Derivative: '<S5>/Derivative' */
    ACDS443_Prop_DW.TimeStampA_i = (rtInf);
    ACDS443_Prop_DW.TimeStampB_j = (rtInf);

    /* InitializeConditions for Derivative: '<S4>/Derivative' */
    ACDS443_Prop_DW.TimeStampA_j = (rtInf);
    ACDS443_Prop_DW.TimeStampB_m = (rtInf);

    /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
    ACDS443_Prop_X.phi[0] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_k[0];

    /* InitializeConditions for Integrator: '<S3>/p,q,r ' */
    ACDS443_Prop_X.p[0] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_a[0];

    /* InitializeConditions for Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE[0] = ACDS443_Prop_P.Integrator1_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_f[0] = ACDS443_Prop_P.Integrator1_IC_g;

    /* InitializeConditions for Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_a[0] = ACDS443_Prop_P.Integrator1_IC_i;

    /* InitializeConditions for RandomNumber: '<Root>/Aerodynamic Drag' */
    u0 = ACDS443_Prop_P.AerodynamicDrag_Seed[0];
    u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[0];
    u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[0];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput[0] = u1;
    ACDS443_Prop_DW.RandSeed[0] = tseed;

    /* InitializeConditions for Integrator: '<Root>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE_j[0] = ACDS443_Prop_P.Integrator_IC_d;

    /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
    ACDS443_Prop_X.phi[1] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_k[1];

    /* InitializeConditions for Integrator: '<S3>/p,q,r ' */
    ACDS443_Prop_X.p[1] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_a[1];

    /* InitializeConditions for Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE[1] = ACDS443_Prop_P.Integrator1_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_f[1] = ACDS443_Prop_P.Integrator1_IC_g;

    /* InitializeConditions for Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_a[1] = ACDS443_Prop_P.Integrator1_IC_i;

    /* InitializeConditions for RandomNumber: '<Root>/Aerodynamic Drag' */
    u0 = ACDS443_Prop_P.AerodynamicDrag_Seed[1];
    u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[1];
    u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[1];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput[1] = u1;
    ACDS443_Prop_DW.RandSeed[1] = tseed;

    /* InitializeConditions for Integrator: '<Root>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE_j[1] = ACDS443_Prop_P.Integrator_IC_d;

    /* InitializeConditions for Integrator: '<S7>/phi theta psi' */
    ACDS443_Prop_X.phi[2] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_k[2];

    /* InitializeConditions for Integrator: '<S3>/p,q,r ' */
    ACDS443_Prop_X.p[2] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_a[2];

    /* InitializeConditions for Integrator: '<S6>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE[2] = ACDS443_Prop_P.Integrator1_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_f[2] = ACDS443_Prop_P.Integrator1_IC_g;

    /* InitializeConditions for Integrator: '<S4>/Integrator1' */
    ACDS443_Prop_X.Integrator1_CSTATE_a[2] = ACDS443_Prop_P.Integrator1_IC_i;

    /* InitializeConditions for RandomNumber: '<Root>/Aerodynamic Drag' */
    u0 = ACDS443_Prop_P.AerodynamicDrag_Seed[2];
    u1 = ACDS443_Prop_P.AerodynamicDrag_StdDev[2];
    u2 = ACDS443_Prop_P.AerodynamicDrag_Mean[2];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput[2] = u1;
    ACDS443_Prop_DW.RandSeed[2] = tseed;

    /* InitializeConditions for Integrator: '<Root>/Integrator' */
    ACDS443_Prop_X.Integrator_CSTATE_j[2] = ACDS443_Prop_P.Integrator_IC_d;

    /* InitializeConditions for Integrator: '<S9>/mass  ' */
    ACDS443_Prop_X.mass = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_g;

    /* InitializeConditions for RandomNumber: '<Root>/SRP and Gravity ' */
    u0 = ACDS443_Prop_P.SRPandGravity_Seed[0];
    u1 = ACDS443_Prop_P.SRPandGravity_StdDev[0];
    u2 = ACDS443_Prop_P.SRPandGravity_Mean[0];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput_k[0] = u1;
    ACDS443_Prop_DW.RandSeed_c[0] = tseed;

    /* InitializeConditions for Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_X.U[0] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAngl[0];

    /* InitializeConditions for Integrator: '<S3>/xe,ye,ze' */
    ACDS443_Prop_X.Xe[0] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_e[0];

    /* InitializeConditions for RandomNumber: '<Root>/SRP and Gravity ' */
    u0 = ACDS443_Prop_P.SRPandGravity_Seed[1];
    u1 = ACDS443_Prop_P.SRPandGravity_StdDev[1];
    u2 = ACDS443_Prop_P.SRPandGravity_Mean[1];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput_k[1] = u1;
    ACDS443_Prop_DW.RandSeed_c[1] = tseed;

    /* InitializeConditions for Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_X.U[1] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAngl[1];

    /* InitializeConditions for Integrator: '<S3>/xe,ye,ze' */
    ACDS443_Prop_X.Xe[1] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_e[1];

    /* InitializeConditions for RandomNumber: '<Root>/SRP and Gravity ' */
    u0 = ACDS443_Prop_P.SRPandGravity_Seed[2];
    u1 = ACDS443_Prop_P.SRPandGravity_StdDev[2];
    u2 = ACDS443_Prop_P.SRPandGravity_Mean[2];
    u0 = floor(u0);
    if (rtIsNaN(u0) || rtIsInf(u0)) {
      u0 = 0.0;
    } else {
      u0 = fmod(u0, 4.294967296E+9);
    }

    tseed = u0 < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-u0 : (uint32_T)u0;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    u1 = rt_nrand_Upu32_Yd_f_pw_snf(&tseed) * u1 + u2;
    ACDS443_Prop_DW.NextOutput_k[2] = u1;
    ACDS443_Prop_DW.RandSeed_c[2] = tseed;

    /* InitializeConditions for Integrator: '<S3>/ub,vb,wb' */
    ACDS443_Prop_X.U[2] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAngl[2];

    /* InitializeConditions for Integrator: '<S3>/xe,ye,ze' */
    ACDS443_Prop_X.Xe[2] = ACDS443_Prop_P.SimpleVariableMass6DOFEulerAn_e[2];

    /* InitializeConditions for root-level periodic continuous states */
    {
      int_T rootPeriodicContStateIndices[3] = { 0, 1, 2 };

      real_T rootPeriodicContStateRanges[6] = { -3.1415926535897931,
        3.1415926535897931, -3.1415926535897931, 3.1415926535897931,
        -3.1415926535897931, 3.1415926535897931 };

      (void) memcpy((void*)ACDS443_Prop_PeriodicIndX,
                    rootPeriodicContStateIndices,
                    3*sizeof(int_T));
      (void) memcpy((void*)ACDS443_Prop_PeriodicRngX,
                    rootPeriodicContStateRanges,
                    6*sizeof(real_T));
    }
  }
}

/* Model terminate function */
void ACDS443_Prop_terminate(void)
{
  /* (no terminate code required) */
}
