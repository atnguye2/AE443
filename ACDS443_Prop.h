/*
 * ACDS443_Prop.h
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

#ifndef RTW_HEADER_ACDS443_Prop_h_
#define RTW_HEADER_ACDS443_Prop_h_
#include <string.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#ifndef ACDS443_Prop_COMMON_INCLUDES_
# define ACDS443_Prop_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#endif                                 /* ACDS443_Prop_COMMON_INCLUDES_ */

#include "ACDS443_Prop_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T phithetapsi[3];               /* '<S7>/phi theta psi' */
  real_T pqr[3];                       /* '<S3>/p,q,r ' */
  real_T Integrator;                   /* '<S4>/Integrator' */
  real_T Integrator_a;                 /* '<S5>/Integrator' */
  real_T Integrator_a1;                /* '<S6>/Integrator' */
  real_T Integrator1[3];               /* '<S6>/Integrator1' */
  real_T MomentumSaturation[3];        /* '<S6>/Momentum Saturation' */
  real_T Derivative[3];                /* '<S6>/Derivative' */
  real_T MomentInputSaturation[3];     /* '<S6>/Moment Input Saturation' */
  real_T Multiply[3];                  /* '<S6>/Multiply' */
  real_T Integrator1_b[3];             /* '<S5>/Integrator1' */
  real_T MomentumSaturation_c[3];      /* '<S5>/Momentum Saturation' */
  real_T Derivative_p[3];              /* '<S5>/Derivative' */
  real_T MomentInputSaturation_o[3];   /* '<S5>/Moment Input Saturation' */
  real_T Multiply_b[3];                /* '<S5>/Multiply' */
  real_T Integrator1_f[3];             /* '<S4>/Integrator1' */
  real_T MomentumSaturation_h[3];      /* '<S4>/Momentum Saturation' */
  real_T Derivative_m[3];              /* '<S4>/Derivative' */
  real_T MomentInputSaturation_g[3];   /* '<S4>/Moment Input Saturation' */
  real_T Multiply_b4[3];               /* '<S4>/Multiply' */
  real_T Multiply1[3];                 /* '<S1>/Multiply1' */
  real_T Multiply_n[3];                /* '<S1>/Multiply' */
  real_T Add1[3];                      /* '<Root>/Add1' */
  real_T Add[3];                       /* '<Root>/Add' */
  real_T AerodynamicDrag[3];           /* '<Root>/Aerodynamic Drag' */
  real_T Gain[3];                      /* '<Root>/Gain' */
  real_T Gain1[3];                     /* '<Root>/Gain1' */
  real_T Integrator_g[3];              /* '<Root>/Integrator' */
  real_T Gain2[3];                     /* '<Root>/Gain2' */
  real_T DataTypeConversion[3];        /* '<S1>/Data Type Conversion' */
  real_T DataTypeConversion1[3];       /* '<S1>/Data Type Conversion1' */
  real_T Sum2[3];                      /* '<S1>/Sum2' */
  real_T DotProduct;                   /* '<S4>/Dot Product' */
  real_T DotProduct_j;                 /* '<S5>/Dot Product' */
  real_T DotProduct_f;                 /* '<S6>/Dot Product' */
  real_T SRPandGravity[3];             /* '<Root>/SRP and Gravity ' */
  real_T TmpSignalConversionAtsincosInpo[3];
  real_T sincos_o1[3];                 /* '<S15>/sincos' */
  real_T sincos_o2[3];                 /* '<S15>/sincos' */
  real_T VectorConcatenate[9];         /* '<S17>/Vector Concatenate' */
  real_T sincos_o1_m[3];               /* '<S16>/sincos' */
  real_T sincos_o2_i[3];               /* '<S16>/sincos' */
  real_T phidot;                       /* '<S16>/phidot' */
  real_T thetadot;                     /* '<S16>/thetadot' */
  real_T psidot;                       /* '<S16>/psidot' */
  real_T TmpSignalConversionAtphithetaps[3];/* '<S7>/phidot thetadot psidot' */
  real_T mass_o1;                      /* '<S9>/mass  ' */
  real_T mass_o2;                      /* '<S9>/mass  ' */
  real_T TmpSignalConversionAtPreLookUpI[2];
  real_T PreLookUpIndexSearch_o2;      /* '<S23>/PreLook-Up Index Search' */
  real_T u[9];                         /* '<S25>/[0]' */
  real_T Fcn;                          /* '<S25>/Fcn' */
  real_T ulambda_x[9];                 /* '<S25>/1-lambda_x' */
  real_T u_o[9];                       /* '<S25>/[1]' */
  real_T lambda_x[9];                  /* '<S25>/lambda_x.' */
  real_T Sum3[9];                      /* '<S25>/Sum3' */
  real_T Sum3_p;                       /* '<Root>/Sum3' */
  real_T Product[9];                   /* '<S23>/Product' */
  real_T MatrixConcatenation[18];      /* '<S9>/Matrix Concatenation' */
  real_T Selector[9];                  /* '<S8>/Selector' */
  real_T Product_h[3];                 /* '<S19>/Product' */
  real_T ixj;                          /* '<S21>/i x j' */
  real_T jxk;                          /* '<S21>/j x k' */
  real_T kxi;                          /* '<S21>/k x i' */
  real_T ixk;                          /* '<S22>/i x k' */
  real_T jxi;                          /* '<S22>/j x i' */
  real_T kxj;                          /* '<S22>/k x j' */
  real_T Sum[3];                       /* '<S18>/Sum' */
  real_T Selector1[9];                 /* '<S8>/Selector1' */
  real_T Product_g[3];                 /* '<S20>/Product' */
  real_T Sum_e[3];                     /* '<Root>/Sum' */
  real_T Subtract[3];                  /* '<Root>/Subtract' */
  real_T Sum2_k[3];                    /* '<S8>/Sum2' */
  real_T Selector2[9];                 /* '<S8>/Selector2' */
  real_T Product2[3];                  /* '<S8>/Product2' */
  real_T Switch;                       /* '<S9>/Switch' */
  real_T Sum_c[3];                     /* '<S9>/Sum' */
  real_T Product_d[3];                 /* '<S3>/Product' */
  real_T ubvbwb[3];                    /* '<S3>/ub,vb,wb' */
  real_T jxk_e;                        /* '<S26>/j x k' */
  real_T kxi_k;                        /* '<S26>/k x i' */
  real_T ixj_f;                        /* '<S26>/i x j' */
  real_T kxj_e;                        /* '<S27>/k x j' */
  real_T ixk_a;                        /* '<S27>/i x k' */
  real_T jxi_m;                        /* '<S27>/j x i' */
  real_T Sum_d[3];                     /* '<S10>/Sum' */
  real_T Sum_i[3];                     /* '<S3>/Sum' */
  real_T Transpose[9];                 /* '<S3>/Transpose' */
  real_T Product_o[3];                 /* '<S14>/Product' */
  real_T xeyeze[3];                    /* '<S3>/xe,ye,ze' */
  real_T Sum1[3];                      /* '<Root>/Sum1' */
  real_T Product_b;                    /* '<S9>/Product' */
  int32_T PreLookUpIndexSearch_o1;     /* '<S23>/PreLook-Up Index Search' */
  int32_T Sum_dq;                      /* '<S25>/Sum' */
  boolean_T LessThan[3];               /* '<S1>/Less Than' */
  boolean_T ZeroOrderHold1[3];         /* '<S1>/Zero-Order Hold1' */
  boolean_T GreaterThan[3];            /* '<S1>/GreaterThan' */
  boolean_T ZeroOrderHold[3];          /* '<S1>/Zero-Order Hold' */
} B_ACDS443_Prop_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T TimeStampA;                   /* '<S6>/Derivative' */
  real_T LastUAtTimeA[3];              /* '<S6>/Derivative' */
  real_T TimeStampB;                   /* '<S6>/Derivative' */
  real_T LastUAtTimeB[3];              /* '<S6>/Derivative' */
  real_T TimeStampA_i;                 /* '<S5>/Derivative' */
  real_T LastUAtTimeA_f[3];            /* '<S5>/Derivative' */
  real_T TimeStampB_j;                 /* '<S5>/Derivative' */
  real_T LastUAtTimeB_a[3];            /* '<S5>/Derivative' */
  real_T TimeStampA_j;                 /* '<S4>/Derivative' */
  real_T LastUAtTimeA_b[3];            /* '<S4>/Derivative' */
  real_T TimeStampB_m;                 /* '<S4>/Derivative' */
  real_T LastUAtTimeB_c[3];            /* '<S4>/Derivative' */
  real_T NextOutput[3];                /* '<Root>/Aerodynamic Drag' */
  real_T NextOutput_k[3];              /* '<Root>/SRP and Gravity ' */
  real_T Product2_DWORK1[9];           /* '<S8>/Product2' */
  real_T Product2_DWORK3[9];           /* '<S8>/Product2' */
  real_T Product2_DWORK4[9];           /* '<S8>/Product2' */
  real_T Product2_DWORK5[9];           /* '<S8>/Product2' */
  struct {
    void *LoggedData;
  } BodyEulerAngles_PWORK;             /* '<Root>/Body Euler Angles' */

  struct {
    void *LoggedData;
  } BodySpinRates_PWORK;               /* '<Root>/Body Spin Rates' */

  struct {
    void *LoggedData;
  } StoredMomentum_PWORK;              /* '<Root>/Stored Momentum' */

  struct {
    void *LoggedData;
  } Wheel1ControlMoment_PWORK;         /* '<Root>/Wheel 1Control Moment' */

  struct {
    void *LoggedData;
  } Wheel2ControlMoment_PWORK;         /* '<Root>/Wheel 2 Control Moment ' */

  struct {
    void *LoggedData;
  } Wheel3ControlMoment_PWORK;         /* '<Root>/Wheel 3 Control Moment' */

  int32_T PreLookUpIndexSearch_DWORK1; /* '<S23>/PreLook-Up Index Search' */
  int32_T Sum_DWORK1;                  /* '<S25>/Sum' */
  int32_T Product2_DWORK2[3];          /* '<S8>/Product2' */
  uint32_T RandSeed[3];                /* '<Root>/Aerodynamic Drag' */
  uint32_T RandSeed_c[3];              /* '<Root>/SRP and Gravity ' */
} DW_ACDS443_Prop_T;

/* Continuous states (default storage) */
typedef struct {
  real_T phi[3];                       /* '<S7>/phi theta psi' */
  real_T p[3];                         /* '<S3>/p,q,r ' */
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T Integrator_CSTATE_i;          /* '<S5>/Integrator' */
  real_T Integrator_CSTATE_k;          /* '<S6>/Integrator' */
  real_T Integrator1_CSTATE[3];        /* '<S6>/Integrator1' */
  real_T Integrator1_CSTATE_f[3];      /* '<S5>/Integrator1' */
  real_T Integrator1_CSTATE_a[3];      /* '<S4>/Integrator1' */
  real_T Integrator_CSTATE_j[3];       /* '<Root>/Integrator' */
  real_T mass;                         /* '<S9>/mass  ' */
  real_T U[3];                         /* '<S3>/ub,vb,wb' */
  real_T Xe[3];                        /* '<S3>/xe,ye,ze' */
} X_ACDS443_Prop_T;

/* Periodic continuous state vector (global) */
typedef int_T PeriodicIndX_ACDS443_Prop_T[3];
typedef real_T PeriodicRngX_ACDS443_Prop_T[6];

/* State derivatives (default storage) */
typedef struct {
  real_T phi[3];                       /* '<S7>/phi theta psi' */
  real_T p[3];                         /* '<S3>/p,q,r ' */
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T Integrator_CSTATE_i;          /* '<S5>/Integrator' */
  real_T Integrator_CSTATE_k;          /* '<S6>/Integrator' */
  real_T Integrator1_CSTATE[3];        /* '<S6>/Integrator1' */
  real_T Integrator1_CSTATE_f[3];      /* '<S5>/Integrator1' */
  real_T Integrator1_CSTATE_a[3];      /* '<S4>/Integrator1' */
  real_T Integrator_CSTATE_j[3];       /* '<Root>/Integrator' */
  real_T mass;                         /* '<S9>/mass  ' */
  real_T U[3];                         /* '<S3>/ub,vb,wb' */
  real_T Xe[3];                        /* '<S3>/xe,ye,ze' */
} XDot_ACDS443_Prop_T;

/* State disabled  */
typedef struct {
  boolean_T phi[3];                    /* '<S7>/phi theta psi' */
  boolean_T p[3];                      /* '<S3>/p,q,r ' */
  boolean_T Integrator_CSTATE;         /* '<S4>/Integrator' */
  boolean_T Integrator_CSTATE_i;       /* '<S5>/Integrator' */
  boolean_T Integrator_CSTATE_k;       /* '<S6>/Integrator' */
  boolean_T Integrator1_CSTATE[3];     /* '<S6>/Integrator1' */
  boolean_T Integrator1_CSTATE_f[3];   /* '<S5>/Integrator1' */
  boolean_T Integrator1_CSTATE_a[3];   /* '<S4>/Integrator1' */
  boolean_T Integrator_CSTATE_j[3];    /* '<Root>/Integrator' */
  boolean_T mass;                      /* '<S9>/mass  ' */
  boolean_T U[3];                      /* '<S3>/ub,vb,wb' */
  boolean_T Xe[3];                     /* '<S3>/xe,ye,ze' */
} XDis_ACDS443_Prop_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Parameters (default storage) */
struct P_ACDS443_Prop_T_ {
  real_T MomentumWheel3_NormalVector[3];/* Mask Parameter: MomentumWheel3_NormalVector
                                         * Referenced by: '<S6>/Wheel Orientation'
                                         */
  real_T MomentumWheel2_NormalVector[3];/* Mask Parameter: MomentumWheel2_NormalVector
                                         * Referenced by: '<S5>/Wheel Orientation'
                                         */
  real_T MomentumWheel1_NormalVector[3];/* Mask Parameter: MomentumWheel1_NormalVector
                                         * Referenced by: '<S4>/Wheel Orientation'
                                         */
  real_T SimpleVariableMass6DOFEulerAngl[3];/* Mask Parameter: SimpleVariableMass6DOFEulerAngl
                                             * Referenced by: '<S3>/ub,vb,wb'
                                             */
  real_T SimpleVariableMass6DOFEulerAn_k[3];/* Mask Parameter: SimpleVariableMass6DOFEulerAn_k
                                             * Referenced by: '<S7>/phi theta psi'
                                             */
  real_T SimpleVariableMass6DOFEulerAn_d[9];/* Mask Parameter: SimpleVariableMass6DOFEulerAn_d
                                             * Referenced by: '<S23>/slope'
                                             */
  real_T SimpleVariableMass6DOFEulerAn_j[9];/* Mask Parameter: SimpleVariableMass6DOFEulerAn_j
                                             * Referenced by: '<S23>/slope'
                                             */
  real_T SimpleVariableMass6DOFEulerAn_g;/* Mask Parameter: SimpleVariableMass6DOFEulerAn_g
                                          * Referenced by: '<S9>/mass  '
                                          */
  real_T SimpleVariableMass6DOFEulerAn_l;/* Mask Parameter: SimpleVariableMass6DOFEulerAn_l
                                          * Referenced by:
                                          *   '<S9>/mass  '
                                          *   '<S23>/Constant'
                                          *   '<S23>/slope'
                                          */
  real_T SimpleVariableMass6DOFEulerAn_i;/* Mask Parameter: SimpleVariableMass6DOFEulerAn_i
                                          * Referenced by:
                                          *   '<S9>/mass  '
                                          *   '<S23>/Constant1'
                                          *   '<S23>/slope'
                                          */
  real_T InterpolateInertia_matrix[18];/* Mask Parameter: InterpolateInertia_matrix
                                        * Referenced by:
                                        *   '<S25>/[0]'
                                        *   '<S25>/[1]'
                                        */
  real_T SimpleVariableMass6DOFEulerAn_a[3];/* Mask Parameter: SimpleVariableMass6DOFEulerAn_a
                                             * Referenced by: '<S3>/p,q,r '
                                             */
  real_T SimpleVariableMass6DOFEulerAn_e[3];/* Mask Parameter: SimpleVariableMass6DOFEulerAn_e
                                             * Referenced by: '<S3>/xe,ye,ze'
                                             */
  real_T phithetapsi_WrappedStateUpperVa;/* Expression: pi
                                          * Referenced by: '<S7>/phi theta psi'
                                          */
  real_T phithetapsi_WrappedStateLowerVa;/* Expression: -pi
                                          * Referenced by: '<S7>/phi theta psi'
                                          */
  real_T Integrator_IC;                /* Expression: .037
                                        * Referenced by: '<S4>/Integrator'
                                        */
  real_T Integrator_IC_p;              /* Expression: 0
                                        * Referenced by: '<S5>/Integrator'
                                        */
  real_T Integrator_IC_l;              /* Expression: 0
                                        * Referenced by: '<S6>/Integrator'
                                        */
  real_T Integrator1_IC;               /* Expression: 0
                                        * Referenced by: '<S6>/Integrator1'
                                        */
  real_T MomentumSaturation_UpperSat;  /* Expression: .05
                                        * Referenced by: '<S6>/Momentum Saturation'
                                        */
  real_T MomentumSaturation_LowerSat;  /* Expression: -.05
                                        * Referenced by: '<S6>/Momentum Saturation'
                                        */
  real_T MomentInputSaturation_UpperSat;/* Expression: 20e-3
                                         * Referenced by: '<S6>/Moment Input Saturation'
                                         */
  real_T MomentInputSaturation_LowerSat;/* Expression: -20e-3
                                         * Referenced by: '<S6>/Moment Input Saturation'
                                         */
  real_T Integrator1_IC_g;             /* Expression: 0
                                        * Referenced by: '<S5>/Integrator1'
                                        */
  real_T MomentumSaturation_UpperSat_b;/* Expression: .05
                                        * Referenced by: '<S5>/Momentum Saturation'
                                        */
  real_T MomentumSaturation_LowerSat_l;/* Expression: -.05
                                        * Referenced by: '<S5>/Momentum Saturation'
                                        */
  real_T MomentInputSaturation_UpperSa_e;/* Expression: 20e-3
                                          * Referenced by: '<S5>/Moment Input Saturation'
                                          */
  real_T MomentInputSaturation_LowerSa_o;/* Expression: -20e-3
                                          * Referenced by: '<S5>/Moment Input Saturation'
                                          */
  real_T Integrator1_IC_i;             /* Expression: 0
                                        * Referenced by: '<S4>/Integrator1'
                                        */
  real_T MomentumSaturation_UpperSat_g;/* Expression: .05
                                        * Referenced by: '<S4>/Momentum Saturation'
                                        */
  real_T MomentumSaturation_LowerSat_m;/* Expression: -.05
                                        * Referenced by: '<S4>/Momentum Saturation'
                                        */
  real_T MomentInputSaturation_UpperSa_h;/* Expression: 20e-3
                                          * Referenced by: '<S4>/Moment Input Saturation'
                                          */
  real_T MomentInputSaturation_LowerSa_a;/* Expression: -20e-3
                                          * Referenced by: '<S4>/Moment Input Saturation'
                                          */
  real_T Saturation1_Value;            /* Expression: -0.04
                                        * Referenced by: '<S1>/Saturation1'
                                        */
  real_T Constant3_Value[3];           /* Expression: ones(1,3)*5e-3
                                        * Referenced by: '<S1>/Constant3'
                                        */
  real_T Saturation_Value;             /* Expression: 0.04
                                        * Referenced by: '<S1>/Saturation'
                                        */
  real_T Constant2_Value[3];           /* Expression: -ones(1,3)*5e-3
                                        * Referenced by: '<S1>/Constant2'
                                        */
  real_T AerodynamicDrag_Mean[3];      /* Expression: [0 0 0]
                                        * Referenced by: '<Root>/Aerodynamic Drag'
                                        */
  real_T AerodynamicDrag_StdDev[3];    /* Computed Parameter: AerodynamicDrag_StdDev
                                        * Referenced by: '<Root>/Aerodynamic Drag'
                                        */
  real_T AerodynamicDrag_Seed[3];      /* Expression: 1:3
                                        * Referenced by: '<Root>/Aerodynamic Drag'
                                        */
  real_T Constant_Value[3];            /* Expression: [0 0 0]
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Gain1_Gain;                   /* Expression: .1
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T Integrator_IC_d;              /* Expression: 0
                                        * Referenced by: '<Root>/Integrator'
                                        */
  real_T Gain2_Gain;                   /* Expression: 3.3754e-06

                                        * Referenced by: '<Root>/Gain2'
                                        */
  real_T SRPandGravity_Mean[3];        /* Expression: [0 0 0]
                                        * Referenced by: '<Root>/SRP and Gravity '
                                        */
  real_T SRPandGravity_StdDev[3];      /* Computed Parameter: SRPandGravity_StdDev
                                        * Referenced by: '<Root>/SRP and Gravity '
                                        */
  real_T SRPandGravity_Seed[3];        /* Expression: 1:3
                                        * Referenced by: '<Root>/SRP and Gravity '
                                        */
  int32_T offsetforupperindex1_Value;  /* Computed Parameter: offsetforupperindex1_Value
                                        * Referenced by: '<S25>/offset for upper index 1'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_ACDS443_Prop_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_ACDS443_Prop_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[28];
  real_T odeF[3][28];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block parameters (default storage) */
extern P_ACDS443_Prop_T ACDS443_Prop_P;

/* Block signals (default storage) */
extern B_ACDS443_Prop_T ACDS443_Prop_B;

/* Continuous states (default storage) */
extern X_ACDS443_Prop_T ACDS443_Prop_X;

/* Block states (default storage) */
extern DW_ACDS443_Prop_T ACDS443_Prop_DW;

/* Model entry point functions */
extern void ACDS443_Prop_initialize(void);
extern void ACDS443_Prop_step(void);
extern void ACDS443_Prop_terminate(void);

/* Real-time Model object */
extern RT_MODEL_ACDS443_Prop_T *const ACDS443_Prop_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/Mass Propellant used' : Unused code path elimination
 * Block '<S11>/Unit Conversion' : Unused code path elimination
 * Block '<S12>/Unit Conversion' : Unused code path elimination
 * Block '<Root>/Sum2' : Unused code path elimination
 * Block '<Root>/Total Time Fired' : Unused code path elimination
 * Block '<S17>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S19>/Reshape1' : Reshape block reduction
 * Block '<S19>/Reshape2' : Reshape block reduction
 * Block '<S20>/Reshape1' : Reshape block reduction
 * Block '<S20>/Reshape2' : Reshape block reduction
 * Block '<S8>/Reshape' : Reshape block reduction
 * Block '<S8>/Reshape1' : Reshape block reduction
 * Block '<S13>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S14>/Reshape1' : Reshape block reduction
 * Block '<S14>/Reshape2' : Reshape block reduction
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'ACDS443_Prop'
 * '<S1>'   : 'ACDS443_Prop/Ideal Propulsion System'
 * '<S2>'   : 'ACDS443_Prop/Orthogonal Reaction Wheels'
 * '<S3>'   : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)'
 * '<S4>'   : 'ACDS443_Prop/Orthogonal Reaction Wheels/Momentum Wheel 1'
 * '<S5>'   : 'ACDS443_Prop/Orthogonal Reaction Wheels/Momentum Wheel 2'
 * '<S6>'   : 'ACDS443_Prop/Orthogonal Reaction Wheels/Momentum Wheel 3'
 * '<S7>'   : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate DCM & Euler Angles'
 * '<S8>'   : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot'
 * '<S9>'   : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Determine Force,  Mass & Inertia'
 * '<S10>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Vbxw'
 * '<S11>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Velocity Conversion'
 * '<S12>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Velocity Conversion1'
 * '<S13>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Velocity Conversion2'
 * '<S14>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/transform to Inertial axes '
 * '<S15>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S16>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate DCM & Euler Angles/phidot thetadot psidot'
 * '<S17>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S18>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product'
 * '<S19>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot/I x w'
 * '<S20>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot/I x w1'
 * '<S21>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product/Subsystem'
 * '<S22>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product/Subsystem1'
 * '<S23>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Determine Force,  Mass & Inertia/Estimate Inertia Tensor'
 * '<S24>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Determine Force,  Mass & Inertia/Estimate Inertia Tensor/Interpolate Inertia '
 * '<S25>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Determine Force,  Mass & Inertia/Estimate Inertia Tensor/Interpolate Inertia /Matrix interpolation'
 * '<S26>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Vbxw/Subsystem'
 * '<S27>'  : 'ACDS443_Prop/Simple Variable Mass 6DOF (Euler Angles)/Vbxw/Subsystem1'
 */
#endif                                 /* RTW_HEADER_ACDS443_Prop_h_ */
