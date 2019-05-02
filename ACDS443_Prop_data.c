/*
 * ACDS443_Prop_data.c
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

/* Block parameters (default storage) */
P_ACDS443_Prop_T ACDS443_Prop_P = {
  /* Mask Parameter: MomentumWheel3_NormalVector
   * Referenced by: '<S6>/Wheel Orientation'
   */
  { 0.0, 0.0, 1.0 },

  /* Mask Parameter: MomentumWheel2_NormalVector
   * Referenced by: '<S5>/Wheel Orientation'
   */
  { 0.0, 1.0, 0.0 },

  /* Mask Parameter: MomentumWheel1_NormalVector
   * Referenced by: '<S4>/Wheel Orientation'
   */
  { 1.0, 0.0, 0.0 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAngl
   * Referenced by: '<S3>/ub,vb,wb'
   */
  { 0.0, 0.0, 0.0 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_k
   * Referenced by: '<S7>/phi theta psi'
   */
  { 0.0, 0.0, 0.0 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_d
   * Referenced by: '<S23>/slope'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_j
   * Referenced by: '<S23>/slope'
   */
  { 0.0541, 0.0, 0.0, 0.0, 0.1054, 0.0, 0.0, 0.0, 0.0914 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_g
   * Referenced by: '<S9>/mass  '
   */
  1.0,

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_l
   * Referenced by:
   *   '<S9>/mass  '
   *   '<S23>/Constant'
   *   '<S23>/slope'
   */
  0.5,

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_i
   * Referenced by:
   *   '<S9>/mass  '
   *   '<S23>/Constant1'
   *   '<S23>/slope'
   */
  2.0,

  /* Mask Parameter: InterpolateInertia_matrix
   * Referenced by:
   *   '<S25>/[0]'
   *   '<S25>/[1]'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0541, 0.0, 0.0, 0.0, 0.1054,
    0.0, 0.0, 0.0, 0.0914 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_a
   * Referenced by: '<S3>/p,q,r '
   */
  { 0.174533, 0.174533, 0.174533 },

  /* Mask Parameter: SimpleVariableMass6DOFEulerAn_e
   * Referenced by: '<S3>/xe,ye,ze'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: pi
   * Referenced by: '<S7>/phi theta psi'
   */
  3.1415926535897931,

  /* Expression: -pi
   * Referenced by: '<S7>/phi theta psi'
   */
  -3.1415926535897931,

  /* Expression: .037
   * Referenced by: '<S4>/Integrator'
   */
  0.037,

  /* Expression: 0
   * Referenced by: '<S5>/Integrator'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S6>/Integrator'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S6>/Integrator1'
   */
  0.0,

  /* Expression: .05
   * Referenced by: '<S6>/Momentum Saturation'
   */
  0.05,

  /* Expression: -.05
   * Referenced by: '<S6>/Momentum Saturation'
   */
  -0.05,

  /* Expression: 20e-3
   * Referenced by: '<S6>/Moment Input Saturation'
   */
  0.02,

  /* Expression: -20e-3
   * Referenced by: '<S6>/Moment Input Saturation'
   */
  -0.02,

  /* Expression: 0
   * Referenced by: '<S5>/Integrator1'
   */
  0.0,

  /* Expression: .05
   * Referenced by: '<S5>/Momentum Saturation'
   */
  0.05,

  /* Expression: -.05
   * Referenced by: '<S5>/Momentum Saturation'
   */
  -0.05,

  /* Expression: 20e-3
   * Referenced by: '<S5>/Moment Input Saturation'
   */
  0.02,

  /* Expression: -20e-3
   * Referenced by: '<S5>/Moment Input Saturation'
   */
  -0.02,

  /* Expression: 0
   * Referenced by: '<S4>/Integrator1'
   */
  0.0,

  /* Expression: .05
   * Referenced by: '<S4>/Momentum Saturation'
   */
  0.05,

  /* Expression: -.05
   * Referenced by: '<S4>/Momentum Saturation'
   */
  -0.05,

  /* Expression: 20e-3
   * Referenced by: '<S4>/Moment Input Saturation'
   */
  0.02,

  /* Expression: -20e-3
   * Referenced by: '<S4>/Moment Input Saturation'
   */
  -0.02,

  /* Expression: -0.04
   * Referenced by: '<S1>/Saturation1'
   */
  -0.04,

  /* Expression: ones(1,3)*5e-3
   * Referenced by: '<S1>/Constant3'
   */
  { 0.005, 0.005, 0.005 },

  /* Expression: 0.04
   * Referenced by: '<S1>/Saturation'
   */
  0.04,

  /* Expression: -ones(1,3)*5e-3
   * Referenced by: '<S1>/Constant2'
   */
  { -0.005, -0.005, -0.005 },

  /* Expression: [0 0 0]
   * Referenced by: '<Root>/Aerodynamic Drag'
   */
  { 0.0, 0.0, 0.0 },

  /* Computed Parameter: AerodynamicDrag_StdDev
   * Referenced by: '<Root>/Aerodynamic Drag'
   */
  { 7.9624116949577533E-6, 7.9624116949577533E-6, 7.9624116949577533E-6 },

  /* Expression: 1:3
   * Referenced by: '<Root>/Aerodynamic Drag'
   */
  { 1.0, 2.0, 3.0 },

  /* Expression: [0 0 0]
   * Referenced by: '<Root>/Constant'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: 1
   * Referenced by: '<Root>/Gain'
   */
  1.0,

  /* Expression: .1
   * Referenced by: '<Root>/Gain1'
   */
  0.1,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator'
   */
  0.0,

  /* Expression: 3.3754e-06

   * Referenced by: '<Root>/Gain2'
   */
  3.3754E-6,

  /* Expression: [0 0 0]
   * Referenced by: '<Root>/SRP and Gravity '
   */
  { 0.0, 0.0, 0.0 },

  /* Computed Parameter: SRPandGravity_StdDev
   * Referenced by: '<Root>/SRP and Gravity '
   */
  { 0.0001140175425099138, 0.0001140175425099138, 0.0001140175425099138 },

  /* Expression: 1:3
   * Referenced by: '<Root>/SRP and Gravity '
   */
  { 1.0, 2.0, 3.0 },

  /* Computed Parameter: offsetforupperindex1_Value
   * Referenced by: '<S25>/offset for upper index 1'
   */
  1
};
