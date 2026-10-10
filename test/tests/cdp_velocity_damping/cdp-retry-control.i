# 3D formula diagnostic; not a dam model or an Abaqus reference.
!include common.i
!include cdp-material.i
!include damage-observables.i
!include rate-observables.i

[Physics/SolidMechanics/CDPVelocityDynamic/all]
  add_variables = true
  strain = SMALL
  incremental = true
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.7
  stiffness_damping_coefficient = 0.03
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
  enable_cdp_velocity_damping = true
  diagnostic_3d_model = true
  expected_youngs_modulus = 3.04e10
  expected_poissons_ratio = 0.2
  expected_stiffness_damping = 0.03
  require_cdp_material = true
[]

[Executioner]
  [TimeStepper]
    type = GroundMotionValidationStepper
    validation_only = true
    time_sequence = '0 0.005 0.01 0.015 0.02 0.025 0.03 0.0325 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.105 0.11 0.115 0.12'
    reject_time = -1
  []
[]
