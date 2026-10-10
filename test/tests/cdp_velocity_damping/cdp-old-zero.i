# 3D formula diagnostic; not a dam model or an Abaqus reference.
!include common.i
!include cdp-load.i
!include cdp-material.i
!include damage-observables.i

[Physics/SolidMechanics/CDPVelocityDynamic/all]
  add_variables = true
  strain = SMALL
  incremental = true
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.7
  stiffness_damping_coefficient = 0
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
  enable_cdp_velocity_damping = false
[]
