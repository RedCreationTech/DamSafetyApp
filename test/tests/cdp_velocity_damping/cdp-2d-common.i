# 2D damaged QUAD4 formula diagnostic, not the dam or an Abaqus counterpart.
[GlobalParams]
  displacements = 'q_x q_y'
  out_of_plane_strain = strain_zz
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  elem_type = QUAD4
[]
[Variables]
  [q_x]
  []
  [q_y]
  []
  [strain_zz]
    order = FIRST
  []
[]
[AuxVariables]
  [vel_x]
  []
  [vel_y]
  []
  [accel_x]
  []
  [accel_y]
  []
[]
[Functions/load]
  type = PiecewiseLinear
  x = '0 0.02 0.04 0.08 0.12'
  y = '0 0.0002 0 -0.0015 0'
[]
[BCs]
  [fix_x]
    type = DirichletBC
    variable = q_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = q_y
    boundary = bottom
    value = 0
  []
  [load_x]
    type = FunctionDirichletBC
    variable = q_x
    boundary = right
    function = load
  []
[]
[Physics/SolidMechanics/CDPVelocityDynamic/all]
  add_variables = false
  velocities = 'vel_x vel_y'
  accelerations = 'accel_x accel_y'
  strain = SMALL
  incremental = true
  planar_formulation = WEAK_PLANE_STRESS
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 1.95
  stiffness_damping_coefficient = 0.00113
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_xy'
  enable_cdp_velocity_damping = false
  diagnostic_2d_plane_stress_projection = true
  expected_youngs_modulus = 3.04e10
  expected_poissons_ratio = 0.2
  expected_stiffness_damping = 0.00113
  require_cdp_material = true
[]
[Materials]
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2260.74
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 3.04e10
    poissons_ratio = 0.2
  []
[]
!include cdp-material.i
!include damage-observables.i
[Materials/cdp]
  maximum_substeps = 256
  compression_recovery = 1
  tension_recovery = 1
  use_stable_principal_stress = false
[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 0.12
  dt = 0.005
  dtmin = 5e-7
  dtmax = 0.01
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11
  timestep_tolerance = 1e-9
  automatic_scaling = true
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  petsc_options = '-snes_test_jacobian'
  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]
[Outputs]
  exodus = true
  [history]
    type = CSV
    precision = 17
  []
[]

[Postprocessors]
  [stress_xx_mean]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_yy_mean]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_zz_mean]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_xy_mean]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xx_mean]
    type = ElementAverageValue
    variable = strain_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_yy_mean]
    type = ElementAverageValue
    variable = strain_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_mean]
    type = ElementAverageValue
    variable = strain_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xy_mean]
    type = ElementAverageValue
    variable = strain_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_0]
    type = NodalVariableValue
    variable = q_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_0]
    type = NodalVariableValue
    variable = q_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_0]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_0]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_0]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_0]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_1]
    type = NodalVariableValue
    variable = q_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_1]
    type = NodalVariableValue
    variable = q_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_1]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_1]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_1]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_1]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_2]
    type = NodalVariableValue
    variable = q_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_2]
    type = NodalVariableValue
    variable = q_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_2]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_2]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_2]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_2]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_3]
    type = NodalVariableValue
    variable = q_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_3]
    type = NodalVariableValue
    variable = q_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_3]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_3]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_3]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_3]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

!include cdp-2d-state-observables.i
