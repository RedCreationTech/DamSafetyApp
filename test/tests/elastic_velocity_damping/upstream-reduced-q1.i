# Diagnostic only: underintegrated Q1 thickness has null modes; do not claim a pass from RunApp exit alone.
# Elastic velocity-rate diagnostic, NOT a dam model or an Abaqus reference.
[GlobalParams]
  displacements = 'q_x q_y'
  out_of_plane_strain = strain_zz
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = 2
  ymax = 1
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
  [v_x]
  []
  [v_y]
  []
  [a_x]
    initial_condition = 0.004
  []
  [a_y]
  []

[]
[Functions]
  [ground_u]
    type = IntegratedAcceleration
    data_file = acceleration.csv
    component = displacement
  []
  [ground_v]
    type = IntegratedAcceleration
    data_file = acceleration.csv
    component = velocity
  []
  [ground_a]
    type = IntegratedAcceleration
    data_file = acceleration.csv
    component = acceleration
  []
[]
[ICs]
  [base_initial_acceleration]
    type = ConstantIC
    variable = a_x
    boundary = bottom
    value = 0
  []
[]
[Physics/SolidMechanics/Dynamic/elastic]
  add_variables = false
  velocities = 'v_x v_y'
  accelerations = 'a_x a_y'
  strain = SMALL
  incremental = true
  planar_formulation = WEAK_PLANE_STRESS
  out_of_plane_strain = strain_zz
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.7
  stiffness_damping_coefficient = 0.03
  generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy'

[]
[Kernels]
  [base_load]
    type = GroundMotionInertia
    variable = q_x
    alpha = -0.05
    eta = 0.7
    ground_velocity = ground_v
    ground_acceleration = ground_a
  []
[]

[BCs]
  [base_x]
    type = DirichletBC
    variable = q_x
    boundary = bottom
    value = 0
  []
  [base_y]
    type = DirichletBC
    variable = q_y
    boundary = bottom
    value = 0
  []
[]
[Materials]
  [density_material]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2.5
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 100
    poissons_ratio = 0.2
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]
[Postprocessors]
  [accepted_dt]
    type = TimestepSize
  []
  [ground_acceleration]
    type = FunctionValuePostprocessor
    function = ground_a
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
  [v_x_0]
    type = NodalVariableValue
    variable = v_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_y_0]
    type = NodalVariableValue
    variable = v_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_0]
    type = NodalVariableValue
    variable = a_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_y_0]
    type = NodalVariableValue
    variable = a_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_0]
    type = NodalVariableValue
    variable = strain_zz
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
  [v_x_1]
    type = NodalVariableValue
    variable = v_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_y_1]
    type = NodalVariableValue
    variable = v_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_1]
    type = NodalVariableValue
    variable = a_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_y_1]
    type = NodalVariableValue
    variable = a_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_1]
    type = NodalVariableValue
    variable = strain_zz
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
  [v_x_2]
    type = NodalVariableValue
    variable = v_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_y_2]
    type = NodalVariableValue
    variable = v_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_2]
    type = NodalVariableValue
    variable = a_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_y_2]
    type = NodalVariableValue
    variable = a_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_2]
    type = NodalVariableValue
    variable = strain_zz
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
  [v_x_3]
    type = NodalVariableValue
    variable = v_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_y_3]
    type = NodalVariableValue
    variable = v_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_3]
    type = NodalVariableValue
    variable = a_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_y_3]
    type = NodalVariableValue
    variable = a_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_3]
    type = NodalVariableValue
    variable = strain_zz
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 0.2
  dtmin = 1e-8
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
  l_tol = 1e-12
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  [Quadrature]
    type = GAUSS
    order = FIRST
  []
  dt = 0.01
[]
[Outputs]
  exodus = true
  csv = true
  [history]
    type = CSV
    precision = 17
  []
[]
