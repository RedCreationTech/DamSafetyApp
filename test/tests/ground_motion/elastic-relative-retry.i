# Small free-solve elastic validation; all y constrained, x bottom fixed relative to ground.
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmin = 0
  xmax = 2
  ymin = 0
  ymax = 1
[]
[Problem]
  type = FEProblem
[]
[Variables]
  [q_x]
  []
  [q_y]
  []
[]
[AuxVariables]
  [v_x]
  []
  [a_x]
    initial_condition = 0.004
  []
  [u_absolute]
  []
  [v_absolute]
  []
  [a_absolute]
  []
  [resid_x]
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
  [base_relative_acceleration]
    type = ConstantIC
    variable = a_x
    boundary = bottom
    value = 0
  []
[]
[Kernels]
  [stress_x]
    type = DynamicStressDivergenceTensors
    variable = q_x
    displacements = 'q_x q_y'
    component = 0
    alpha = -0.05
    zeta = 0.03
    save_in = resid_x
  []
  [stress_y]
    type = StressDivergenceTensors
    variable = q_y
    displacements = 'q_x q_y'
    component = 1
  []
  [inertia_x]
    type = InertialForce
    variable = q_x
    velocity = v_x
    acceleration = a_x
    beta = 0.275625
    gamma = 0.55
    alpha = -0.05
    eta = 0.7
    save_in = resid_x
  []
  [ground_load]
    type = GroundMotionInertia
    variable = q_x
    ground_acceleration = ground_a
    ground_velocity = ground_v
    alpha = -0.05
    eta = 0.7
    save_in = resid_x
  []
  [mass_correction]
    type = ADQuad4ConsistentMassCorrection
    variable = q_x
    velocity = v_x
    acceleration = a_x
    beta = 0.275625
    gamma = 0.55
    alpha = -0.05
    eta = 0.7
    ground_acceleration = ground_a
    ground_velocity = ground_v
    save_in = resid_x
  []
  [hourglass]
    type = ADQuad4HourglassControl
    variable = q_x
    coefficient = 1
    shear_modulus = 41.6666666666666667
    alpha = -0.05
    zeta = 0.03
    save_in = resid_x
  []
[]
[AuxKernels]
  [relative_acceleration]
    type = NewmarkAccelAux
    variable = a_x
    displacement = q_x
    velocity = v_x
    beta = 0.275625
    execute_on = TIMESTEP_END
  []
  [relative_velocity]
    type = NewmarkVelAux
    variable = v_x
    acceleration = a_x
    gamma = 0.55
    execute_on = TIMESTEP_END
  []
  [u_absolute]
    type = AbsoluteGroundMotionAux
    variable = u_absolute
    relative = q_x
    ground = ground_u
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_absolute]
    type = AbsoluteGroundMotionAux
    variable = v_absolute
    relative = v_x
    ground = ground_v
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_absolute]
    type = AbsoluteGroundMotionAux
    variable = a_absolute
    relative = a_x
    ground = ground_a
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[BCs]
  [base_x]
    type = DirichletBC
    variable = q_x
    boundary = bottom
    value = 0
  []
  [all_y]
    type = DirichletBC
    variable = q_y
    boundary = 'bottom top left right'
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
  [strain]
    type = ComputeIncrementalStrain
    displacements = 'q_x q_y'
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]
[Postprocessors]
  [q_x_0]
    type = NodalVariableValue
    variable = q_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_x_0]
    type = NodalVariableValue
    variable = v_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_0]
    type = NodalVariableValue
    variable = a_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [u_absolute_0]
    type = NodalVariableValue
    variable = u_absolute
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_absolute_0]
    type = NodalVariableValue
    variable = v_absolute
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_absolute_0]
    type = NodalVariableValue
    variable = a_absolute
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [resid_x_0]
    type = NodalVariableValue
    variable = resid_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_1]
    type = NodalVariableValue
    variable = q_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_x_1]
    type = NodalVariableValue
    variable = v_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_1]
    type = NodalVariableValue
    variable = a_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [u_absolute_1]
    type = NodalVariableValue
    variable = u_absolute
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_absolute_1]
    type = NodalVariableValue
    variable = v_absolute
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_absolute_1]
    type = NodalVariableValue
    variable = a_absolute
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [resid_x_1]
    type = NodalVariableValue
    variable = resid_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_2]
    type = NodalVariableValue
    variable = q_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_x_2]
    type = NodalVariableValue
    variable = v_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_2]
    type = NodalVariableValue
    variable = a_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [u_absolute_2]
    type = NodalVariableValue
    variable = u_absolute
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_absolute_2]
    type = NodalVariableValue
    variable = v_absolute
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_absolute_2]
    type = NodalVariableValue
    variable = a_absolute
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [resid_x_2]
    type = NodalVariableValue
    variable = resid_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_3]
    type = NodalVariableValue
    variable = q_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_x_3]
    type = NodalVariableValue
    variable = v_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_x_3]
    type = NodalVariableValue
    variable = a_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [u_absolute_3]
    type = NodalVariableValue
    variable = u_absolute
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [v_absolute_3]
    type = NodalVariableValue
    variable = v_absolute
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_absolute_3]
    type = NodalVariableValue
    variable = a_absolute
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [resid_x_3]
    type = NodalVariableValue
    variable = resid_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ground_displacement]
    type = FunctionValuePostprocessor
    function = ground_u
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ground_velocity]
    type = FunctionValuePostprocessor
    function = ground_v
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ground_acceleration]
    type = FunctionValuePostprocessor
    function = ground_a
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accepted_dt]
    type = TimestepSize
  []
[]
[Executioner]
  type = Transient
  [Quadrature]
    type = GAUSS
    order = FIRST
  []
  solve_type = NEWTON
  end_time = 0.2
  dtmin = 1e-8
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
  l_tol = 1e-12
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  [TimeStepper]
    type = GroundMotionValidationStepper
    validation_only = true
    time_sequence = '0 0.025 0.05 0.08 0.11 0.16 0.2'
    reject_time = 0.08
  []
[]
[Outputs]
  console = false
  [validation_console]
    type = Console
  []
  [field]
    type = Exodus
  []
  [history]
    type = CSV
    precision = 17
  []
[]
