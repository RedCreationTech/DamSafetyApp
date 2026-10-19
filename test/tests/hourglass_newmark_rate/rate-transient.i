# Bounded kernel development fixture; not a dam or CDP validation result.
[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    elem_type = QUAD4
  []
  [node_0]
    type = ExtraNodesetGenerator
    input = base
    new_boundary = node_0
    coord = '0 0'
  []
  [node_1]
    type = ExtraNodesetGenerator
    input = node_0
    new_boundary = node_1
    coord = '1 0'
  []
  [node_2]
    type = ExtraNodesetGenerator
    input = node_1
    new_boundary = node_2
    coord = '1 1'
  []
  [node_3]
    type = ExtraNodesetGenerator
    input = node_2
    new_boundary = node_3
    coord = '0 1'
  []
[]

[Variables/u]
[]
[AuxVariables/dummy_y]
[]
[AuxVariables/velocity_u]
[]
[AuxVariables/acceleration_u]
[]
[AuxKernels/acceleration_u]
 type = NewmarkAccelAux
 variable = acceleration_u
 displacement = u
 velocity = velocity_u
 beta = 0.275625
 execute_on = TIMESTEP_END
[]
[AuxKernels/velocity_u]
 type = NewmarkVelAux
 variable = velocity_u
 acceleration = acceleration_u
 gamma = 0.55
 execute_on = TIMESTEP_END
[]

[Kernels/hourglass]
  type = ADQuad4HourglassControl
  variable = u
  coefficient = 1
  shear_modulus = 10
  zeta = 0.00113
  alpha = -0.05
  elastic_operator = gradient_raw_conjugate
  use_previous_accepted_interval = true
  velocity = velocity_u
  acceleration = acceleration_u
  newmark_beta = 0.275625
  newmark_gamma = 0.55
[]

[BCs]
  [node_0]
    type = DirichletBC
    variable = u
    boundary = node_0
    value = 0
  []
  [node_1]
    type = DirichletBC
    variable = u
    boundary = node_1
    value = 0
  []
  [node_3]
    type = DirichletBC
    variable = u
    boundary = node_3
    value = 0
  []
[]

[DiracKernels/load]
  type = ConstantPointSource
  variable = u
  point = '1 1 0'
  value = 1
[]

[Postprocessors/node_2_u]
  type = PointValue
  variable = u
  point = '1 1 0'
[]

[Preconditioning/smp]
  type = SMP
  full = true
[]



[Kernels]
  [time]
    type = ADTimeDerivative
    variable = u
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_abs_tol = 1e-12
  [Quadrature]
    type = GAUSS
    order = FIRST
  []
  start_time = 0
  end_time = 0.03
  dt = 0.01
  dtmin = 1e-8
  [TimeStepper]
    type = GroundMotionValidationStepper
    validation_only = true
    time_sequence = '0 0.01 0.02 0.03'
    reject_time = -1
  []
[]
[Outputs]
  csv = false
  [history]
    type = CSV
    precision = 17
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors/node_2_velocity]
 type = PointValue
 variable = velocity_u
 point = '1 1 0'
[]
[Postprocessors/node_2_acceleration]
 type = PointValue
 variable = acceleration_u
 point = '1 1 0'
[]
