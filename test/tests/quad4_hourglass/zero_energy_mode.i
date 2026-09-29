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

[Kernels/hourglass]
  type = ADQuad4HourglassControl
  variable = u
  coefficient = 1
  shear_modulus = 10
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
  type = NodalVariableValue
  variable = u
  nodeid = 3
[]

[Preconditioning/smp]
  type = SMP
  full = true
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-12
  [Quadrature]
    type = GAUSS
    order = FIRST
  []
[]

[Outputs]
  csv = true
[]
