# Two unit bars of stiffness1000, F_right=1, u0=0, u2=.5*u0+.5*u1.
# Independent stationary energy gives u1=.0004, u2=.0002, lambda=1.2.
[Mesh]
  [line]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 2
    xmax = 2
  []
  [anchors]
    type = AqueductConstraintNodesGenerator
    input = line
    nodes = 2
    block_id = 10
    block_name = constraint_anchors
  []
  allow_renumbering = false
[]
[Problem]
  kernel_coverage_check = SKIP_LIST
  kernel_coverage_block_list = 10
  material_coverage_check = SKIP_LIST
  material_coverage_block_list = 10
[]
[Variables]
  [u]
    block = 0
  []
  [lambda]
    block = constraint_anchors
  []
[]
[Kernels]
  [multiplier]
    type = NullKernel
    variable = lambda
    block = constraint_anchors
  []
  [bar]
    type = MatDiffusion
    variable = u
    block = 0
    diffusivity = stiffness
  []
[]
[Materials]
  [bar]
    type = GenericConstantMaterial
    block = 0
    prop_names = stiffness
    prop_values = 1000
  []
[]
[NodalKernels]
  [load]
    type = UserForcingFunctorNodalKernel
    variable = u
    boundary = right
    functor = 1
  []
[]
[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
[]
[Constraints]
  [interpolation]
    type = AqueductLinearMPC
    variable = u
    anchor_node = 2
    term_nodes = '2 0 1'
    terms = 'u u u'
    coefficients = '1 -.5 -.5'
    multiplier = lambda
  []
[]
[Postprocessors]
  [middle_u]
    type = NodalVariableValue
    variable = u
    nodeid = 1
  []
  [right_u]
    type = NodalVariableValue
    variable = u
    nodeid = 2
  []
  [multiplier]
    type = NodalVariableValue
    variable = lambda
    nodeid = 2
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = lu
[]
[Outputs]
  csv = true
  exodus = true
[]
