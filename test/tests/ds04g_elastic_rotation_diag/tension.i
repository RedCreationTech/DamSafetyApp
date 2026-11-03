[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmax = 0.001
    ymax = 0.001
  []
  [n1]
    type = BoundingBoxNodeSetGenerator
    input = gmg
    new_boundary = n1
    bottom_left = '-1e-8 -1e-8 -1e-8'
    top_right = '1e-8 1e-8 1e-8'
  []
  [n2]
    type = BoundingBoxNodeSetGenerator
    input = n1
    new_boundary = n2
    bottom_left = '0.00099999 -1e-8 -1e-8'
    top_right = '0.00100001 1e-8 1e-8'
  []
  [n3]
    type = BoundingBoxNodeSetGenerator
    input = n2
    new_boundary = n3
    bottom_left = '0.00099999 0.00099999 -1e-8'
    top_right = '0.00100001 0.00100001 1e-8'
  []
  [n4]
    type = BoundingBoxNodeSetGenerator
    input = n3
    new_boundary = n4
    bottom_left = '-1e-8 0.00099999 -1e-8'
    top_right = '1e-8 0.00100001 1e-8'
  []
[]
[Physics/SolidMechanics/FiniteNewmark/all]
  add_variables = true
  sp1_rotation_elastic_diagnostic = true
  enable_finite_newmark_damping = false
  spatial_discretization = sp1_center_hg
  sp1_consistent_tangent = false
  sp1_density = 2400
  sp1_static_end = 0.040001
  sp1_gravity = 0
  expected_youngs_modulus = 2.97915e10
  expected_poissons_ratio = 0.2
  expected_stiffness_damping = 0.003604
  velocities = 'vel_x vel_y'
  accelerations = 'accel_x accel_y'
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.05
  stiffness_damping_coefficient = 0.003604
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
[]
[Kernels]
  [static_x]
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [static_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [gravity_y]
    type = Gravity
    variable = disp_y
    value = 0
    alpha = -0.05
    use_displaced_mesh = true
  []
[]
[Controls/phase]
  type = TimePeriod
  execute_on = 'INITIAL TIMESTEP_BEGIN'
  start_time = 0
  end_time = 0.040001
  enable_objects = '*/static_x */static_y'
  disable_objects = '*/TM_all0 */TM_all1 */TM_all_inertia_x */TM_all_inertia_y */TM_all_accel_x */TM_all_accel_y */TM_all_vel_x */TM_all_vel_y'
[]
[Materials]
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2400
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.97915e10
    poissons_ratio = 0.2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = ''
    perform_finite_strain_rotations = false
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 0.04
  dt = 0.002
  dtmin = 1e-7
[]
[VectorPostprocessors]
  [fields]
    type = ElementValueSampler
    variable = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
    sort_by = id
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Outputs]
  csv = true
  exodus = true
  file_base = mixed_elastic
[]
[Functions]
  [eps]
    type = PiecewiseLinear
    x = '0 0.01 0.02 0.03 0.04'
    y = '0 2e-7 2e-7 2e-7 2e-7'
  []
  [g]
    type = PiecewiseLinear
    x = '0 0.01 0.02 0.03 0.04'
    y = '0 0 0 0 0'
  []
  [comb]
    type = PiecewiseLinear
    x = '0 0.01 0.02 0.03 0.04'
    y = '0 2e-7 2e-7 2e-7 2e-7'
  []
  [vy]
    type = PiecewiseLinear
    x = '0 0.01 0.02 0.03 0.04'
    y = '0 -5e-8 -5e-8 -5e-8 -5e-8'
  []
[]
[BCs]
  [n1_disp_x]
    type = DirichletBC
    boundary = n1
    variable = disp_x
    value = 0
  []
  [n1_disp_y]
    type = DirichletBC
    boundary = n1
    variable = disp_y
    value = 0
  []
  [n2_disp_x]
    type = FunctionDirichletBC
    boundary = n2
    variable = disp_x
    function = eps
  []
  [n2_disp_y]
    type = DirichletBC
    boundary = n2
    variable = disp_y
    value = 0
  []
  [n3_disp_x]
    type = FunctionDirichletBC
    boundary = n3
    variable = disp_x
    function = comb
  []
  [n3_disp_y]
    type = FunctionDirichletBC
    boundary = n3
    variable = disp_y
    function = vy
  []
  [n4_disp_x]
    type = FunctionDirichletBC
    boundary = n4
    variable = disp_x
    function = g
  []
  [n4_disp_y]
    type = FunctionDirichletBC
    boundary = n4
    variable = disp_y
    function = vy
  []
[]
