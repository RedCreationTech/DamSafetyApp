[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]
[Physics/SolidMechanics/QuasiStatic/all]
  add_variables = true
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  decomposition_method = HughesWinget
  generate_output = 'stress_xx stress_yy stress_xy'
[]
[Functions]
  [varying_stiffness]
    type = ParsedFunction
    expression = '1+t'
  []
  [ux]
    type = ParsedFunction
    expression = '0.0001*t*x'
  []
[]
[BCs]
  [x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right top bottom'
    function = ux
  []
  [y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left right top bottom'
    value = 0
  []
[]
[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.97915e10
    poissons_ratio = 0.2
  []
  [z_stress]
    type = DS04GHWSpatialElasticStress
    acknowledge_research_limits = true
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 0.5
  end_time = 1
[]
[Outputs]
  csv = true
[]

[Materials/a_extra]
  type = ComputeExtraStressConstant
  extra_stress_tensor = '1e-200 0 0 0 0 0 0 0 0'
[]
