[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX8
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = SMALL
  add_variables = true
  generate_output = 'stress_xx strain_xx'
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [pull_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 1e-4
  []
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 3.04e10
    poissons_ratio = 0.2
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Postprocessors]
  [average_stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [average_strain_xx]
    type = ElementAverageValue
    variable = strain_xx
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
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
  file_base = elastic_out
[]
