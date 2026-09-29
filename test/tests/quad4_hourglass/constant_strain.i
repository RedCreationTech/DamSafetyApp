[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  elem_type = QUAD4
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Variables/strain_zz]
  family = MONOMIAL
  order = CONSTANT
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = SMALL
  add_variables = true
  planar_formulation = WEAK_PLANE_STRESS
  generate_output = 'stress_xx strain_xx'
[]

[Kernels]
  [hourglass_x]
    type = ADQuad4HourglassControl
    variable = disp_x
    coefficient = 1
    shear_modulus = 1.2666666666666666e10
  []
  [hourglass_y]
    type = ADQuad4HourglassControl
    variable = disp_y
    coefficient = 1
    shear_modulus = 1.2666666666666666e10
  []
[]

[BCs]
  [left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [right_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 1e-4
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
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
  [Quadrature]
    type = GAUSS
    order = FIRST
  []
[]

[Outputs]
  csv = true
[]
