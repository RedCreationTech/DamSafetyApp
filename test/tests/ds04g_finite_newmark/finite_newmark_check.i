# DS04G S2 平面应变动力冒烟：2D QUAD4 PLANE_STRAIN + 有限应变 + HHT/Newmark
# + 双系数 Rayleigh（本例值 0.05 / 0.003604）+ 基底 X 向 PresetAcceleration。
# 解析对照见 analytic_check.py（刚体一致性：顶部位移≈基底位移二重积分）。
# 不是大坝算例，不对应任何 Abaqus 参考。
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 4
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 2
  elem_type = QUAD4
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

# Synthetic function used only for input construction checks; no time advance.
[Functions/base_accel]
  type = ParsedFunction
  expression = '2*sin(12.566370614359172*t)'
[]

[BCs]
  [base_y_fixed]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [base_x_accel]
    type = PresetAcceleration
    variable = disp_x
    boundary = bottom
    function = base_accel
    velocity = vel_x
    acceleration = accel_x
    beta = 0.275625
  []
[]

[Physics/SolidMechanics/FiniteNewmark/all]
  add_variables = true
  enable_finite_newmark_damping = true
  expected_youngs_modulus = 2.97915e10
  expected_poissons_ratio = 0.2
  expected_stiffness_damping = 0.003604
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.05
  stiffness_damping_coefficient = 0.003604
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
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
    type = ComputeFiniteStrainElasticStress
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 1.0
  dt = 0.01
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
[]

[Postprocessors]
  [top_disp_x]
    type = PointValue
    variable = disp_x
    point = '0.5 2 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [top_vel_x]
    type = PointValue
    variable = vel_x
    point = '0.5 2 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [top_accel_x]
    type = PointValue
    variable = accel_x
    point = '0.5 2 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_disp_x]
    type = PointValue
    variable = disp_x
    point = '0.5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [szz_mean]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  csv = true
  exodus = true
  file_base = finite_newmark_check_out
[]
