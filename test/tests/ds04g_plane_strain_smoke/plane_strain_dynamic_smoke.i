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

# a(t) = 2*sin(2*pi*2*t) m/s^2, 0.02 s 采样（与参考波表同间隔）
[Functions/base_accel]
  type = PiecewiseLinear
  x = '0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58 0.60 0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80 0.82 0.84 0.86 0.88 0.90 0.92 0.94 0.96 0.98 1.00'
  y = '0 0.49738 0.963507 1.36909 1.68866 1.90211 1.99605 1.96457 1.80965 1.54103 1.17557 0.736249 0.250666 -0.250666 -0.736249 -1.17557 -1.54103 -1.80965 -1.96457 -1.99605 -1.90211 -1.68866 -1.36909 -0.963507 -0.49738 -4.89859e-16 0.49738 0.963507 1.36909 1.68866 1.90211 1.99605 1.96457 1.80965 1.54103 1.17557 0.736249 0.250666 -0.250666 -0.736249 -1.17557 -1.54103 -1.80965 -1.96457 -1.99605 -1.90211 -1.68866 -1.36909 -0.963507 -0.49738 -9.79717e-16'
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
    beta = 0.25
  []
[]

[Physics/SolidMechanics/Dynamic/all]
  add_variables = true
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  newmark_beta = 0.25
  newmark_gamma = 0.5
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
  file_base = plane_strain_dynamic_smoke_out
[]
