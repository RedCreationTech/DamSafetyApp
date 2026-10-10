# 3D single HEX8 formula diagnostic, NOT the 2D dam or an Abaqus comparison.
[GlobalParams]
  displacements = 'q_x q_y q_z'
[]
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX8
[]
[Functions]
  [load]
    type = PiecewiseLinear
    x = '0 0.02 0.04 0.08 0.12'
    y = '0 0.0002 0 -0.0015 0'
  []
[]
[BCs]
  [fix_x]
    type = DirichletBC
    variable = q_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = q_y
    boundary = bottom
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = q_z
    boundary = back
    value = 0
  []
  [load_x]
    type = FunctionDirichletBC
    variable = q_x
    boundary = right
    function = load
  []
[]
[Materials]
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2500
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 3.04e10
    poissons_ratio = 0.2
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
  end_time = 0.12
  dt = 0.005
  dtmin = 1e-8
  nl_max_its = 40
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  l_tol = 1e-12
  automatic_scaling = true
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
[]
[Outputs]
  exodus = true
  [history]
    type = CSV
    precision = 17
  []
[]


[Postprocessors]
  [q_x_0]
    type = NodalVariableValue
    variable = q_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_0]
    type = NodalVariableValue
    variable = q_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_0]
    type = NodalVariableValue
    variable = q_z
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_0]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_0]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_0]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_0]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_0]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_0]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_1]
    type = NodalVariableValue
    variable = q_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_1]
    type = NodalVariableValue
    variable = q_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_1]
    type = NodalVariableValue
    variable = q_z
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_1]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_1]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_1]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_1]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_1]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_1]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_2]
    type = NodalVariableValue
    variable = q_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_2]
    type = NodalVariableValue
    variable = q_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_2]
    type = NodalVariableValue
    variable = q_z
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_2]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_2]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_2]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_2]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_2]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_2]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_3]
    type = NodalVariableValue
    variable = q_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_3]
    type = NodalVariableValue
    variable = q_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_3]
    type = NodalVariableValue
    variable = q_z
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_3]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_3]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_3]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_3]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_3]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_3]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_4]
    type = NodalVariableValue
    variable = q_x
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_4]
    type = NodalVariableValue
    variable = q_y
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_4]
    type = NodalVariableValue
    variable = q_z
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_4]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_4]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_4]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_4]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_4]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_4]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_5]
    type = NodalVariableValue
    variable = q_x
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_5]
    type = NodalVariableValue
    variable = q_y
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_5]
    type = NodalVariableValue
    variable = q_z
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_5]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_5]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_5]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_5]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_5]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_5]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_6]
    type = NodalVariableValue
    variable = q_x
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_6]
    type = NodalVariableValue
    variable = q_y
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_6]
    type = NodalVariableValue
    variable = q_z
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_6]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_6]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_6]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_6]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_6]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_6]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 6
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_x_7]
    type = NodalVariableValue
    variable = q_x
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_y_7]
    type = NodalVariableValue
    variable = q_y
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [q_z_7]
    type = NodalVariableValue
    variable = q_z
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_x_7]
    type = NodalVariableValue
    variable = vel_x
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_y_7]
    type = NodalVariableValue
    variable = vel_y
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [vel_z_7]
    type = NodalVariableValue
    variable = vel_z
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_x_7]
    type = NodalVariableValue
    variable = accel_x
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_y_7]
    type = NodalVariableValue
    variable = accel_y
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accel_z_7]
    type = NodalVariableValue
    variable = accel_z
    nodeid = 7
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [stress_xx_mean]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_yy_mean]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_zz_mean]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_xy_mean]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xx_mean]
    type = ElementAverageValue
    variable = strain_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_yy_mean]
    type = ElementAverageValue
    variable = strain_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz_mean]
    type = ElementAverageValue
    variable = strain_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xy_mean]
    type = ElementAverageValue
    variable = strain_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [accepted_dt]
    type = TimestepSize
  []
[]
