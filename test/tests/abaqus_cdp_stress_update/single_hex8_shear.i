!include single_hex8_common.i

[Functions]
  [load_path]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 5e-4'
  []
[]

[BCs]
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [load_top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top
    function = load_path
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
  file_base = single_hex8_shear_out
[]
