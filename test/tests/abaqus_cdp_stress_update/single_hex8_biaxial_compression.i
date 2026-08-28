!include single_hex8_common.i

[Functions]
  [load_path]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 -1.5e-3'
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [load_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = load_path
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [load_y]
    type = FunctionDirichletBC
    variable = disp_y
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
  file_base = single_hex8_biaxial_compression_out
[]
