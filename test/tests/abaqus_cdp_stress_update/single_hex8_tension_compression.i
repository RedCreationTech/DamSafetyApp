!include single_hex8_common.i

[Functions]
  [load_path]
    type = PiecewiseLinear
    x = '0 0.4 1'
    y = '0 2e-4 -1.5e-3'
  []
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
  [load_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = load_path
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'initial timestep_end'
  file_base = single_hex8_tension_compression_out
[]
