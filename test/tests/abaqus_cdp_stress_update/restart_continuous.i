!include single_hex8_common.i

# P1 restart capability check, continuous reference run.
# The same absolute-in-time load path is used by the two-phase runs, so the three
# solutions must agree step by step without any time shifting.

[Functions]
  [load_path]
    type = PiecewiseLinear
    x = '0 0.5 1'
    y = '0 1.6e-4 3.0e-4'
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

[Executioner]
  end_time = 1
  dt = 0.05
[]

[Outputs]
  exodus = false
  csv = true
  execute_on = 'initial timestep_end'
  file_base = restart_continuous_out
[]
