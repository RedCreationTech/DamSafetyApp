!include single_hex8_common.i

# P1 restart capability check, phase 1 of the two-phase run.
# Stops at the same time where the load path kinks, and writes a MOOSE checkpoint so
# that phase 2 can continue from the restored state instead of re-solving this half.

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
  end_time = 0.5
  dt = 0.05
[]

[Outputs]
  exodus = false
  csv = true
  execute_on = 'initial timestep_end'
  file_base = restart_phase1_out

  [checkpoint]
    type = Checkpoint
    file_base = restart_phase1
  []
[]
