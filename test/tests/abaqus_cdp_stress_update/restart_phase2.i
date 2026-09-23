!include single_hex8_common.i

# P1 restart capability check, phase 2 of the two-phase run.
# __RESTART_BASE__ is replaced by scripts/cdp/verify_checkpoint_restart.py with the
# checkpoint written by restart_phase1.i, so the only difference from the continuous
# run is how the first half was obtained. Keeping the parameter inside the input file
# (rather than "-r" on the command line) is deliberate: it is the form that the C06
# command whitelist can already accept, so no agent change is needed to run this route.

[Functions]
  [load_path]
    type = PiecewiseLinear
    x = '0 0.5 1'
    y = '0 1.6e-4 3.0e-4'
  []
[]

[Problem]
  restart_file_base = __RESTART_BASE__
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
  start_time = 0.5
[]

# end_time and dt come from single_hex8_common.i; the checkpoint restores the clock.

[Outputs]
  exodus = false
  csv = true
  execute_on = 'initial timestep_end'
  file_base = restart_phase2_out
[]
