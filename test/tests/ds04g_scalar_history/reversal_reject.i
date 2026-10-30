!include reversal_common.i
[Functions]
  [path_x]
    type = PiecewiseLinear
    x = '0 .002 .004 .006 .007'
    y = '0 .00000020 -.00000004 .00000020 .00000020'
  []
  [path_y]
    type = PiecewiseLinear
    x = '0 .002 .004 .006 .007'
    y = '0 -.00000005 .00000001 -.00000005 -.00000005'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = .007
  dt = .002
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  [TimeStepper]
    type = DS04GAcceptedIntervalTestStepper
    dt = .002
  []
[]
[Materials/cdp]
  use_combined_damage_relaxation = true
[]
