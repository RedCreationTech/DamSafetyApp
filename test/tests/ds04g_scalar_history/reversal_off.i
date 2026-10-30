!include reversal_common.i
[Functions]
  [path_x]
    type = PiecewiseLinear
    x = '0 .010 .030 .032 .040 .042 .050'
    y = '0 .00000020 .00000020 -.00000004 -.00000004 .00000020 .00000020'
  []
  [path_y]
    type = PiecewiseLinear
    x = '0 .010 .030 .032 .040 .042 .050'
    y = '0 -.00000005 -.00000005 .00000001 .00000001 -.00000005 -.00000005'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = .050
  dt = .002
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
[]
