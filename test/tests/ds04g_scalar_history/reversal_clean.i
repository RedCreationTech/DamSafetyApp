!include reversal_on.i
[Functions/path_x]
  x = '0 .002 .004 .006 .007'
  y = '0 .00000020 -.00000004 .00000020 .00000020'
[]
[Functions/path_y]
  x = '0 .002 .004 .006 .007'
  y = '0 -.00000005 .00000001 -.00000005 -.00000005'
[]
[Executioner]
  end_time = .007
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '0 .002 .003 .005 .007'
  []
[]
