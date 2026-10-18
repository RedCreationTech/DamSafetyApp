# Bounded kernel development fixture; not a dam or CDP validation result.
!include ../quad4_hourglass/zero_energy_mode.i

[Kernels]
  [time]
    type = ADTimeDerivative
    variable = u
  []
[]
[Executioner]
  type = Transient
  start_time = 0
  end_time = 0.03
  dt = 0.01
  dtmin = 1e-8
  [TimeStepper]
    type = GroundMotionValidationStepper
    validation_only = true
    time_sequence = '0 0.01 0.02 0.03'
    reject_time = -1
  []
[]
[Kernels/hourglass]
  zeta = 0.00113
  alpha = -0.05
[]
[Outputs]
  csv = false
  [history]
    type = CSV
    precision = 17
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
