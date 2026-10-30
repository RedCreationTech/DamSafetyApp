# Runtime fixture only: uniform unit-mass/unit-stiffness oscillator.
# Quasi-static u ramps0->1 and holds. Constant load after2s must keep u=1,v=a=0.
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]
[Variables/u]
[]
[AuxVariables/v]
[]
[AuxVariables/a]
[]
[Functions/load]
  type = ParsedFunction
  expression = 'if(t<1,t,1)'
[]
[Functions/clock]
  type = PiecewiseConstant
  direction = LEFT_INCLUSIVE
  x = '0 1 2 3'
  y = '.2 .1 .05 .05'
[]
[Materials/density]
  type = GenericConstantMaterial
  prop_names = density
  prop_values = 1
[]
[Kernels/spring]
  type = Reaction
  variable = u
[]
[Kernels/force]
  type = BodyForce
  variable = u
  value = 1
  function = load
[]
[Kernels/mass]
  type = AqueductHHTInertia
  variable = u
  velocity = v
  acceleration = a
  beta = .275625
  gamma = .55
  hht_alpha = -.05
  hht_start_time = 2
[]
[AuxKernels/acceleration]
  type = AqueductNewmarkAcceleration
  variable = a
  displacement = u
  velocity = v
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/velocity]
  type = AqueductNewmarkVelocity
  variable = v
  acceleration = a
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Postprocessors/clock]
  type = TimePostprocessor
[]
[Postprocessors/u]
  type = NodalVariableValue
  variable = u
  nodeid = 0
[]
[Postprocessors/v]
  type = NodalExtremeValue
  variable = v
  value_type = max
[]
[Postprocessors/a]
  type = NodalExtremeValue
  variable = a
  value_type = max
[]
[UserObjects/assert_equilibrium_and_zero_rates]
  type = Terminator
  expression = 'if(clock>=1,abs(u-1)>1e-9,0) | abs(v)>1e-8 | abs(a)>1e-7'
  message = 'Gravity handoff fixture failed: nonzero pseudo-static rates or lost equilibrium'
  error_level = ERROR
  execute_on = TIMESTEP_END
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 2.15
  nl_abs_tol = 1e-12
  [TimeIntegrator]
    type = AqueductHHT
    hht_alpha = -.05
    beta = .275625
    gamma = .55
    hht_start_time = 2
    quasistatic_initialization = true
    inactive_tsteps = 0
  []
  [TimeStepper]
    type = FunctionDT
    function = clock
  []
[]
[Outputs/history]
  type = AqueductScheduledCSV
  file_base = gravity_handoff_history
  phase_origin = 2
  phase_end = 2.15
  mark_interval = .05
  show = 'clock u v a'
  execute_on = 'TIMESTEP_END FINAL'
[]
