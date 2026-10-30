# AD shell/SPOS gravity handoff and output fixture; not the full aqueduct.
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = .2
  ymax = .2
  elem_type = QUAD4
[]
[GlobalParams]
  displacements = 'ux uy uz'
  rotations = 'rx ry rz'
[]
[Variables]
[]
[Variables/ux]
[]
[Variables/uy]
[]
[Variables/uz]
[]
[Variables/rx]
[]
[Variables/ry]
[]
[Variables/rz]
[]
[Kernels]
[]
[Kernels/ux]
  type = AqueductFlatShell
  variable = ux
  component = 0
[]
[Kernels/uy]
  type = AqueductFlatShell
  variable = uy
  component = 1
[]
[Kernels/uz]
  type = AqueductFlatShell
  variable = uz
  component = 2
[]
[Kernels/rx]
  type = AqueductFlatShell
  variable = rx
  component = 3
[]
[Kernels/ry]
  type = AqueductFlatShell
  variable = ry
  component = 4
[]
[Kernels/rz]
  type = AqueductFlatShell
  variable = rz
  component = 5
[]
[Materials]
[]
[Materials/shell]
  type = AqueductFlatShellStress
  youngs_modulus = 200e9
  poissons_ratio = 0
  yield_stress = 400e6
  ultimate_stress = 696e6
  plastic_strain_at_ultimate = .1
  thickness = .1
  offset = .5
  shear_correction = .8333333333333333
  drilling_factor = .001
[]
[BCs]
[]
[BCs/ux]
  type = DirichletBC
  variable = ux
  boundary = left
  value = 0
[]
[BCs/uy]
  type = DirichletBC
  variable = uy
  boundary = left
  value = 0
[]
[BCs/uz]
  type = DirichletBC
  variable = uz
  boundary = left
  value = 0
[]
[BCs/rx]
  type = DirichletBC
  variable = rx
  boundary = left
  value = 0
[]
[BCs/ry]
  type = DirichletBC
  variable = ry
  boundary = left
  value = 0
[]
[BCs/rz]
  type = DirichletBC
  variable = rz
  boundary = left
  value = 0
[]
[Postprocessors]
[]
[Postprocessors/tip_ux]
  type = PointValue
  variable = ux
  point = '.2 .1 0'
[]
[Postprocessors/tip_uz]
  type = PointValue
  variable = uz
  point = '.2 .1 0'
[]
[Postprocessors/tip_ry]
  type = PointValue
  variable = ry
  point = '.2 .1 0'
[]
[Preconditioning]
[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Functions/gravity]
  type = ParsedFunction
  expression = 'if(t<1,t,1)'
[]
[Functions/clock]
  type = PiecewiseConstant
  direction = LEFT_INCLUSIVE
  x = '0 1 2 3'
  y = '.2 .1 .05 .05'
[]
[Kernels/gravity]
  type = BodyForce
  variable = uz
  value = -7693
  function = gravity
  use_displaced_mesh = false
[]
[AuxVariables/vux]
[]
[AuxVariables/aux]
[]
[AuxKernels/aux]
  type = AqueductNewmarkAcceleration
  variable = aux
  displacement = ux
  velocity = vux
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vux]
  type = AqueductNewmarkVelocity
  variable = vux
  acceleration = aux
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_ux]
  type = AqueductFlatShellInertia
  variable = ux
  component = 0
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vux_min]
  type = NodalExtremeValue
  variable = vux
  value_type = min
[]
[Postprocessors/vux_max]
  type = NodalExtremeValue
  variable = vux
  value_type = max
[]
[Postprocessors/aux_min]
  type = NodalExtremeValue
  variable = aux
  value_type = min
[]
[Postprocessors/aux_max]
  type = NodalExtremeValue
  variable = aux
  value_type = max
[]
[AuxVariables/vuy]
[]
[AuxVariables/auy]
[]
[AuxKernels/auy]
  type = AqueductNewmarkAcceleration
  variable = auy
  displacement = uy
  velocity = vuy
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vuy]
  type = AqueductNewmarkVelocity
  variable = vuy
  acceleration = auy
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_uy]
  type = AqueductFlatShellInertia
  variable = uy
  component = 1
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vuy_min]
  type = NodalExtremeValue
  variable = vuy
  value_type = min
[]
[Postprocessors/vuy_max]
  type = NodalExtremeValue
  variable = vuy
  value_type = max
[]
[Postprocessors/auy_min]
  type = NodalExtremeValue
  variable = auy
  value_type = min
[]
[Postprocessors/auy_max]
  type = NodalExtremeValue
  variable = auy
  value_type = max
[]
[AuxVariables/vuz]
[]
[AuxVariables/auz]
[]
[AuxKernels/auz]
  type = AqueductNewmarkAcceleration
  variable = auz
  displacement = uz
  velocity = vuz
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vuz]
  type = AqueductNewmarkVelocity
  variable = vuz
  acceleration = auz
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_uz]
  type = AqueductFlatShellInertia
  variable = uz
  component = 2
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vuz_min]
  type = NodalExtremeValue
  variable = vuz
  value_type = min
[]
[Postprocessors/vuz_max]
  type = NodalExtremeValue
  variable = vuz
  value_type = max
[]
[Postprocessors/auz_min]
  type = NodalExtremeValue
  variable = auz
  value_type = min
[]
[Postprocessors/auz_max]
  type = NodalExtremeValue
  variable = auz
  value_type = max
[]
[AuxVariables/vrx]
[]
[AuxVariables/arx]
[]
[AuxKernels/arx]
  type = AqueductNewmarkAcceleration
  variable = arx
  displacement = rx
  velocity = vrx
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vrx]
  type = AqueductNewmarkVelocity
  variable = vrx
  acceleration = arx
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_rx]
  type = AqueductFlatShellInertia
  variable = rx
  component = 3
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vrx_min]
  type = NodalExtremeValue
  variable = vrx
  value_type = min
[]
[Postprocessors/vrx_max]
  type = NodalExtremeValue
  variable = vrx
  value_type = max
[]
[Postprocessors/arx_min]
  type = NodalExtremeValue
  variable = arx
  value_type = min
[]
[Postprocessors/arx_max]
  type = NodalExtremeValue
  variable = arx
  value_type = max
[]
[AuxVariables/vry]
[]
[AuxVariables/ary]
[]
[AuxKernels/ary]
  type = AqueductNewmarkAcceleration
  variable = ary
  displacement = ry
  velocity = vry
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vry]
  type = AqueductNewmarkVelocity
  variable = vry
  acceleration = ary
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_ry]
  type = AqueductFlatShellInertia
  variable = ry
  component = 4
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vry_min]
  type = NodalExtremeValue
  variable = vry
  value_type = min
[]
[Postprocessors/vry_max]
  type = NodalExtremeValue
  variable = vry
  value_type = max
[]
[Postprocessors/ary_min]
  type = NodalExtremeValue
  variable = ary
  value_type = min
[]
[Postprocessors/ary_max]
  type = NodalExtremeValue
  variable = ary
  value_type = max
[]
[AuxVariables/vrz]
[]
[AuxVariables/arz]
[]
[AuxKernels/arz]
  type = AqueductNewmarkAcceleration
  variable = arz
  displacement = rz
  velocity = vrz
  beta = .275625
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[AuxKernels/vrz]
  type = AqueductNewmarkVelocity
  variable = vrz
  acceleration = arz
  gamma = .55
  static_end_time = 2
  execute_on = TIMESTEP_END
[]
[Kernels/mass_rz]
  type = AqueductFlatShellInertia
  variable = rz
  component = 5
  density = 7850
  thickness = .1
  offset = .5
  hht_alpha = -.05
  hht_start_time = 2
[]
[Postprocessors/vrz_min]
  type = NodalExtremeValue
  variable = vrz
  value_type = min
[]
[Postprocessors/vrz_max]
  type = NodalExtremeValue
  variable = vrz
  value_type = max
[]
[Postprocessors/arz_min]
  type = NodalExtremeValue
  variable = arz
  value_type = min
[]
[Postprocessors/arz_max]
  type = NodalExtremeValue
  variable = arz
  value_type = max
[]
[AuxVariables/peq]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/peq]
  type = MaterialStdVectorAux
  variable = peq
  property = aq_shell_state
  index = 4
  selected_qp = 0
  execute_on = TIMESTEP_END
[]
[VectorPostprocessors/native]
  type = AqueductElementSampler
  variable = 'peq'
  element_ids = '0'
  sort_by = id
  execute_on = TIMESTEP_END
[]
[Postprocessors/clock]
  type = TimePostprocessor
[]
[UserObjects/assert_rates]
  type = Terminator
  expression = 'abs(vux_min)>1e-8 | abs(vux_max)>1e-8 | abs(aux_min)>1e-7 | abs(aux_max)>1e-7 | abs(vuy_min)>1e-8 | abs(vuy_max)>1e-8 | abs(auy_min)>1e-7 | abs(auy_max)>1e-7 | abs(vuz_min)>1e-8 | abs(vuz_max)>1e-8 | abs(auz_min)>1e-7 | abs(auz_max)>1e-7 | abs(vrx_min)>1e-8 | abs(vrx_max)>1e-8 | abs(arx_min)>1e-7 | abs(arx_max)>1e-7 | abs(vry_min)>1e-8 | abs(vry_max)>1e-8 | abs(ary_min)>1e-7 | abs(ary_max)>1e-7 | abs(vrz_min)>1e-8 | abs(vrz_max)>1e-8 | abs(arz_min)>1e-7 | abs(arz_max)>1e-7'
  error_level = ERROR
  execute_on = TIMESTEP_END
  message = 'Shell handoff generated spurious physical rates'
[]
[UserObjects/assert_loaded]
  type = Terminator
  expression = 'clock>=1 & tip_uz>=-1e-12'
  error_level = ERROR
  execute_on = TIMESTEP_END
  message = 'Shell gravity fixture did not deform'
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 2.15
  dtmin = 1e-10
  dtmax = .2
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[Executioner/TimeIntegrator]
  type = AqueductHHT
  hht_alpha = -.05
  beta = .275625
  gamma = .55
  hht_start_time = 2
  quasistatic_initialization = true
  inactive_tsteps = 0
[]
[Executioner/TimeStepper]
  type = AqueductPhaseDT
  function = clock
  static_end_time = 2
  static_dtmin = 1e-10
  static_dtmax = .2
  dynamic_dtmin = 1e-8
  dynamic_dtmax = .05
[]
[Outputs/history]
  type = AqueductScheduledCSV
  file_base = shell_handoff_history
  mark_interval = .01
  show = 'clock tip_uz tip_ux tip_ry vux_min vux_max aux_min aux_max vuy_min vuy_max auy_min auy_max vuz_min vuz_max auz_min auz_max vrx_min vrx_max arx_min arx_max vry_min vry_max ary_min ary_max vrz_min vrz_max arz_min arz_max'
  phase_origin = 2
  phase_end = 2.15
  execute_on = 'TIMESTEP_END FINAL'
[]
[Outputs/overview]
  type = AqueductScheduledExodus
  file_base = shell_overview
  mark_interval = .1
  show = 'ux uy uz peq'
  phase_origin = 2
  phase_end = 2.15
  execute_on = 'TIMESTEP_END FINAL'
[]
[Outputs/damage]
  type = AqueductScheduledExodus
  file_base = shell_damage
  mark_interval = .5
  show = 'peq'
  phase_origin = 2
  phase_end = 2.15
  execute_on = 'TIMESTEP_END FINAL'
[]
[Outputs/region]
  type = AqueductScheduledExodus
  file_base = shell_region
  mark_interval = .05
  file = shell-output.msh
  show = 'peq'
  phase_origin = 2
  phase_end = 2.15
  execute_on = 'TIMESTEP_END FINAL'
[]
[Outputs/native]
  type = AqueductScheduledCSV
  file_base = shell_native
  mark_interval = .05
  show = 'native'
  time_data = true
  phase_origin = 2
  phase_end = 2.15
  execute_on = 'TIMESTEP_END FINAL'
[]
