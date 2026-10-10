# TASK-DAM2B-CDP-RATE-053: real local Newton observation only.
!include cdp-new-5ms.i
[Problem]
  type = CDPTrialProblem
  trial_start = 0.045
  trial_end = 0.085
  trial_capture_times = '0.045 0.085'
  trial_prefix = cdp_jac
  capture_dof_map = true
  capture_actual_jacobian = true
  capture_material_jacobian = true
  capture_material_substeps = true
  capture_material_local_newton = true
[]
[Materials]
  [actual_material_observer]
    type = CDPAssemblyProbe
    probe_elements = '0'
  []
[]
[Executioner]
  petsc_options = '-snes_test_jacobian'
[]
