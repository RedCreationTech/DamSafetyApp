!include single_hex8_compression.i

[Materials/cdp]
  maximum_tensile_history_increment = 2e-7
[]

[Postprocessors/material_state_dt]
  type = MaterialTimeStepPostprocessor
  maximum_value = 1
[]

[Executioner/TimeStepper]
  type = IterationAdaptiveDT
  dt = 0.05
  optimal_iterations = 8
  iteration_window = 3
  growth_factor = 1.15
  cutback_factor = 0.5
  timestep_limiting_postprocessor = material_state_dt
  reject_large_step = true
  reject_large_step_threshold = 0.999
[]
