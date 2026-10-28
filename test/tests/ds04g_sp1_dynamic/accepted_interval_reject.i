# TEST ONLY: actual Transient dt bookkeeping, no dam/material model.
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]
[Variables/u]
[]
[Kernels/time]
  type = TimeDerivative
  variable = u
[]
[Kernels/reaction]
  type = Reaction
  variable = u
[]
[BCs/zero]
  type = DirichletBC
  variable = u
  boundary = 'left right'
  value = 0
[]
[Executioner]
  type = Transient
  end_time = 0.007
  num_steps = 8
  [TimeStepper]
    type = DS04GAcceptedIntervalTestStepper
    dt = 0.002
  []
[]
[Outputs]
  csv = true
[]
