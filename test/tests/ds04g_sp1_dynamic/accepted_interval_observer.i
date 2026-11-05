# TEST ONLY: no material; stepper retry/accept plus independent execute/output identities.
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]
[Problem]
  type = CDPObservationProblem
  observation_enabled = true
  observation_start = 0
  observation_end = 0.01
  observation_labels = '1'
  observation_vertices = '0 0 1 0 1 1 0 1'
  observation_prefix = 'accepted-output-observer'
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
  boundary = 'left right bottom top'
  value = 0
[]
[Executioner]
  type = Transient
  end_time = 0.007
  num_steps = 8
  [TimeStepper]
    type = DS04GAcceptedIntervalTestStepper
    dt = 0.002
    check_observation = true
  []
[]
[Outputs]
  csv = true
[]
