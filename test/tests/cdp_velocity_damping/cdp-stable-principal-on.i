# TASK-DAM2B-CDP-2D-055: no tolerance, partition, damping or loading change.
!include cdp-local-newton-trace.i
[Materials]
  [cdp]
    use_stable_principal_stress = true
  []
[]
