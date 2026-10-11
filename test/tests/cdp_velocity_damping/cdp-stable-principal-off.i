# TASK-DAM2B-CDP-2D-055: same damping, only the principal-spectrum switch differs.
!include cdp-local-newton-trace.i
[Materials]
  [cdp]
    use_stable_principal_stress = false
  []
[]
