# Combined scalar damage viscosity candidate

Status: opt-in diagnostic candidate. The default remains `false`.

The published Abaqus CDP Duvaut-Lions regularization advances one scalar
stiffness degradation state,

```text
dot(d_v) = (d - d_v) / mu
sigma = (1 - d_v) D0 : (epsilon - epsilon_vpl)
```

The previous prototype independently relaxed `DamageT` and `DamageC`, then
combined the relaxed branches with the current stress weight. Branch relaxation
and nonlinear damage combination do not commute for mixed principal stress and
stiffness recovery.

`use_scalar_damage_viscosity = true` keeps the existing branch histories for
`DamageT` and `DamageC` output, but advances a separate backward-Euler
`cdp_combined_damage` state from the inviscid combined damage. This scalar alone
sets the stress stiffness factor and its algorithmic tangent. Plastic strain,
hardening variables, material tables, yield surface, flow potential, recovery
parameters, local tolerances, and substepping are unchanged.

The candidate must pass unit tangent checks, default-off regressions, opt-in
tension/compression/mixed TestHarness cases, and paired LIMS tension/compression
validation before acceptance.
