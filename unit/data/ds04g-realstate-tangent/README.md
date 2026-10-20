# DS04G recorded material-point fixtures

Source: completed LIMS Job `job_20261020_223540_ds33p4`, solver
`8a0114cff87eac51b5db1a796434906301388d0e`. These are recorded inputs, not
synthetic expected physics or an Abaqus local-integrator reference.

`source-identity.json` records the seven attempt/evaluation/call/element/QP keys
and hashes of the numeric fixtures and the unchanged material tables. QP is zero
based. F1/F2 have different actual dt (1ms/2ms), not just different final N.

Each whitespace fixture contains id, old strain[6], new strain[6], old state[16],
dt and N, recorded Cauchy stress[6], new state[16], column-major tangent[36], and
a local-snapshot flag. Local snapshots additionally contain target[6], old
backbone[8], stress/strain scales, z[9], residual[9], row-major AD[81]/FD[81].
State order is backbone plastic[6], kappa_t/c, viscous plastic[6], DT/DC. Tensor
order is xx/yy/zz/xy/yz/xz, with tensor shear; stress is in Pa.

Generate from the integration repository's frozen selection JSON with
`tools/diagnostics/prepare_ds04g_tangent_fixture.py --selected <json> --input
<original-input-directory> --output unit/data/ds04g-realstate-tangent`.

Run from the application root, in the locked MOOSE environment:

```
mkdir -p realstate-results
DS04G_TANGENT_OUTPUT="$PWD/realstate-results" ./unit/DamSafetyApp-unit-opt --gtest_filter='DS04G/AbaqusCDPRealState.*'
```

The seven tests assert frozen input replay, repeated evaluation, immutable old
state and recorded local residual/matrices. Derivative differences are diagnostic
data, not assertions that every real state is differentiable. A green unit test
does not assert that G1/G2/G3 pass. The integration report evaluates stable
step-size windows and branch/partition changes without altering any tolerances.

Output CSV columns identify layer/mode/direction/h/side/N/signature/status.
`sample`/`substate` values are stress[6] then state[16]; side +/-2 supplies the
second point for a second-order one-sided difference. Both sides are retained,
including exceptions. `path` holds partitions/cutbacks/iterations/AD-FD fallbacks.
An adaptive sample uses the original primal integrator, then reproduces its
selected fixed-N path to expose substep states; the two returned results must be
bitwise equal. Failed N attempts in baseline linearized evaluation are retained
by the existing passive observation sink. Perturbed primal results expose final
N and aggregate cutbacks, not a complete per-attempt failure trace.

`local-*` rows expose recorded trial residual derivatives separately from roots
captured under `trace-*`. The only production header change is a friend access
declaration. No production formula, default, selector or runtime entry is added.
