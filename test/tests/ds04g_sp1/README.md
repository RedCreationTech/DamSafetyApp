# SP1 admission and build preparation

Default `spatial_discretization=legacy_q4` delegates to the897 path with no new kernels/quadrature changes. `sp1_center_hg` is an opt-in prototype. Add `sp1_density=2400` with the retained expected E/nu/stiffness-damping fields, PLANE_STRAIN/FINITE and `enable_finite_newmark_damping=false`. It owns volume quadrature; no separate Quadrature override is admitted. Single subdomain, FIRST LAGRANGE/QUAD4 only.

It creates reference-assembly HG/mass corrections and current-geometry uniform Y-gravity correction. `sp1_static_end=1.000001` MUST match the original TimePeriod phase, `sp1_gravity=-9.8` the original uniform Gravity. Those configuration invariants must be verified against the full frozen input before a LIMS trial. Do not convert the absolute-coordinate baseline to relative ground loading. No4QP restart/projection.

The unit tests check geometry/mode/AD suboperators and replay six frozen material-core inputs. They do NOT replace actual FE element assembly, FINITE prestress rotation, full coupled Jacobian, seven frozen LIMS cases or the same-candidate DEFAULT3.50s control. No new finite-element run is authorized by this file. Never promote this candidate or register a production solver solely from unit tests or --check-input.

No production constitutive formula/table is modified. Abaqus native matrices are initial elasticity only; dynamic/CDP equivalence remains open.

## Opt-in rotation-false tangent (work package A)

`sp1_consistent_tangent=false` is the default. True is admitted only with
`spatial_discretization=sp1_center_hg`, OFF and the frozen TaylorExpansion/false
stress path. It selects `DS04GSP1StressDivergence` for the generated dynamic
kernel. Explicit static kernels must use that type with `sp1_dynamic=false`
(its default). Both residuals directly delegate to their original parent.
The Jacobian contracts the unchanged material D_alg with the exact Taylor
increment derivative and adds current volume/gradient terms. Dynamic geometry
uses the complete A*stress+B*old+C*older stress, not just A*current stress.
The original reference HG and its own history weights remain separate.

The frozen audit can request `Problem/audit_kinematics=true`,
`Problem/refine_compression=true`, and `Materials/probe/capture_kinematics=true`.
It captures same-call Jacobian material C plus nine-entry F/Fold/rotation,
and a field explicitly named `T_reconstructed` (not a directly observed CDP
entrance). Compression adds h3=3e-8 and h4=1e-8 m; old h0/h1/h2 stay unchanged.
These flags do not execute time steps; actually running this FE assembly
fixture still requires LIMS and belongs to the subsequent combined gate.

`Problem/refine_nonlinear=true` also adds h3/h4 to affine and tension. This is
needed to show two small-h scales after large-h truncation/partition effects;
it does not alter the base U, material parameters, or the original h0/h1/h2.
The consistent mode rejects `tangent_operator=elastic` rather than silently
treating C0 as the full material update derivative.

The added C++ tests independently finite-difference Taylor kinematics and the
material/geometry force chain for six frozen real material inputs in both
test geometries, including finite synthetic old geometry and prescribed old
stress weights. They are not an actual FE/Newmark/rollback test or a new dam
trajectory. Same-binary seven frozen LIMS regressions remain required.

## SP1 gravity Jacobian (default off)

`sp1_gravity_jacobian=false` retains the original Gravity/default LEGACY behavior.
To opt in, set it true in `Physics/SolidMechanics/FiniteNewmark/all` and change
only the existing uniform Y-gravity block type to `DS04GSP1Gravity`. The wrapper
inherits the original `Gravity` residual and adds its current-volume X/Y
derivatives. It uses `GlobalParams/displacements`; the original block name,
controls, value, function, alpha, and density remain. The SP1 AD gravity
correction is unchanged. True in LEGACY or true without the paired wrapper is
rejected. No constitutive, history-weight, or upstream default is changed.

`../ds04g_sp1_dynamic/e*_gravity_jacobian.i` preserve the previous three meshes
and four prescribed histories. Formal runs go through LIMS; these are zero-step
manufactured-state assembly audits, not full-dam response or CDP regressions.
