# SP1 admission and build preparation

Default `spatial_discretization=legacy_q4` delegates to the897 path with no new kernels/quadrature changes. `sp1_center_hg` is an opt-in prototype. Add `sp1_density=2400` with the retained expected E/nu/stiffness-damping fields, PLANE_STRAIN/FINITE and `enable_finite_newmark_damping=false`. It owns volume quadrature; no separate Quadrature override is admitted. Single subdomain, FIRST LAGRANGE/QUAD4 only.

It creates reference-assembly HG/mass corrections and current-geometry uniform Y-gravity correction. `sp1_static_end=1.000001` MUST match the original TimePeriod phase, `sp1_gravity=-9.8` the original uniform Gravity. Those configuration invariants must be verified against the full frozen input before a LIMS trial. Do not convert the absolute-coordinate baseline to relative ground loading. No4QP restart/projection.

The unit tests check geometry/mode/AD suboperators and replay six frozen material-core inputs. They do NOT replace actual FE element assembly, FINITE prestress rotation, full coupled Jacobian, seven frozen LIMS cases or the same-candidate DEFAULT3.50s control. No new finite-element run is authorized by this file. Never promote this candidate or register a production solver solely from unit tests or --check-input.

No production constitutive formula/table is modified. Abaqus native matrices are initial elasticity only; dynamic/CDP equivalence remains open.
