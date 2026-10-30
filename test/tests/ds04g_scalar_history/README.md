# OPT-04F bounded material history fixtures

One QUAD4/center quadrature/FINITE plane strain/rotations false; prescribed affine path, no gravity/EQ. 17 adapter histories are exported. OFF omits new flag for a9 comparison; ON adds only flag. Part ends at .030, same-mode restart starts .030 from actual checkpoint name supplied after listing archive. Never enable ON using an old a9 damaged checkpoint.

Reject injects one global failure after a converged trial, using existing serial .01s/8step-bounded tester. Clean time sequence replays its accepted grid; differences in physical histories must be zero. Source defaults/LEGACY, backbone, mu and branch outputs unchanged. These are not the seven frozen CDP regression cases. All real FE checks through LIMS.
