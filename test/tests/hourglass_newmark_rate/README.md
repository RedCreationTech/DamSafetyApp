# Independent HG damping time rate (development only)

`damping_rate` defaults to `legacy_secant`. Opt-in `newmark_velocity` needs
nodal `velocity` and `acceleration` histories and matching Newmark beta/gamma.
Elastic stiffness, spatial C, material Rayleigh coefficients and HHT are unchanged.

These bounded fixtures test actual kernel AD and a forced rejection. CSV clean/retry
and old binary/default/explicit legacy comparisons are required by the integration
preparation. They do not prove full CDP state rollback or dam response accuracy.
