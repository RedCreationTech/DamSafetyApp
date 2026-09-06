# Captured local failure regression data

Material tables from LIMS Jobs job_20260906_090429_k8cuhe (D01) and job_20260906_090600_7vltup (D02), generated from the existing public-validation Abaqus-compatible CDP cases. E=29791500000 Pa, nu=0.2, dilation=36 degrees, eccentricity=0.1, fb0/fc0=1.16, Kc=0.667.

Five trial strains and prior local-substep states are embedded in AbaqusCDPLocalIntegratorTest.C. These are failed Newton trial inputs, not accepted global material histories. The tests exercise the local backbone equations, not finite-element equilibrium.
