# TEST ONLY dynamic assembly; original snapshot U, scaled manufactured histories.
# Test-only, actual frozen-state assembly; no equilibrium solve or time advance.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Mesh]
  type = FileMesh
  file = e1908.msh
[]
[Problem]
  type = DS04GSP1DynamicAuditProblem
  snapshot_u = '-0.17968058665763284 0.00075390720977902852 -0.17980894196250521 0.00077284601911782396 -0.17974079686732117 0.00070225505345020709 -0.17962126863831249 0.00069130034353451924'
  snapshot_v = '-0.22240763980807951 0.00016979798881875579 -0.22517708161297528 0.0040798984681508031 -0.22357586852191419 0.0033801733512463794 -0.22280323691764459 0.0017346168157602671'
  snapshot_a = '4.5395029169119896 -2.8708875920298609 1.0218357571186241 2.2918955344233076 3.1312722292728856 0.93031296497793337 3.9296799425167732 -0.5873125157641147'
  snapshot_stress = '598700.07941166195 1465857.5410785573 433101.55453604308 -910709.51575746667'
[]
[Physics/SolidMechanics/FiniteNewmark/all]
  add_variables = true
  spatial_discretization = sp1_center_hg
  sp1_consistent_tangent = true
  sp1_gravity_jacobian = true
  sp1_two_interval_history = true
  sp1_density = 2400
  sp1_static_end = 1.000001
  sp1_gravity = -9.8
  enable_finite_newmark_damping = false
  expected_youngs_modulus = 2.97915e10
  expected_poissons_ratio = 0.2
  expected_stiffness_damping = 0.003604
  velocities = 'vel_x vel_y'
  accelerations = 'accel_x accel_y'
  newmark_beta = 0.275625
  newmark_gamma = 0.55
  hht_alpha = -0.05
  mass_damping_coefficient = 0.05
  stiffness_damping_coefficient = 0.003604
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
[]

[Kernels]
  [static_x]
    type = DS04GSP1StressDivergence
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [static_y]
    type = DS04GSP1StressDivergence
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [gravity_y]
    type = DS04GSP1Gravity
    variable = disp_y
    value = -9.8
    alpha = -0.05
  []
[]

[Materials]
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2400
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.97915e10
    poissons_ratio = 0.2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = cdp_stress_update
    perform_finite_strain_rotations = false
  []
  [cdp_stress_update]
    type = AbaqusCDPStressUpdate
    maximum_substeps = 256
    maximum_strain_increment = 2.5e-5
    youngs_modulus = 2.97915e10
    poissons_ratio = 0.2
    dilation_angle = 30
    eccentricity = 0.1
    biaxial_to_uniaxial_compression_ratio = 1.16
    tensile_meridian_ratio = 0.6667
    viscosity = 0.0005
    tension_recovery = 0
    compression_recovery = 1
    compression_hardening_file = compression_hardening.csv
    compression_damage_file = compression_damage.csv
    tension_stiffening_file = tension_stiffening.csv
    tension_damage_file = tension_damage.csv
  []
[]


[Materials/probe]
  type = CDPAssemblyProbe
  probe_elements = 0
  capture_kinematics = true
[]
[Controls/phase]
  type = TimePeriod
  execute_on = 'INITIAL TIMESTEP_BEGIN'
  start_time = 0
  end_time = 1.000001
  enable_objects = '*/static_x */static_y'
  disable_objects = '*/TM_all0 */TM_all1 */TM_all_inertia_x */TM_all_inertia_y */TM_all_accel_x */TM_all_accel_y */TM_all_vel_x */TM_all_vel_y */gravity_y */ds04g_sp1_gravity_1'
[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 0.002
  dt = 0.002
  num_steps = 0
[]
[Outputs]
  csv = true
  file_base = audit_out
[]
