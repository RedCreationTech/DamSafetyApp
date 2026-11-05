# TEST ONLY: one QUAD4, real CDP, one injected reject, and passive accepted/output recording.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  elem_type = QUAD4
[]
[Problem]
  type = CDPObservationProblem
  observation_enabled = true
  trace_state_swap = true
  observation_start = 0
  observation_end = 0.01
  observation_labels = '1'
  observation_vertices = '0 0 1 0 1 1 0 1'
  observation_prefix = 'accepted-output-observer'
[]
[Physics/SolidMechanics/QuasiStatic/all]
  strain = SMALL
  incremental = true
  add_variables = true
  planar_formulation = PLANE_STRAIN
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
[]
[Functions/load_path]
  type = PiecewiseLinear
  x = '0 0.007'
  y = '0 2e-4'
[]
[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [load_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = load_path
  []
[]
[AuxVariables]
  [damage_t]
    order = CONSTANT
    family = MONOMIAL
  []
  [accepted_substeps]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [damage_t]
    type = MaterialRealAux
    variable = damage_t
    property = DamageT
    execute_on = 'initial timestep_end'
  []
  [accepted_substeps]
    type = MaterialRealAux
    variable = accepted_substeps
    property = cdp_accepted_substeps
    execute_on = 'initial timestep_end'
  []
[]
[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 3.04e10
    poissons_ratio = 0.2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = cdp
    perform_finite_strain_rotations = false
  []
  [cdp]
    type = AbaqusCDPStressUpdate
    compression_hardening_file = ../cdp_material_table/data/compression_hardening.csv
    compression_damage_file = ../cdp_material_table/data/compression_damage.csv
    tension_stiffening_file = ../cdp_material_table/data/tension_stiffening.csv
    tension_damage_file = ../cdp_material_table/data/tension_damage.csv
    youngs_modulus = 3.04e10
    poissons_ratio = 0.2
    dilation_angle = 36.31
    eccentricity = 0.1
    biaxial_to_uniaxial_compression_ratio = 1.16
    tensile_meridian_ratio = 0.667
    viscosity = 0.0005
    maximum_strain_increment = 2.5e-5
  []
[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 0.007
  num_steps = 8
  nl_max_its = 40
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  automatic_scaling = true
  [TimeStepper]
    type = DS04GAcceptedIntervalTestStepper
    dt = 0.002
    check_observation = true
    check_material_observation = true
  []
[]
[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = cdp_accepted_output_observer_out
[]
