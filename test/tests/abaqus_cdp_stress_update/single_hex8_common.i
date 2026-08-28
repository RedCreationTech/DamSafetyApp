[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX8
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = SMALL
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy
                     strain_xx strain_yy strain_zz strain_xy vonmises_stress'
[]

[AuxVariables]
  [damage_t]
    order = CONSTANT
    family = MONOMIAL
  []
  [damage_c]
    order = CONSTANT
    family = MONOMIAL
  []
  [combined_damage]
    order = CONSTANT
    family = MONOMIAL
  []
  [stiffness_factor]
    order = CONSTANT
    family = MONOMIAL
  []
  [kappa_t]
    order = CONSTANT
    family = MONOMIAL
  []
  [kappa_c]
    order = CONSTANT
    family = MONOMIAL
  []
  [local_iterations]
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
  [damage_c]
    type = MaterialRealAux
    variable = damage_c
    property = DamageC
    execute_on = 'initial timestep_end'
  []
  [combined_damage]
    type = MaterialRealAux
    variable = combined_damage
    property = cdp_combined_damage
    execute_on = 'initial timestep_end'
  []
  [stiffness_factor]
    type = MaterialRealAux
    variable = stiffness_factor
    property = cdp_stiffness_factor
    execute_on = 'initial timestep_end'
  []
  [kappa_t]
    type = MaterialRealAux
    variable = kappa_t
    property = cdp_kappa_t
    execute_on = 'initial timestep_end'
  []
  [kappa_c]
    type = MaterialRealAux
    variable = kappa_c
    property = cdp_kappa_c
    execute_on = 'initial timestep_end'
  []
  [local_iterations]
    type = MaterialRealAux
    variable = local_iterations
    property = cdp_local_iterations
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

[Postprocessors]
  [average_stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [average_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [average_stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [average_strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  []
  [average_strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  []
  [average_strain_xy]
    type = ElementAverageValue
    variable = strain_xy
  []
  [average_damage_t]
    type = ElementAverageValue
    variable = damage_t
  []
  [average_damage_c]
    type = ElementAverageValue
    variable = damage_c
  []
  [average_stiffness_factor]
    type = ElementAverageValue
    variable = stiffness_factor
  []
  [average_kappa_t]
    type = ElementAverageValue
    variable = kappa_t
  []
  [average_kappa_c]
    type = ElementAverageValue
    variable = kappa_c
  []
  [average_local_iterations]
    type = ElementAverageValue
    variable = local_iterations
  []
  [average_accepted_substeps]
    type = ElementAverageValue
    variable = accepted_substeps
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1
  dt = 0.05
  dtmin = 1e-5
  nl_max_its = 40
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  automatic_scaling = true
[]
