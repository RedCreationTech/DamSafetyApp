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
  [jacobian_fallbacks]
    order = CONSTANT
    family = MONOMIAL
  []
  [failed_material_calls]
    order = CONSTANT
    family = MONOMIAL
  []
  [attempted_partitions]
    order = CONSTANT
    family = MONOMIAL
  []
  [maximum_partition_depth]
    order = CONSTANT
    family = MONOMIAL
  []
  [local_factorizations]
    order = CONSTANT
    family = MONOMIAL
  []
  [local_backsolves]
    order = CONSTANT
    family = MONOMIAL
  []
  [integration_microseconds]
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
  [jacobian_fallbacks]
    type = MaterialRealAux
    variable = jacobian_fallbacks
    property = cdp_jacobian_fallbacks
    execute_on = 'initial timestep_end'
  []
  [failed_material_calls]
    type = MaterialRealAux
    variable = failed_material_calls
    property = cdp_failed_material_calls
    execute_on = 'initial timestep_end'
  []
  [attempted_partitions]
    type = MaterialRealAux
    variable = attempted_partitions
    property = cdp_attempted_partitions
    execute_on = 'initial timestep_end'
  []
  [maximum_partition_depth]
    type = MaterialRealAux
    variable = maximum_partition_depth
    property = cdp_maximum_partition_depth
    execute_on = 'initial timestep_end'
  []
  [local_factorizations]
    type = MaterialRealAux
    variable = local_factorizations
    property = cdp_local_factorizations
    execute_on = 'initial timestep_end'
  []
  [local_backsolves]
    type = MaterialRealAux
    variable = local_backsolves
    property = cdp_local_backsolves
    execute_on = 'initial timestep_end'
  []
  [integration_microseconds]
    type = MaterialRealAux
    variable = integration_microseconds
    property = cdp_integration_microseconds
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
    enable_performance_diagnostics = true
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
  [maximum_jacobian_fallbacks]
    type = ElementExtremeValue
    variable = jacobian_fallbacks
    value_type = max
  []
  [maximum_failed_material_calls]
    type = ElementExtremeValue
    variable = failed_material_calls
    value_type = max
  []
  [maximum_attempted_partitions]
    type = ElementExtremeValue
    variable = attempted_partitions
    value_type = max
  []
  [maximum_partition_depth]
    type = ElementExtremeValue
    variable = maximum_partition_depth
    value_type = max
  []
  [maximum_local_factorizations]
    type = ElementExtremeValue
    variable = local_factorizations
    value_type = max
  []
  [maximum_local_backsolves]
    type = ElementExtremeValue
    variable = local_backsolves
    value_type = max
  []
  [maximum_integration_microseconds]
    type = ElementExtremeValue
    variable = integration_microseconds
    value_type = max
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
