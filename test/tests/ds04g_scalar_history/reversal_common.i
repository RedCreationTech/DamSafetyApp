# Test-only one QUAD4 material history. No gravity, earthquake or dam solve.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmin = 0
  xmax = 0.001
  ymin = 0
  ymax = 0.001
  elem_type = QUAD4
[]
[Executioner/Quadrature]
  type = GAUSS
  order = CONSTANT
[]
[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  planar_formulation = PLANE_STRAIN
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy'
[]
[BCs]
  [left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = path_x
  []
  [top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = path_y
  []
[]
[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.97915e10
    poissons_ratio = .2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = cdp
    perform_finite_strain_rotations = false
  []
  [cdp]
    type = AbaqusCDPStressUpdate
    youngs_modulus = 2.97915e10
    poissons_ratio = .2
    dilation_angle = 30
    eccentricity = .1
    biaxial_to_uniaxial_compression_ratio = 1.16
    tensile_meridian_ratio = .6667
    tension_recovery = 0
    compression_recovery = 1
    viscosity = .0005
    maximum_substeps = 256
    maximum_strain_increment = 2.5e-5
    compression_hardening_file = compression_hardening.csv
    compression_damage_file = compression_damage.csv
    tension_stiffening_file = tension_stiffening.csv
    tension_damage_file = tension_damage.csv
  []
[]
[AuxVariables]
  [damage_t]
    family = MONOMIAL
    order = CONSTANT
  []
  [damage_c]
    family = MONOMIAL
    order = CONSTANT
  []
  [combined]
    family = MONOMIAL
    order = CONSTANT
  []
  [kappa_t]
    family = MONOMIAL
    order = CONSTANT
  []
  [kappa_c]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_0]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_1]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_2]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_3]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_4]
    family = MONOMIAL
    order = CONSTANT
  []
  [backbone_5]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_0]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_1]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_2]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_3]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_4]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscous_5]
    family = MONOMIAL
    order = CONSTANT
  []
[]
[AuxKernels]
  [damage_t]
    type = MaterialRealAux
    variable = damage_t
    property = DamageT
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damage_c]
    type = MaterialRealAux
    variable = damage_c
    property = DamageC
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [combined]
    type = MaterialRealAux
    variable = combined
    property = cdp_combined_damage
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_t]
    type = MaterialRealAux
    variable = kappa_t
    property = cdp_kappa_t
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_c]
    type = MaterialRealAux
    variable = kappa_c
    property = cdp_kappa_c
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_0]
    type = RankTwoAux
    variable = backbone_0
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 0
    index_j = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_1]
    type = RankTwoAux
    variable = backbone_1
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 1
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_2]
    type = RankTwoAux
    variable = backbone_2
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 2
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_3]
    type = RankTwoAux
    variable = backbone_3
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 0
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_4]
    type = RankTwoAux
    variable = backbone_4
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 1
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_5]
    type = RankTwoAux
    variable = backbone_5
    rank_two_tensor = cdp_backbone_plastic_strain
    index_i = 0
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_0]
    type = RankTwoAux
    variable = viscous_0
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 0
    index_j = 0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_1]
    type = RankTwoAux
    variable = viscous_1
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 1
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_2]
    type = RankTwoAux
    variable = viscous_2
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 2
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_3]
    type = RankTwoAux
    variable = viscous_3
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 0
    index_j = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_4]
    type = RankTwoAux
    variable = viscous_4
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 1
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_5]
    type = RankTwoAux
    variable = viscous_5
    rank_two_tensor = cdp_viscous_plastic_strain
    index_i = 0
    index_j = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Postprocessors]
  [damage_t]
    type = ElementAverageValue
    variable = damage_t
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damage_c]
    type = ElementAverageValue
    variable = damage_c
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [combined]
    type = ElementAverageValue
    variable = combined
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_t]
    type = ElementAverageValue
    variable = kappa_t
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_c]
    type = ElementAverageValue
    variable = kappa_c
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_0]
    type = ElementAverageValue
    variable = backbone_0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_1]
    type = ElementAverageValue
    variable = backbone_1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_2]
    type = ElementAverageValue
    variable = backbone_2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_3]
    type = ElementAverageValue
    variable = backbone_3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_4]
    type = ElementAverageValue
    variable = backbone_4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [backbone_5]
    type = ElementAverageValue
    variable = backbone_5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_0]
    type = ElementAverageValue
    variable = viscous_0
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_1]
    type = ElementAverageValue
    variable = viscous_1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_2]
    type = ElementAverageValue
    variable = viscous_2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_3]
    type = ElementAverageValue
    variable = viscous_3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_4]
    type = ElementAverageValue
    variable = viscous_4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [viscous_5]
    type = ElementAverageValue
    variable = viscous_5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xx]
    type = ElementAverageValue
    variable = strain_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_yy]
    type = ElementAverageValue
    variable = strain_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_zz]
    type = ElementAverageValue
    variable = strain_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [strain_xy]
    type = ElementAverageValue
    variable = strain_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Outputs]
  file_base = reversal_out
  csv = true
  exodus = true
  [cp]
    type = Checkpoint
    num_files = 2
  []
[]
