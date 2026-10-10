[Materials]
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
