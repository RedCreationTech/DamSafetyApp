# DS04G S4: replay S3 gravity 0..1s, then uninterrupted EQ 0..1.5s.
# Same binary/material/mesh; no Exodus-only initialization or history reset.
# Dynamic kernels exist from INITIAL to allocate older stress, but are disabled through t=1.
# Abaqus MSG confirms HHT alpha=-0.05, beta=0.275625, gamma=0.55.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = FileMesh
  file = DamSeismic_04g_si.e
[]

[Physics/SolidMechanics/Dynamic/all]
  add_variables = true
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
    type = StressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [static_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [gravity_y]
    type = Gravity
    variable = disp_y
    value = -9.8
    alpha = -0.05
  []
[]

[BCs]
  [base_acceleration_x]
    type = PresetAcceleration
    variable = disp_x
    boundary = 'DAM_BASE__Dam_1'
    function = base_acceleration
    velocity = vel_x
    acceleration = accel_x
    beta = 0.275625
  []
  [base_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'DAM_BASE__Dam_1'
    value = 0
  []
  [base_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'DAM_BASE__Dam_1'
    value = 0
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

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 2.5
  dt = 0.05
  dtmin = 1e-7
  [TimeStepper]
    type = TimeSequenceStepper
    use_last_t_for_end_time = true
    time_sequence = '0.000000 0.050000 0.100000 0.150000 0.200000 0.250000 0.300000 0.350000 0.400000 0.450000 0.500000 0.550000 0.600000 0.650000 0.700000 0.750000 0.800000 0.850000 0.900000 0.950000 1.000000 1.002000 1.004000 1.006000 1.008000 1.010000 1.012000 1.014000 1.016000 1.018000 1.020000 1.022000 1.024000 1.026000 1.028000 1.030000 1.032000 1.034000 1.036000 1.038000 1.040000 1.042000 1.044000 1.046000 1.048000 1.050000 1.052000 1.054000 1.056000 1.058000 1.060000 1.062000 1.064000 1.066000 1.068000 1.070000 1.072000 1.074000 1.076000 1.078000 1.080000 1.082000 1.084000 1.086000 1.088000 1.090000 1.092000 1.094000 1.096000 1.098000 1.100000 1.102000 1.104000 1.106000 1.108000 1.110000 1.112000 1.114000 1.116000 1.118000 1.120000 1.122000 1.124000 1.126000 1.128000 1.130000 1.132000 1.134000 1.136000 1.138000 1.140000 1.142000 1.144000 1.146000 1.148000 1.150000 1.152000 1.154000 1.156000 1.158000 1.160000 1.162000 1.164000 1.166000 1.168000 1.170000 1.172000 1.174000 1.176000 1.178000 1.180000 1.182000 1.184000 1.186000 1.188000 1.190000 1.192000 1.194000 1.196000 1.198000 1.200000 1.202000 1.204000 1.206000 1.208000 1.210000 1.212000 1.214000 1.216000 1.218000 1.220000 1.222000 1.224000 1.226000 1.228000 1.230000 1.232000 1.234000 1.236000 1.238000 1.240000 1.242000 1.244000 1.246000 1.248000 1.250000 1.252000 1.254000 1.256000 1.258000 1.260000 1.262000 1.264000 1.266000 1.268000 1.270000 1.272000 1.274000 1.276000 1.278000 1.280000 1.282000 1.284000 1.286000 1.288000 1.290000 1.292000 1.294000 1.296000 1.298000 1.300000 1.302000 1.304000 1.306000 1.308000 1.310000 1.312000 1.314000 1.316000 1.318000 1.320000 1.322000 1.324000 1.326000 1.328000 1.330000 1.332000 1.334000 1.336000 1.338000 1.340000 1.342000 1.344000 1.346000 1.348000 1.350000 1.352000 1.354000 1.356000 1.358000 1.360000 1.362000 1.364000 1.366000 1.368000 1.370000 1.372000 1.374000 1.376000 1.378000 1.380000 1.382000 1.384000 1.386000 1.388000 1.390000 1.392000 1.394000 1.396000 1.398000 1.400000 1.402000 1.404000 1.406000 1.408000 1.410000 1.412000 1.414000 1.416000 1.418000 1.420000 1.422000 1.424000 1.426000 1.428000 1.430000 1.432000 1.434000 1.436000 1.438000 1.440000 1.442000 1.444000 1.446000 1.448000 1.450000 1.452000 1.454000 1.456000 1.458000 1.460000 1.462000 1.464000 1.466000 1.468000 1.470000 1.472000 1.474000 1.476000 1.478000 1.480000 1.482000 1.484000 1.486000 1.488000 1.490000 1.492000 1.494000 1.496000 1.498000 1.500000 1.502000 1.504000 1.506000 1.508000 1.510000 1.512000 1.514000 1.516000 1.518000 1.520000 1.522000 1.524000 1.526000 1.528000 1.530000 1.532000 1.534000 1.536000 1.538000 1.540000 1.542000 1.544000 1.546000 1.548000 1.550000 1.552000 1.554000 1.556000 1.558000 1.560000 1.562000 1.564000 1.566000 1.568000 1.570000 1.572000 1.574000 1.576000 1.578000 1.580000 1.582000 1.584000 1.586000 1.588000 1.590000 1.592000 1.594000 1.596000 1.598000 1.600000 1.602000 1.604000 1.606000 1.608000 1.610000 1.612000 1.614000 1.616000 1.618000 1.620000 1.622000 1.624000 1.626000 1.628000 1.630000 1.632000 1.634000 1.636000 1.638000 1.640000 1.642000 1.644000 1.646000 1.648000 1.650000 1.652000 1.654000 1.656000 1.658000 1.660000 1.662000 1.664000 1.666000 1.668000 1.670000 1.672000 1.674000 1.676000 1.678000 1.680000 1.682000 1.684000 1.686000 1.688000 1.690000 1.692000 1.694000 1.696000 1.698000 1.700000 1.702000 1.704000 1.706000 1.708000 1.710000 1.712000 1.714000 1.716000 1.718000 1.720000 1.722000 1.724000 1.726000 1.728000 1.730000 1.732000 1.734000 1.736000 1.738000 1.740000 1.742000 1.744000 1.746000 1.748000 1.750000 1.752000 1.754000 1.756000 1.758000 1.760000 1.762000 1.764000 1.766000 1.768000 1.770000 1.772000 1.774000 1.776000 1.778000 1.780000 1.782000 1.784000 1.786000 1.788000 1.790000 1.792000 1.794000 1.796000 1.798000 1.800000 1.802000 1.804000 1.806000 1.808000 1.810000 1.812000 1.814000 1.816000 1.818000 1.820000 1.822000 1.824000 1.826000 1.828000 1.830000 1.832000 1.834000 1.836000 1.838000 1.840000 1.842000 1.844000 1.846000 1.848000 1.850000 1.852000 1.854000 1.856000 1.858000 1.860000 1.862000 1.864000 1.866000 1.868000 1.870000 1.872000 1.874000 1.876000 1.878000 1.880000 1.882000 1.884000 1.886000 1.888000 1.890000 1.892000 1.894000 1.896000 1.898000 1.900000 1.902000 1.904000 1.906000 1.908000 1.910000 1.912000 1.914000 1.916000 1.918000 1.920000 1.922000 1.924000 1.926000 1.928000 1.930000 1.932000 1.934000 1.936000 1.938000 1.940000 1.942000 1.944000 1.946000 1.948000 1.950000 1.952000 1.954000 1.956000 1.958000 1.960000 1.962000 1.964000 1.966000 1.968000 1.970000 1.972000 1.974000 1.976000 1.978000 1.980000 1.982000 1.984000 1.986000 1.988000 1.990000 1.992000 1.994000 1.996000 1.998000 2.000000 2.002000 2.004000 2.006000 2.008000 2.010000 2.012000 2.014000 2.016000 2.018000 2.020000 2.022000 2.024000 2.026000 2.028000 2.030000 2.032000 2.034000 2.036000 2.038000 2.040000 2.042000 2.044000 2.046000 2.048000 2.050000 2.052000 2.054000 2.056000 2.058000 2.060000 2.062000 2.064000 2.066000 2.068000 2.070000 2.072000 2.074000 2.076000 2.078000 2.080000 2.082000 2.084000 2.086000 2.088000 2.090000 2.092000 2.094000 2.096000 2.098000 2.100000 2.102000 2.104000 2.106000 2.108000 2.110000 2.112000 2.114000 2.116000 2.118000 2.120000 2.122000 2.124000 2.126000 2.128000 2.130000 2.132000 2.134000 2.136000 2.138000 2.140000 2.142000 2.144000 2.146000 2.148000 2.150000 2.152000 2.154000 2.156000 2.158000 2.160000 2.162000 2.164000 2.166000 2.168000 2.170000 2.172000 2.174000 2.176000 2.178000 2.180000 2.182000 2.184000 2.186000 2.188000 2.190000 2.192000 2.194000 2.196000 2.198000 2.200000 2.202000 2.204000 2.206000 2.208000 2.210000 2.212000 2.214000 2.216000 2.218000 2.220000 2.222000 2.224000 2.226000 2.228000 2.230000 2.232000 2.234000 2.236000 2.238000 2.240000 2.242000 2.244000 2.246000 2.248000 2.250000 2.252000 2.254000 2.256000 2.258000 2.260000 2.262000 2.264000 2.266000 2.268000 2.270000 2.272000 2.274000 2.276000 2.278000 2.280000 2.282000 2.284000 2.286000 2.288000 2.290000 2.292000 2.294000 2.296000 2.298000 2.300000 2.302000 2.304000 2.306000 2.308000 2.310000 2.312000 2.314000 2.316000 2.318000 2.320000 2.322000 2.324000 2.326000 2.328000 2.330000 2.332000 2.334000 2.336000 2.338000 2.340000 2.342000 2.344000 2.346000 2.348000 2.350000 2.352000 2.354000 2.356000 2.358000 2.360000 2.362000 2.364000 2.366000 2.368000 2.370000 2.372000 2.374000 2.376000 2.378000 2.380000 2.382000 2.384000 2.386000 2.388000 2.390000 2.392000 2.394000 2.396000 2.398000 2.400000 2.402000 2.404000 2.406000 2.408000 2.410000 2.412000 2.414000 2.416000 2.418000 2.420000 2.422000 2.424000 2.426000 2.428000 2.430000 2.432000 2.434000 2.436000 2.438000 2.440000 2.442000 2.444000 2.446000 2.448000 2.450000 2.452000 2.454000 2.456000 2.458000 2.460000 2.462000 2.464000 2.466000 2.468000 2.470000 2.472000 2.474000 2.476000 2.478000 2.480000 2.482000 2.484000 2.486000 2.488000 2.490000 2.492000 2.494000 2.496000 2.498000 2.500000'
  []
  nl_max_its = 50
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-3
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
[]

[Postprocessors]
  [rel_disp_crest_base_x]
    type = DifferencePostprocessor
    value1 = crest_disp_x
    value2 = base_disp_x
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_disp_x]
    type = NodalSum
    variable = disp_x
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_disp_x]
    type = NodalSum
    variable = disp_x
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [volume_total]
    type = VolumePostprocessor
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  [ip_probes]
    type = ElementValueSampler
    variable = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_yy strain_zz strain_xy damage_t damage_c kappa_t kappa_c'
    sort_by = id
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  file_base = ds04g_s4_eq1p5_out
  [history]
    type = CSV
    execute_vector_postprocessors_on = NONE
  []
  [fields_csv]
    type = CSV
    file_base = ds04g_s4_fields
    execute_postprocessors_on = NONE
    time_step_interval = 10
  []
  [fields_exodus]
    type = Exodus
    time_step_interval = 10
  []
  [cp]
    type = Checkpoint
    time_step_interval = 100
    num_files = 2
  []
[]

[Functions/base_acceleration]
  type = PiecewiseLinear
  data_file = elcentro_total_si.csv
  format = columns
[]

[Controls/phase]
  type = TimePeriod
  execute_on = 'INITIAL TIMESTEP_BEGIN'
  start_time = 0
  end_time = 1.000001
  enable_objects = '*/static_x */static_y */base_fix_x'
  disable_objects = '*/TM_all0 */TM_all1 */TM_all_inertia_x */TM_all_inertia_y */TM_all_accel_x */TM_all_accel_y */TM_all_vel_x */TM_all_vel_y */base_acceleration_x'
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
  [kappa_t]
    family = MONOMIAL
    order = CONSTANT
  []
  [kappa_c]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [damage_t_output]
    type = MaterialRealAux
    variable = damage_t
    property = DamageT
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damage_c_output]
    type = MaterialRealAux
    variable = damage_c
    property = DamageC
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_t_output]
    type = MaterialRealAux
    variable = kappa_t
    property = cdp_kappa_t
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [kappa_c_output]
    type = MaterialRealAux
    variable = kappa_c
    property = cdp_kappa_c
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [base_disp_y]
    type = NodalSum
    variable = disp_y
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_vel_x]
    type = NodalSum
    variable = vel_x
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_vel_y]
    type = NodalSum
    variable = vel_y
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_accel_x]
    type = NodalSum
    variable = accel_x
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [base_accel_y]
    type = NodalSum
    variable = accel_y
    boundary = 'BASE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_disp_y]
    type = NodalSum
    variable = disp_y
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_vel_x]
    type = NodalSum
    variable = vel_x
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_vel_y]
    type = NodalSum
    variable = vel_y
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_accel_x]
    type = NodalSum
    variable = accel_x
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [crest_accel_y]
    type = NodalSum
    variable = accel_y
    boundary = 'CREST__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_disp_x]
    type = NodalSum
    variable = disp_x
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_disp_y]
    type = NodalSum
    variable = disp_y
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_vel_x]
    type = NodalSum
    variable = vel_x
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_vel_y]
    type = NodalSum
    variable = vel_y
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_accel_x]
    type = NodalSum
    variable = accel_x
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [heel_accel_y]
    type = NodalSum
    variable = accel_y
    boundary = 'HEEL__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_disp_x]
    type = NodalSum
    variable = disp_x
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_disp_y]
    type = NodalSum
    variable = disp_y
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_vel_x]
    type = NodalSum
    variable = vel_x
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_vel_y]
    type = NodalSum
    variable = vel_y
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_accel_x]
    type = NodalSum
    variable = accel_x
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [toe_accel_y]
    type = NodalSum
    variable = accel_y
    boundary = 'TOE__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_disp_x]
    type = NodalSum
    variable = disp_x
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_disp_y]
    type = NodalSum
    variable = disp_y
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_vel_x]
    type = NodalSum
    variable = vel_x
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_vel_y]
    type = NodalSum
    variable = vel_y
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_accel_x]
    type = NodalSum
    variable = accel_x
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [corner_accel_y]
    type = NodalSum
    variable = accel_y
    boundary = 'CORNER__Dam_1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_damage_t]
    type = ElementExtremeValue
    variable = damage_t
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_damage_c]
    type = ElementExtremeValue
    variable = damage_c
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
