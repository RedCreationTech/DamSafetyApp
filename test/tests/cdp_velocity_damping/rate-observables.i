[AuxVariables]
  [rate_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [rate_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [rate_zz]
    family = MONOMIAL
    order = CONSTANT
  []
  [rate_xy]
    family = MONOMIAL
    order = CONSTANT
  []
  [damping_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [damping_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [damping_zz]
    family = MONOMIAL
    order = CONSTANT
  []
  [damping_xy]
    family = MONOMIAL
    order = CONSTANT
  []
[]
[AuxKernels]
  [rate_xx]
    type = RankTwoAux
    variable = rate_xx
    rank_two_tensor = cdp_total_velocity_rate
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [rate_yy]
    type = RankTwoAux
    variable = rate_yy
    rank_two_tensor = cdp_total_velocity_rate
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [rate_zz]
    type = RankTwoAux
    variable = rate_zz
    rank_two_tensor = cdp_total_velocity_rate
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [rate_xy]
    type = RankTwoAux
    variable = rate_xy
    rank_two_tensor = cdp_total_velocity_rate
    index_i = 0
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [damping_xx]
    type = RankTwoAux
    variable = damping_xx
    rank_two_tensor = cdp_velocity_damping_stress
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  []
  [damping_yy]
    type = RankTwoAux
    variable = damping_yy
    rank_two_tensor = cdp_velocity_damping_stress
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  []
  [damping_zz]
    type = RankTwoAux
    variable = damping_zz
    rank_two_tensor = cdp_velocity_damping_stress
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  []
  [damping_xy]
    type = RankTwoAux
    variable = damping_xy
    rank_two_tensor = cdp_velocity_damping_stress
    index_i = 0
    index_j = 1
    execute_on = TIMESTEP_END
  []
[]
[Postprocessors]
  [rate_xx_mean]
    type = ElementAverageValue
    variable = rate_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rate_yy_mean]
    type = ElementAverageValue
    variable = rate_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rate_zz_mean]
    type = ElementAverageValue
    variable = rate_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [rate_xy_mean]
    type = ElementAverageValue
    variable = rate_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damping_xx_mean]
    type = ElementAverageValue
    variable = damping_xx
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damping_yy_mean]
    type = ElementAverageValue
    variable = damping_yy
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damping_zz_mean]
    type = ElementAverageValue
    variable = damping_zz
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damping_xy_mean]
    type = ElementAverageValue
    variable = damping_xy
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
