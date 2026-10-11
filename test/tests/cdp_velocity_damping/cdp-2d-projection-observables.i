# Surrogate P: in-plane total rate and separate damping projection.
[AuxVariables/projection_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xx_qp0]
  type = RankTwoAux
  variable = projection_xx_qp0
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xx_qp0]
  type = ElementAverageValue
  variable = projection_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp0]
  type = RankTwoAux
  variable = damping_xx_qp0
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xx_qp0]
  type = ElementAverageValue
  variable = damping_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp0]
  type = RankTwoAux
  variable = rate_xx_qp0
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xx_qp0]
  type = ElementAverageValue
  variable = rate_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_yy_qp0]
  type = RankTwoAux
  variable = projection_yy_qp0
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_yy_qp0]
  type = ElementAverageValue
  variable = projection_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp0]
  type = RankTwoAux
  variable = damping_yy_qp0
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_yy_qp0]
  type = ElementAverageValue
  variable = damping_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp0]
  type = RankTwoAux
  variable = rate_yy_qp0
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_yy_qp0]
  type = ElementAverageValue
  variable = rate_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_zz_qp0]
  type = RankTwoAux
  variable = projection_zz_qp0
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 2
  index_j = 2
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_zz_qp0]
  type = ElementAverageValue
  variable = projection_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp0]
  type = RankTwoAux
  variable = damping_zz_qp0
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_zz_qp0]
  type = ElementAverageValue
  variable = damping_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xy_qp0]
  type = RankTwoAux
  variable = projection_xy_qp0
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xy_qp0]
  type = ElementAverageValue
  variable = projection_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp0]
  type = RankTwoAux
  variable = damping_xy_qp0
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xy_qp0]
  type = ElementAverageValue
  variable = damping_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp0]
  type = RankTwoAux
  variable = rate_xy_qp0
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xy_qp0]
  type = ElementAverageValue
  variable = rate_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/static_thickness_increment_over_dt_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/static_thickness_increment_over_dt_qp0]
  type = MaterialRealAux
  variable = static_thickness_increment_over_dt_qp0
  property = cdp_static_thickness_increment_over_dt
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/static_thickness_increment_over_dt_qp0]
  type = ElementAverageValue
  variable = static_thickness_increment_over_dt_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xx_qp1]
  type = RankTwoAux
  variable = projection_xx_qp1
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xx_qp1]
  type = ElementAverageValue
  variable = projection_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp1]
  type = RankTwoAux
  variable = damping_xx_qp1
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xx_qp1]
  type = ElementAverageValue
  variable = damping_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp1]
  type = RankTwoAux
  variable = rate_xx_qp1
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xx_qp1]
  type = ElementAverageValue
  variable = rate_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_yy_qp1]
  type = RankTwoAux
  variable = projection_yy_qp1
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_yy_qp1]
  type = ElementAverageValue
  variable = projection_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp1]
  type = RankTwoAux
  variable = damping_yy_qp1
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_yy_qp1]
  type = ElementAverageValue
  variable = damping_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp1]
  type = RankTwoAux
  variable = rate_yy_qp1
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_yy_qp1]
  type = ElementAverageValue
  variable = rate_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_zz_qp1]
  type = RankTwoAux
  variable = projection_zz_qp1
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 2
  index_j = 2
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_zz_qp1]
  type = ElementAverageValue
  variable = projection_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp1]
  type = RankTwoAux
  variable = damping_zz_qp1
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_zz_qp1]
  type = ElementAverageValue
  variable = damping_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xy_qp1]
  type = RankTwoAux
  variable = projection_xy_qp1
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xy_qp1]
  type = ElementAverageValue
  variable = projection_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp1]
  type = RankTwoAux
  variable = damping_xy_qp1
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xy_qp1]
  type = ElementAverageValue
  variable = damping_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp1]
  type = RankTwoAux
  variable = rate_xy_qp1
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xy_qp1]
  type = ElementAverageValue
  variable = rate_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/static_thickness_increment_over_dt_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/static_thickness_increment_over_dt_qp1]
  type = MaterialRealAux
  variable = static_thickness_increment_over_dt_qp1
  property = cdp_static_thickness_increment_over_dt
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/static_thickness_increment_over_dt_qp1]
  type = ElementAverageValue
  variable = static_thickness_increment_over_dt_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xx_qp2]
  type = RankTwoAux
  variable = projection_xx_qp2
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xx_qp2]
  type = ElementAverageValue
  variable = projection_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp2]
  type = RankTwoAux
  variable = damping_xx_qp2
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xx_qp2]
  type = ElementAverageValue
  variable = damping_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp2]
  type = RankTwoAux
  variable = rate_xx_qp2
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xx_qp2]
  type = ElementAverageValue
  variable = rate_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_yy_qp2]
  type = RankTwoAux
  variable = projection_yy_qp2
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_yy_qp2]
  type = ElementAverageValue
  variable = projection_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp2]
  type = RankTwoAux
  variable = damping_yy_qp2
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_yy_qp2]
  type = ElementAverageValue
  variable = damping_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp2]
  type = RankTwoAux
  variable = rate_yy_qp2
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_yy_qp2]
  type = ElementAverageValue
  variable = rate_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_zz_qp2]
  type = RankTwoAux
  variable = projection_zz_qp2
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 2
  index_j = 2
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_zz_qp2]
  type = ElementAverageValue
  variable = projection_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp2]
  type = RankTwoAux
  variable = damping_zz_qp2
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_zz_qp2]
  type = ElementAverageValue
  variable = damping_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xy_qp2]
  type = RankTwoAux
  variable = projection_xy_qp2
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xy_qp2]
  type = ElementAverageValue
  variable = projection_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp2]
  type = RankTwoAux
  variable = damping_xy_qp2
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xy_qp2]
  type = ElementAverageValue
  variable = damping_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp2]
  type = RankTwoAux
  variable = rate_xy_qp2
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xy_qp2]
  type = ElementAverageValue
  variable = rate_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/static_thickness_increment_over_dt_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/static_thickness_increment_over_dt_qp2]
  type = MaterialRealAux
  variable = static_thickness_increment_over_dt_qp2
  property = cdp_static_thickness_increment_over_dt
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/static_thickness_increment_over_dt_qp2]
  type = ElementAverageValue
  variable = static_thickness_increment_over_dt_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xx_qp3]
  type = RankTwoAux
  variable = projection_xx_qp3
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xx_qp3]
  type = ElementAverageValue
  variable = projection_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp3]
  type = RankTwoAux
  variable = damping_xx_qp3
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xx_qp3]
  type = ElementAverageValue
  variable = damping_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp3]
  type = RankTwoAux
  variable = rate_xx_qp3
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xx_qp3]
  type = ElementAverageValue
  variable = rate_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_yy_qp3]
  type = RankTwoAux
  variable = projection_yy_qp3
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_yy_qp3]
  type = ElementAverageValue
  variable = projection_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp3]
  type = RankTwoAux
  variable = damping_yy_qp3
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_yy_qp3]
  type = ElementAverageValue
  variable = damping_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp3]
  type = RankTwoAux
  variable = rate_yy_qp3
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_yy_qp3]
  type = ElementAverageValue
  variable = rate_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_zz_qp3]
  type = RankTwoAux
  variable = projection_zz_qp3
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 2
  index_j = 2
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_zz_qp3]
  type = ElementAverageValue
  variable = projection_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp3]
  type = RankTwoAux
  variable = damping_zz_qp3
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_zz_qp3]
  type = ElementAverageValue
  variable = damping_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/projection_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/projection_xy_qp3]
  type = RankTwoAux
  variable = projection_xy_qp3
  rank_two_tensor = cdp_damping_projected_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/projection_xy_qp3]
  type = ElementAverageValue
  variable = projection_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damping_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp3]
  type = RankTwoAux
  variable = damping_xy_qp3
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damping_xy_qp3]
  type = ElementAverageValue
  variable = damping_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/rate_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp3]
  type = RankTwoAux
  variable = rate_xy_qp3
  rank_two_tensor = cdp_inplane_total_velocity_rate
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/rate_xy_qp3]
  type = ElementAverageValue
  variable = rate_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/static_thickness_increment_over_dt_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/static_thickness_increment_over_dt_qp3]
  type = MaterialRealAux
  variable = static_thickness_increment_over_dt_qp3
  property = cdp_static_thickness_increment_over_dt
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/static_thickness_increment_over_dt_qp3]
  type = ElementAverageValue
  variable = static_thickness_increment_over_dt_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
