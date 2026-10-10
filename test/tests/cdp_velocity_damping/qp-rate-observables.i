# Explicit selected_qp data on a single element; NOT quadrature averages.
[AuxVariables/rate_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp0]
  type = RankTwoAux
  variable = rate_xx_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp0]
  type = ElementAverageValue
  variable = rate_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp1]
  type = RankTwoAux
  variable = rate_xx_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp1]
  type = ElementAverageValue
  variable = rate_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp2]
  type = RankTwoAux
  variable = rate_xx_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp2]
  type = ElementAverageValue
  variable = rate_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp3]
  type = RankTwoAux
  variable = rate_xx_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp3]
  type = ElementAverageValue
  variable = rate_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp4]
  type = RankTwoAux
  variable = rate_xx_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp4]
  type = ElementAverageValue
  variable = rate_xx_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp5]
  type = RankTwoAux
  variable = rate_xx_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp5]
  type = ElementAverageValue
  variable = rate_xx_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp6]
  type = RankTwoAux
  variable = rate_xx_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp6]
  type = ElementAverageValue
  variable = rate_xx_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xx_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xx_qp7]
  type = RankTwoAux
  variable = rate_xx_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 0
[]
[Postprocessors/rate_xx_qp7]
  type = ElementAverageValue
  variable = rate_xx_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp0]
  type = RankTwoAux
  variable = rate_yy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp0]
  type = ElementAverageValue
  variable = rate_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp1]
  type = RankTwoAux
  variable = rate_yy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp1]
  type = ElementAverageValue
  variable = rate_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp2]
  type = RankTwoAux
  variable = rate_yy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp2]
  type = ElementAverageValue
  variable = rate_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp3]
  type = RankTwoAux
  variable = rate_yy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp3]
  type = ElementAverageValue
  variable = rate_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp4]
  type = RankTwoAux
  variable = rate_yy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp4]
  type = ElementAverageValue
  variable = rate_yy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp5]
  type = RankTwoAux
  variable = rate_yy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp5]
  type = ElementAverageValue
  variable = rate_yy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp6]
  type = RankTwoAux
  variable = rate_yy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp6]
  type = ElementAverageValue
  variable = rate_yy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yy_qp7]
  type = RankTwoAux
  variable = rate_yy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 1
[]
[Postprocessors/rate_yy_qp7]
  type = ElementAverageValue
  variable = rate_yy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp0]
  type = RankTwoAux
  variable = rate_zz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp0]
  type = ElementAverageValue
  variable = rate_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp1]
  type = RankTwoAux
  variable = rate_zz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp1]
  type = ElementAverageValue
  variable = rate_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp2]
  type = RankTwoAux
  variable = rate_zz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp2]
  type = ElementAverageValue
  variable = rate_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp3]
  type = RankTwoAux
  variable = rate_zz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp3]
  type = ElementAverageValue
  variable = rate_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp4]
  type = RankTwoAux
  variable = rate_zz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp4]
  type = ElementAverageValue
  variable = rate_zz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp5]
  type = RankTwoAux
  variable = rate_zz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp5]
  type = ElementAverageValue
  variable = rate_zz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp6]
  type = RankTwoAux
  variable = rate_zz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp6]
  type = ElementAverageValue
  variable = rate_zz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_zz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_zz_qp7]
  type = RankTwoAux
  variable = rate_zz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 2
  index_j = 2
[]
[Postprocessors/rate_zz_qp7]
  type = ElementAverageValue
  variable = rate_zz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp0]
  type = RankTwoAux
  variable = rate_xy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp0]
  type = ElementAverageValue
  variable = rate_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp1]
  type = RankTwoAux
  variable = rate_xy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp1]
  type = ElementAverageValue
  variable = rate_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp2]
  type = RankTwoAux
  variable = rate_xy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp2]
  type = ElementAverageValue
  variable = rate_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp3]
  type = RankTwoAux
  variable = rate_xy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp3]
  type = ElementAverageValue
  variable = rate_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp4]
  type = RankTwoAux
  variable = rate_xy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp4]
  type = ElementAverageValue
  variable = rate_xy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp5]
  type = RankTwoAux
  variable = rate_xy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp5]
  type = ElementAverageValue
  variable = rate_xy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp6]
  type = RankTwoAux
  variable = rate_xy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp6]
  type = ElementAverageValue
  variable = rate_xy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xy_qp7]
  type = RankTwoAux
  variable = rate_xy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 1
[]
[Postprocessors/rate_xy_qp7]
  type = ElementAverageValue
  variable = rate_xy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp0]
  type = RankTwoAux
  variable = rate_xz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp0]
  type = ElementAverageValue
  variable = rate_xz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp1]
  type = RankTwoAux
  variable = rate_xz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp1]
  type = ElementAverageValue
  variable = rate_xz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp2]
  type = RankTwoAux
  variable = rate_xz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp2]
  type = ElementAverageValue
  variable = rate_xz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp3]
  type = RankTwoAux
  variable = rate_xz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp3]
  type = ElementAverageValue
  variable = rate_xz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp4]
  type = RankTwoAux
  variable = rate_xz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp4]
  type = ElementAverageValue
  variable = rate_xz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp5]
  type = RankTwoAux
  variable = rate_xz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp5]
  type = ElementAverageValue
  variable = rate_xz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp6]
  type = RankTwoAux
  variable = rate_xz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp6]
  type = ElementAverageValue
  variable = rate_xz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_xz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_xz_qp7]
  type = RankTwoAux
  variable = rate_xz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 0
  index_j = 2
[]
[Postprocessors/rate_xz_qp7]
  type = ElementAverageValue
  variable = rate_xz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp0]
  type = RankTwoAux
  variable = rate_yz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp0]
  type = ElementAverageValue
  variable = rate_yz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp1]
  type = RankTwoAux
  variable = rate_yz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp1]
  type = ElementAverageValue
  variable = rate_yz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp2]
  type = RankTwoAux
  variable = rate_yz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp2]
  type = ElementAverageValue
  variable = rate_yz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp3]
  type = RankTwoAux
  variable = rate_yz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp3]
  type = ElementAverageValue
  variable = rate_yz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp4]
  type = RankTwoAux
  variable = rate_yz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp4]
  type = ElementAverageValue
  variable = rate_yz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp5]
  type = RankTwoAux
  variable = rate_yz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp5]
  type = ElementAverageValue
  variable = rate_yz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp6]
  type = RankTwoAux
  variable = rate_yz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp6]
  type = ElementAverageValue
  variable = rate_yz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/rate_yz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/rate_yz_qp7]
  type = RankTwoAux
  variable = rate_yz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_total_velocity_rate
  index_i = 1
  index_j = 2
[]
[Postprocessors/rate_yz_qp7]
  type = ElementAverageValue
  variable = rate_yz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp0]
  type = RankTwoAux
  variable = damping_xx_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp0]
  type = ElementAverageValue
  variable = damping_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp1]
  type = RankTwoAux
  variable = damping_xx_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp1]
  type = ElementAverageValue
  variable = damping_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp2]
  type = RankTwoAux
  variable = damping_xx_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp2]
  type = ElementAverageValue
  variable = damping_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp3]
  type = RankTwoAux
  variable = damping_xx_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp3]
  type = ElementAverageValue
  variable = damping_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp4]
  type = RankTwoAux
  variable = damping_xx_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp4]
  type = ElementAverageValue
  variable = damping_xx_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp5]
  type = RankTwoAux
  variable = damping_xx_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp5]
  type = ElementAverageValue
  variable = damping_xx_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp6]
  type = RankTwoAux
  variable = damping_xx_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp6]
  type = ElementAverageValue
  variable = damping_xx_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xx_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xx_qp7]
  type = RankTwoAux
  variable = damping_xx_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/damping_xx_qp7]
  type = ElementAverageValue
  variable = damping_xx_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp0]
  type = RankTwoAux
  variable = damping_yy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp0]
  type = ElementAverageValue
  variable = damping_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp1]
  type = RankTwoAux
  variable = damping_yy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp1]
  type = ElementAverageValue
  variable = damping_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp2]
  type = RankTwoAux
  variable = damping_yy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp2]
  type = ElementAverageValue
  variable = damping_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp3]
  type = RankTwoAux
  variable = damping_yy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp3]
  type = ElementAverageValue
  variable = damping_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp4]
  type = RankTwoAux
  variable = damping_yy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp4]
  type = ElementAverageValue
  variable = damping_yy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp5]
  type = RankTwoAux
  variable = damping_yy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp5]
  type = ElementAverageValue
  variable = damping_yy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp6]
  type = RankTwoAux
  variable = damping_yy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp6]
  type = ElementAverageValue
  variable = damping_yy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yy_qp7]
  type = RankTwoAux
  variable = damping_yy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/damping_yy_qp7]
  type = ElementAverageValue
  variable = damping_yy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp0]
  type = RankTwoAux
  variable = damping_zz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp0]
  type = ElementAverageValue
  variable = damping_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp1]
  type = RankTwoAux
  variable = damping_zz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp1]
  type = ElementAverageValue
  variable = damping_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp2]
  type = RankTwoAux
  variable = damping_zz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp2]
  type = ElementAverageValue
  variable = damping_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp3]
  type = RankTwoAux
  variable = damping_zz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp3]
  type = ElementAverageValue
  variable = damping_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp4]
  type = RankTwoAux
  variable = damping_zz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp4]
  type = ElementAverageValue
  variable = damping_zz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp5]
  type = RankTwoAux
  variable = damping_zz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp5]
  type = ElementAverageValue
  variable = damping_zz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp6]
  type = RankTwoAux
  variable = damping_zz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp6]
  type = ElementAverageValue
  variable = damping_zz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_zz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_zz_qp7]
  type = RankTwoAux
  variable = damping_zz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/damping_zz_qp7]
  type = ElementAverageValue
  variable = damping_zz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp0]
  type = RankTwoAux
  variable = damping_xy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp0]
  type = ElementAverageValue
  variable = damping_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp1]
  type = RankTwoAux
  variable = damping_xy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp1]
  type = ElementAverageValue
  variable = damping_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp2]
  type = RankTwoAux
  variable = damping_xy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp2]
  type = ElementAverageValue
  variable = damping_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp3]
  type = RankTwoAux
  variable = damping_xy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp3]
  type = ElementAverageValue
  variable = damping_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp4]
  type = RankTwoAux
  variable = damping_xy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp4]
  type = ElementAverageValue
  variable = damping_xy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp5]
  type = RankTwoAux
  variable = damping_xy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp5]
  type = ElementAverageValue
  variable = damping_xy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp6]
  type = RankTwoAux
  variable = damping_xy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp6]
  type = ElementAverageValue
  variable = damping_xy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xy_qp7]
  type = RankTwoAux
  variable = damping_xy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/damping_xy_qp7]
  type = ElementAverageValue
  variable = damping_xy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp0]
  type = RankTwoAux
  variable = damping_xz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp0]
  type = ElementAverageValue
  variable = damping_xz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp1]
  type = RankTwoAux
  variable = damping_xz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp1]
  type = ElementAverageValue
  variable = damping_xz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp2]
  type = RankTwoAux
  variable = damping_xz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp2]
  type = ElementAverageValue
  variable = damping_xz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp3]
  type = RankTwoAux
  variable = damping_xz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp3]
  type = ElementAverageValue
  variable = damping_xz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp4]
  type = RankTwoAux
  variable = damping_xz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp4]
  type = ElementAverageValue
  variable = damping_xz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp5]
  type = RankTwoAux
  variable = damping_xz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp5]
  type = ElementAverageValue
  variable = damping_xz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp6]
  type = RankTwoAux
  variable = damping_xz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp6]
  type = ElementAverageValue
  variable = damping_xz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_xz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_xz_qp7]
  type = RankTwoAux
  variable = damping_xz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/damping_xz_qp7]
  type = ElementAverageValue
  variable = damping_xz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp0]
  type = RankTwoAux
  variable = damping_yz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp0]
  type = ElementAverageValue
  variable = damping_yz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp1]
  type = RankTwoAux
  variable = damping_yz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp1]
  type = ElementAverageValue
  variable = damping_yz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp2]
  type = RankTwoAux
  variable = damping_yz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp2]
  type = ElementAverageValue
  variable = damping_yz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp3]
  type = RankTwoAux
  variable = damping_yz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp3]
  type = ElementAverageValue
  variable = damping_yz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp4]
  type = RankTwoAux
  variable = damping_yz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp4]
  type = ElementAverageValue
  variable = damping_yz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp5]
  type = RankTwoAux
  variable = damping_yz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp5]
  type = ElementAverageValue
  variable = damping_yz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp6]
  type = RankTwoAux
  variable = damping_yz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp6]
  type = ElementAverageValue
  variable = damping_yz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damping_yz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damping_yz_qp7]
  type = RankTwoAux
  variable = damping_yz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = cdp_velocity_damping_stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/damping_yz_qp7]
  type = ElementAverageValue
  variable = damping_yz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]
