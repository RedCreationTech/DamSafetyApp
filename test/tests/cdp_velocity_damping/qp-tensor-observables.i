# Explicit selected_qp data on a single element; NOT quadrature averages.
[AuxVariables/stress_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp0]
  type = RankTwoAux
  variable = stress_xx_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp0]
  type = ElementAverageValue
  variable = stress_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp1]
  type = RankTwoAux
  variable = stress_xx_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp1]
  type = ElementAverageValue
  variable = stress_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp2]
  type = RankTwoAux
  variable = stress_xx_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp2]
  type = ElementAverageValue
  variable = stress_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp3]
  type = RankTwoAux
  variable = stress_xx_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp3]
  type = ElementAverageValue
  variable = stress_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp4]
  type = RankTwoAux
  variable = stress_xx_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp4]
  type = ElementAverageValue
  variable = stress_xx_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp5]
  type = RankTwoAux
  variable = stress_xx_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp5]
  type = ElementAverageValue
  variable = stress_xx_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp6]
  type = RankTwoAux
  variable = stress_xx_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp6]
  type = ElementAverageValue
  variable = stress_xx_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xx_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp7]
  type = RankTwoAux
  variable = stress_xx_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
[]
[Postprocessors/stress_xx_qp7]
  type = ElementAverageValue
  variable = stress_xx_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp0]
  type = RankTwoAux
  variable = stress_yy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp0]
  type = ElementAverageValue
  variable = stress_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp1]
  type = RankTwoAux
  variable = stress_yy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp1]
  type = ElementAverageValue
  variable = stress_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp2]
  type = RankTwoAux
  variable = stress_yy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp2]
  type = ElementAverageValue
  variable = stress_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp3]
  type = RankTwoAux
  variable = stress_yy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp3]
  type = ElementAverageValue
  variable = stress_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp4]
  type = RankTwoAux
  variable = stress_yy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp4]
  type = ElementAverageValue
  variable = stress_yy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp5]
  type = RankTwoAux
  variable = stress_yy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp5]
  type = ElementAverageValue
  variable = stress_yy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp6]
  type = RankTwoAux
  variable = stress_yy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp6]
  type = ElementAverageValue
  variable = stress_yy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp7]
  type = RankTwoAux
  variable = stress_yy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
[]
[Postprocessors/stress_yy_qp7]
  type = ElementAverageValue
  variable = stress_yy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp0]
  type = RankTwoAux
  variable = stress_zz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp0]
  type = ElementAverageValue
  variable = stress_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp1]
  type = RankTwoAux
  variable = stress_zz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp1]
  type = ElementAverageValue
  variable = stress_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp2]
  type = RankTwoAux
  variable = stress_zz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp2]
  type = ElementAverageValue
  variable = stress_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp3]
  type = RankTwoAux
  variable = stress_zz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp3]
  type = ElementAverageValue
  variable = stress_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp4]
  type = RankTwoAux
  variable = stress_zz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp4]
  type = ElementAverageValue
  variable = stress_zz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp5]
  type = RankTwoAux
  variable = stress_zz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp5]
  type = ElementAverageValue
  variable = stress_zz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp6]
  type = RankTwoAux
  variable = stress_zz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp6]
  type = ElementAverageValue
  variable = stress_zz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_zz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp7]
  type = RankTwoAux
  variable = stress_zz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
[]
[Postprocessors/stress_zz_qp7]
  type = ElementAverageValue
  variable = stress_zz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp0]
  type = RankTwoAux
  variable = stress_xy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp0]
  type = ElementAverageValue
  variable = stress_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp1]
  type = RankTwoAux
  variable = stress_xy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp1]
  type = ElementAverageValue
  variable = stress_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp2]
  type = RankTwoAux
  variable = stress_xy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp2]
  type = ElementAverageValue
  variable = stress_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp3]
  type = RankTwoAux
  variable = stress_xy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp3]
  type = ElementAverageValue
  variable = stress_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp4]
  type = RankTwoAux
  variable = stress_xy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp4]
  type = ElementAverageValue
  variable = stress_xy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp5]
  type = RankTwoAux
  variable = stress_xy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp5]
  type = ElementAverageValue
  variable = stress_xy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp6]
  type = RankTwoAux
  variable = stress_xy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp6]
  type = ElementAverageValue
  variable = stress_xy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp7]
  type = RankTwoAux
  variable = stress_xy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
[]
[Postprocessors/stress_xy_qp7]
  type = ElementAverageValue
  variable = stress_xy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp0]
  type = RankTwoAux
  variable = stress_xz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp0]
  type = ElementAverageValue
  variable = stress_xz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp1]
  type = RankTwoAux
  variable = stress_xz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp1]
  type = ElementAverageValue
  variable = stress_xz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp2]
  type = RankTwoAux
  variable = stress_xz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp2]
  type = ElementAverageValue
  variable = stress_xz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp3]
  type = RankTwoAux
  variable = stress_xz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp3]
  type = ElementAverageValue
  variable = stress_xz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp4]
  type = RankTwoAux
  variable = stress_xz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp4]
  type = ElementAverageValue
  variable = stress_xz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp5]
  type = RankTwoAux
  variable = stress_xz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp5]
  type = ElementAverageValue
  variable = stress_xz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp6]
  type = RankTwoAux
  variable = stress_xz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp6]
  type = ElementAverageValue
  variable = stress_xz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_xz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xz_qp7]
  type = RankTwoAux
  variable = stress_xz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 0
  index_j = 2
[]
[Postprocessors/stress_xz_qp7]
  type = ElementAverageValue
  variable = stress_xz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp0]
  type = RankTwoAux
  variable = stress_yz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp0]
  type = ElementAverageValue
  variable = stress_yz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp1]
  type = RankTwoAux
  variable = stress_yz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp1]
  type = ElementAverageValue
  variable = stress_yz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp2]
  type = RankTwoAux
  variable = stress_yz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp2]
  type = ElementAverageValue
  variable = stress_yz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp3]
  type = RankTwoAux
  variable = stress_yz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp3]
  type = ElementAverageValue
  variable = stress_yz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp4]
  type = RankTwoAux
  variable = stress_yz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp4]
  type = ElementAverageValue
  variable = stress_yz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp5]
  type = RankTwoAux
  variable = stress_yz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp5]
  type = ElementAverageValue
  variable = stress_yz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp6]
  type = RankTwoAux
  variable = stress_yz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp6]
  type = ElementAverageValue
  variable = stress_yz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/stress_yz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yz_qp7]
  type = RankTwoAux
  variable = stress_yz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = stress
  index_i = 1
  index_j = 2
[]
[Postprocessors/stress_yz_qp7]
  type = ElementAverageValue
  variable = stress_yz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp0]
  type = RankTwoAux
  variable = strain_xx_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp0]
  type = ElementAverageValue
  variable = strain_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp1]
  type = RankTwoAux
  variable = strain_xx_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp1]
  type = ElementAverageValue
  variable = strain_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp2]
  type = RankTwoAux
  variable = strain_xx_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp2]
  type = ElementAverageValue
  variable = strain_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp3]
  type = RankTwoAux
  variable = strain_xx_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp3]
  type = ElementAverageValue
  variable = strain_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp4]
  type = RankTwoAux
  variable = strain_xx_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp4]
  type = ElementAverageValue
  variable = strain_xx_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp5]
  type = RankTwoAux
  variable = strain_xx_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp5]
  type = ElementAverageValue
  variable = strain_xx_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp6]
  type = RankTwoAux
  variable = strain_xx_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp6]
  type = ElementAverageValue
  variable = strain_xx_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xx_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp7]
  type = RankTwoAux
  variable = strain_xx_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
[]
[Postprocessors/strain_xx_qp7]
  type = ElementAverageValue
  variable = strain_xx_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp0]
  type = RankTwoAux
  variable = strain_yy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp0]
  type = ElementAverageValue
  variable = strain_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp1]
  type = RankTwoAux
  variable = strain_yy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp1]
  type = ElementAverageValue
  variable = strain_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp2]
  type = RankTwoAux
  variable = strain_yy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp2]
  type = ElementAverageValue
  variable = strain_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp3]
  type = RankTwoAux
  variable = strain_yy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp3]
  type = ElementAverageValue
  variable = strain_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp4]
  type = RankTwoAux
  variable = strain_yy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp4]
  type = ElementAverageValue
  variable = strain_yy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp5]
  type = RankTwoAux
  variable = strain_yy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp5]
  type = ElementAverageValue
  variable = strain_yy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp6]
  type = RankTwoAux
  variable = strain_yy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp6]
  type = ElementAverageValue
  variable = strain_yy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp7]
  type = RankTwoAux
  variable = strain_yy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
[]
[Postprocessors/strain_yy_qp7]
  type = ElementAverageValue
  variable = strain_yy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp0]
  type = RankTwoAux
  variable = strain_zz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp0]
  type = ElementAverageValue
  variable = strain_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp1]
  type = RankTwoAux
  variable = strain_zz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp1]
  type = ElementAverageValue
  variable = strain_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp2]
  type = RankTwoAux
  variable = strain_zz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp2]
  type = ElementAverageValue
  variable = strain_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp3]
  type = RankTwoAux
  variable = strain_zz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp3]
  type = ElementAverageValue
  variable = strain_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp4]
  type = RankTwoAux
  variable = strain_zz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp4]
  type = ElementAverageValue
  variable = strain_zz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp5]
  type = RankTwoAux
  variable = strain_zz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp5]
  type = ElementAverageValue
  variable = strain_zz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp6]
  type = RankTwoAux
  variable = strain_zz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp6]
  type = ElementAverageValue
  variable = strain_zz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_zz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp7]
  type = RankTwoAux
  variable = strain_zz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
[]
[Postprocessors/strain_zz_qp7]
  type = ElementAverageValue
  variable = strain_zz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp0]
  type = RankTwoAux
  variable = strain_xy_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp0]
  type = ElementAverageValue
  variable = strain_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp1]
  type = RankTwoAux
  variable = strain_xy_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp1]
  type = ElementAverageValue
  variable = strain_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp2]
  type = RankTwoAux
  variable = strain_xy_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp2]
  type = ElementAverageValue
  variable = strain_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp3]
  type = RankTwoAux
  variable = strain_xy_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp3]
  type = ElementAverageValue
  variable = strain_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp4]
  type = RankTwoAux
  variable = strain_xy_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp4]
  type = ElementAverageValue
  variable = strain_xy_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp5]
  type = RankTwoAux
  variable = strain_xy_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp5]
  type = ElementAverageValue
  variable = strain_xy_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp6]
  type = RankTwoAux
  variable = strain_xy_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp6]
  type = ElementAverageValue
  variable = strain_xy_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xy_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp7]
  type = RankTwoAux
  variable = strain_xy_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
[]
[Postprocessors/strain_xy_qp7]
  type = ElementAverageValue
  variable = strain_xy_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp0]
  type = RankTwoAux
  variable = strain_xz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp0]
  type = ElementAverageValue
  variable = strain_xz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp1]
  type = RankTwoAux
  variable = strain_xz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp1]
  type = ElementAverageValue
  variable = strain_xz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp2]
  type = RankTwoAux
  variable = strain_xz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp2]
  type = ElementAverageValue
  variable = strain_xz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp3]
  type = RankTwoAux
  variable = strain_xz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp3]
  type = ElementAverageValue
  variable = strain_xz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp4]
  type = RankTwoAux
  variable = strain_xz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp4]
  type = ElementAverageValue
  variable = strain_xz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp5]
  type = RankTwoAux
  variable = strain_xz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp5]
  type = ElementAverageValue
  variable = strain_xz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp6]
  type = RankTwoAux
  variable = strain_xz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp6]
  type = ElementAverageValue
  variable = strain_xz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_xz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xz_qp7]
  type = RankTwoAux
  variable = strain_xz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 2
[]
[Postprocessors/strain_xz_qp7]
  type = ElementAverageValue
  variable = strain_xz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp0]
  type = RankTwoAux
  variable = strain_yz_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp0]
  type = ElementAverageValue
  variable = strain_yz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp1]
  type = RankTwoAux
  variable = strain_yz_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp1]
  type = ElementAverageValue
  variable = strain_yz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp2]
  type = RankTwoAux
  variable = strain_yz_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp2]
  type = ElementAverageValue
  variable = strain_yz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp3]
  type = RankTwoAux
  variable = strain_yz_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp3]
  type = ElementAverageValue
  variable = strain_yz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp4]
  type = RankTwoAux
  variable = strain_yz_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp4]
  type = ElementAverageValue
  variable = strain_yz_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp5]
  type = RankTwoAux
  variable = strain_yz_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp5]
  type = ElementAverageValue
  variable = strain_yz_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp6]
  type = RankTwoAux
  variable = strain_yz_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp6]
  type = ElementAverageValue
  variable = strain_yz_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/strain_yz_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yz_qp7]
  type = RankTwoAux
  variable = strain_yz_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 2
[]
[Postprocessors/strain_yz_qp7]
  type = ElementAverageValue
  variable = strain_yz_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]
