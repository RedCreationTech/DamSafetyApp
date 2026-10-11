# Actual selected Gauss points of a single full-integration 2D QUAD4. No Abaqus IP identity claim.
[AuxVariables/stress_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp0]
  type = RankTwoAux
  variable = stress_xx_qp0
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xx_qp0]
  type = ElementAverageValue
  variable = stress_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp0]
  type = RankTwoAux
  variable = strain_xx_qp0
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xx_qp0]
  type = ElementAverageValue
  variable = strain_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xx_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xx_qp0]
  type = RankTwoAux
  variable = pe_xx_qp0
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 0
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xx_qp0]
  type = ElementAverageValue
  variable = pe_xx_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp0]
  type = RankTwoAux
  variable = stress_yy_qp0
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_yy_qp0]
  type = ElementAverageValue
  variable = stress_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp0]
  type = RankTwoAux
  variable = strain_yy_qp0
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_yy_qp0]
  type = ElementAverageValue
  variable = strain_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_yy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_yy_qp0]
  type = RankTwoAux
  variable = pe_yy_qp0
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 1
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_yy_qp0]
  type = ElementAverageValue
  variable = pe_yy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp0]
  type = RankTwoAux
  variable = stress_zz_qp0
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_zz_qp0]
  type = ElementAverageValue
  variable = stress_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp0]
  type = RankTwoAux
  variable = strain_zz_qp0
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_zz_qp0]
  type = ElementAverageValue
  variable = strain_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_zz_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_zz_qp0]
  type = RankTwoAux
  variable = pe_zz_qp0
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 2
  index_j = 2
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_zz_qp0]
  type = ElementAverageValue
  variable = pe_zz_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp0]
  type = RankTwoAux
  variable = stress_xy_qp0
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xy_qp0]
  type = ElementAverageValue
  variable = stress_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp0]
  type = RankTwoAux
  variable = strain_xy_qp0
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xy_qp0]
  type = ElementAverageValue
  variable = strain_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xy_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xy_qp0]
  type = RankTwoAux
  variable = pe_xy_qp0
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 1
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xy_qp0]
  type = ElementAverageValue
  variable = pe_xy_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_t_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp0]
  type = MaterialRealAux
  variable = damage_t_qp0
  property = DamageT
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_t_qp0]
  type = ElementAverageValue
  variable = damage_t_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_c_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp0]
  type = MaterialRealAux
  variable = damage_c_qp0
  property = DamageC
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_c_qp0]
  type = ElementAverageValue
  variable = damage_c_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_t_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp0]
  type = MaterialRealAux
  variable = kappa_t_qp0
  property = cdp_kappa_t
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_t_qp0]
  type = ElementAverageValue
  variable = kappa_t_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_c_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp0]
  type = MaterialRealAux
  variable = kappa_c_qp0
  property = cdp_kappa_c
  selected_qp = 0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_c_qp0]
  type = ElementAverageValue
  variable = kappa_c_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp1]
  type = RankTwoAux
  variable = stress_xx_qp1
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xx_qp1]
  type = ElementAverageValue
  variable = stress_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp1]
  type = RankTwoAux
  variable = strain_xx_qp1
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xx_qp1]
  type = ElementAverageValue
  variable = strain_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xx_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xx_qp1]
  type = RankTwoAux
  variable = pe_xx_qp1
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 0
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xx_qp1]
  type = ElementAverageValue
  variable = pe_xx_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp1]
  type = RankTwoAux
  variable = stress_yy_qp1
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_yy_qp1]
  type = ElementAverageValue
  variable = stress_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp1]
  type = RankTwoAux
  variable = strain_yy_qp1
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_yy_qp1]
  type = ElementAverageValue
  variable = strain_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_yy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_yy_qp1]
  type = RankTwoAux
  variable = pe_yy_qp1
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 1
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_yy_qp1]
  type = ElementAverageValue
  variable = pe_yy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp1]
  type = RankTwoAux
  variable = stress_zz_qp1
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_zz_qp1]
  type = ElementAverageValue
  variable = stress_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp1]
  type = RankTwoAux
  variable = strain_zz_qp1
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_zz_qp1]
  type = ElementAverageValue
  variable = strain_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_zz_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_zz_qp1]
  type = RankTwoAux
  variable = pe_zz_qp1
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 2
  index_j = 2
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_zz_qp1]
  type = ElementAverageValue
  variable = pe_zz_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp1]
  type = RankTwoAux
  variable = stress_xy_qp1
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xy_qp1]
  type = ElementAverageValue
  variable = stress_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp1]
  type = RankTwoAux
  variable = strain_xy_qp1
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xy_qp1]
  type = ElementAverageValue
  variable = strain_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xy_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xy_qp1]
  type = RankTwoAux
  variable = pe_xy_qp1
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 1
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xy_qp1]
  type = ElementAverageValue
  variable = pe_xy_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_t_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp1]
  type = MaterialRealAux
  variable = damage_t_qp1
  property = DamageT
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_t_qp1]
  type = ElementAverageValue
  variable = damage_t_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_c_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp1]
  type = MaterialRealAux
  variable = damage_c_qp1
  property = DamageC
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_c_qp1]
  type = ElementAverageValue
  variable = damage_c_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_t_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp1]
  type = MaterialRealAux
  variable = kappa_t_qp1
  property = cdp_kappa_t
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_t_qp1]
  type = ElementAverageValue
  variable = kappa_t_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_c_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp1]
  type = MaterialRealAux
  variable = kappa_c_qp1
  property = cdp_kappa_c
  selected_qp = 1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_c_qp1]
  type = ElementAverageValue
  variable = kappa_c_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp2]
  type = RankTwoAux
  variable = stress_xx_qp2
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xx_qp2]
  type = ElementAverageValue
  variable = stress_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp2]
  type = RankTwoAux
  variable = strain_xx_qp2
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xx_qp2]
  type = ElementAverageValue
  variable = strain_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xx_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xx_qp2]
  type = RankTwoAux
  variable = pe_xx_qp2
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 0
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xx_qp2]
  type = ElementAverageValue
  variable = pe_xx_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp2]
  type = RankTwoAux
  variable = stress_yy_qp2
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_yy_qp2]
  type = ElementAverageValue
  variable = stress_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp2]
  type = RankTwoAux
  variable = strain_yy_qp2
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_yy_qp2]
  type = ElementAverageValue
  variable = strain_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_yy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_yy_qp2]
  type = RankTwoAux
  variable = pe_yy_qp2
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 1
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_yy_qp2]
  type = ElementAverageValue
  variable = pe_yy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp2]
  type = RankTwoAux
  variable = stress_zz_qp2
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_zz_qp2]
  type = ElementAverageValue
  variable = stress_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp2]
  type = RankTwoAux
  variable = strain_zz_qp2
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_zz_qp2]
  type = ElementAverageValue
  variable = strain_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_zz_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_zz_qp2]
  type = RankTwoAux
  variable = pe_zz_qp2
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 2
  index_j = 2
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_zz_qp2]
  type = ElementAverageValue
  variable = pe_zz_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp2]
  type = RankTwoAux
  variable = stress_xy_qp2
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xy_qp2]
  type = ElementAverageValue
  variable = stress_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp2]
  type = RankTwoAux
  variable = strain_xy_qp2
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xy_qp2]
  type = ElementAverageValue
  variable = strain_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xy_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xy_qp2]
  type = RankTwoAux
  variable = pe_xy_qp2
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 1
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xy_qp2]
  type = ElementAverageValue
  variable = pe_xy_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_t_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp2]
  type = MaterialRealAux
  variable = damage_t_qp2
  property = DamageT
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_t_qp2]
  type = ElementAverageValue
  variable = damage_t_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_c_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp2]
  type = MaterialRealAux
  variable = damage_c_qp2
  property = DamageC
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_c_qp2]
  type = ElementAverageValue
  variable = damage_c_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_t_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp2]
  type = MaterialRealAux
  variable = kappa_t_qp2
  property = cdp_kappa_t
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_t_qp2]
  type = ElementAverageValue
  variable = kappa_t_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_c_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp2]
  type = MaterialRealAux
  variable = kappa_c_qp2
  property = cdp_kappa_c
  selected_qp = 2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_c_qp2]
  type = ElementAverageValue
  variable = kappa_c_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xx_qp3]
  type = RankTwoAux
  variable = stress_xx_qp3
  rank_two_tensor = stress
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xx_qp3]
  type = ElementAverageValue
  variable = stress_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xx_qp3]
  type = RankTwoAux
  variable = strain_xx_qp3
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xx_qp3]
  type = ElementAverageValue
  variable = strain_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xx_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xx_qp3]
  type = RankTwoAux
  variable = pe_xx_qp3
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 0
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xx_qp3]
  type = ElementAverageValue
  variable = pe_xx_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_yy_qp3]
  type = RankTwoAux
  variable = stress_yy_qp3
  rank_two_tensor = stress
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_yy_qp3]
  type = ElementAverageValue
  variable = stress_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_yy_qp3]
  type = RankTwoAux
  variable = strain_yy_qp3
  rank_two_tensor = total_strain
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_yy_qp3]
  type = ElementAverageValue
  variable = strain_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_yy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_yy_qp3]
  type = RankTwoAux
  variable = pe_yy_qp3
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 1
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_yy_qp3]
  type = ElementAverageValue
  variable = pe_yy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_zz_qp3]
  type = RankTwoAux
  variable = stress_zz_qp3
  rank_two_tensor = stress
  index_i = 2
  index_j = 2
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_zz_qp3]
  type = ElementAverageValue
  variable = stress_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_zz_qp3]
  type = RankTwoAux
  variable = strain_zz_qp3
  rank_two_tensor = total_strain
  index_i = 2
  index_j = 2
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_zz_qp3]
  type = ElementAverageValue
  variable = strain_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_zz_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_zz_qp3]
  type = RankTwoAux
  variable = pe_zz_qp3
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 2
  index_j = 2
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_zz_qp3]
  type = ElementAverageValue
  variable = pe_zz_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/stress_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/stress_xy_qp3]
  type = RankTwoAux
  variable = stress_xy_qp3
  rank_two_tensor = stress
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/stress_xy_qp3]
  type = ElementAverageValue
  variable = stress_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/strain_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/strain_xy_qp3]
  type = RankTwoAux
  variable = strain_xy_qp3
  rank_two_tensor = total_strain
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/strain_xy_qp3]
  type = ElementAverageValue
  variable = strain_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/pe_xy_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/pe_xy_qp3]
  type = RankTwoAux
  variable = pe_xy_qp3
  rank_two_tensor = cdp_backbone_plastic_strain
  index_i = 0
  index_j = 1
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/pe_xy_qp3]
  type = ElementAverageValue
  variable = pe_xy_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_t_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp3]
  type = MaterialRealAux
  variable = damage_t_qp3
  property = DamageT
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_t_qp3]
  type = ElementAverageValue
  variable = damage_t_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/damage_c_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp3]
  type = MaterialRealAux
  variable = damage_c_qp3
  property = DamageC
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/damage_c_qp3]
  type = ElementAverageValue
  variable = damage_c_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_t_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp3]
  type = MaterialRealAux
  variable = kappa_t_qp3
  property = cdp_kappa_t
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_t_qp3]
  type = ElementAverageValue
  variable = kappa_t_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[AuxVariables/kappa_c_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp3]
  type = MaterialRealAux
  variable = kappa_c_qp3
  property = cdp_kappa_c
  selected_qp = 3
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_c_qp3]
  type = ElementAverageValue
  variable = kappa_c_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]
