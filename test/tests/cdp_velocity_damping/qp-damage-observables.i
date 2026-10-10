# Explicit selected_qp data on a single element; NOT quadrature averages.
[AuxVariables/damage_t_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp0]
  type = MaterialRealAux
  variable = damage_t_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp0]
  type = ElementAverageValue
  variable = damage_t_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp1]
  type = MaterialRealAux
  variable = damage_t_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp1]
  type = ElementAverageValue
  variable = damage_t_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp2]
  type = MaterialRealAux
  variable = damage_t_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp2]
  type = ElementAverageValue
  variable = damage_t_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp3]
  type = MaterialRealAux
  variable = damage_t_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp3]
  type = ElementAverageValue
  variable = damage_t_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp4]
  type = MaterialRealAux
  variable = damage_t_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp4]
  type = ElementAverageValue
  variable = damage_t_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp5]
  type = MaterialRealAux
  variable = damage_t_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp5]
  type = ElementAverageValue
  variable = damage_t_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp6]
  type = MaterialRealAux
  variable = damage_t_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp6]
  type = ElementAverageValue
  variable = damage_t_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_t_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_t_qp7]
  type = MaterialRealAux
  variable = damage_t_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  property = DamageT
[]
[Postprocessors/damage_t_qp7]
  type = ElementAverageValue
  variable = damage_t_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp0]
  type = MaterialRealAux
  variable = damage_c_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp0]
  type = ElementAverageValue
  variable = damage_c_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp1]
  type = MaterialRealAux
  variable = damage_c_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp1]
  type = ElementAverageValue
  variable = damage_c_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp2]
  type = MaterialRealAux
  variable = damage_c_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp2]
  type = ElementAverageValue
  variable = damage_c_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp3]
  type = MaterialRealAux
  variable = damage_c_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp3]
  type = ElementAverageValue
  variable = damage_c_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp4]
  type = MaterialRealAux
  variable = damage_c_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp4]
  type = ElementAverageValue
  variable = damage_c_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp5]
  type = MaterialRealAux
  variable = damage_c_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp5]
  type = ElementAverageValue
  variable = damage_c_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp6]
  type = MaterialRealAux
  variable = damage_c_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp6]
  type = ElementAverageValue
  variable = damage_c_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/damage_c_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/damage_c_qp7]
  type = MaterialRealAux
  variable = damage_c_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  property = DamageC
[]
[Postprocessors/damage_c_qp7]
  type = ElementAverageValue
  variable = damage_c_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp0]
  type = MaterialRealAux
  variable = kappa_t_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp0]
  type = ElementAverageValue
  variable = kappa_t_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp1]
  type = MaterialRealAux
  variable = kappa_t_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp1]
  type = ElementAverageValue
  variable = kappa_t_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp2]
  type = MaterialRealAux
  variable = kappa_t_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp2]
  type = ElementAverageValue
  variable = kappa_t_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp3]
  type = MaterialRealAux
  variable = kappa_t_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp3]
  type = ElementAverageValue
  variable = kappa_t_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp4]
  type = MaterialRealAux
  variable = kappa_t_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp4]
  type = ElementAverageValue
  variable = kappa_t_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp5]
  type = MaterialRealAux
  variable = kappa_t_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp5]
  type = ElementAverageValue
  variable = kappa_t_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp6]
  type = MaterialRealAux
  variable = kappa_t_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp6]
  type = ElementAverageValue
  variable = kappa_t_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_t_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t_qp7]
  type = MaterialRealAux
  variable = kappa_t_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  property = cdp_kappa_t
[]
[Postprocessors/kappa_t_qp7]
  type = ElementAverageValue
  variable = kappa_t_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp0]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp0]
  type = MaterialRealAux
  variable = kappa_c_qp0
  selected_qp = 0
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp0]
  type = ElementAverageValue
  variable = kappa_c_qp0
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp1]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp1]
  type = MaterialRealAux
  variable = kappa_c_qp1
  selected_qp = 1
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp1]
  type = ElementAverageValue
  variable = kappa_c_qp1
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp2]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp2]
  type = MaterialRealAux
  variable = kappa_c_qp2
  selected_qp = 2
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp2]
  type = ElementAverageValue
  variable = kappa_c_qp2
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp3]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp3]
  type = MaterialRealAux
  variable = kappa_c_qp3
  selected_qp = 3
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp3]
  type = ElementAverageValue
  variable = kappa_c_qp3
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp4]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp4]
  type = MaterialRealAux
  variable = kappa_c_qp4
  selected_qp = 4
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp4]
  type = ElementAverageValue
  variable = kappa_c_qp4
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp5]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp5]
  type = MaterialRealAux
  variable = kappa_c_qp5
  selected_qp = 5
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp5]
  type = ElementAverageValue
  variable = kappa_c_qp5
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp6]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp6]
  type = MaterialRealAux
  variable = kappa_c_qp6
  selected_qp = 6
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp6]
  type = ElementAverageValue
  variable = kappa_c_qp6
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c_qp7]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c_qp7]
  type = MaterialRealAux
  variable = kappa_c_qp7
  selected_qp = 7
  execute_on = TIMESTEP_END
  property = cdp_kappa_c
[]
[Postprocessors/kappa_c_qp7]
  type = ElementAverageValue
  variable = kappa_c_qp7
  execute_on = 'INITIAL TIMESTEP_END'
[]
