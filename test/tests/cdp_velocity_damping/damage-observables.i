[AuxVariables]
  [damage_t]
    family = MONOMIAL
    order = CONSTANT
  []
  [damage_c]
    family = MONOMIAL
    order = CONSTANT
  []
  [stiffness_factor]
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
  [stiffness_factor]
    type = MaterialRealAux
    variable = stiffness_factor
    property = cdp_stiffness_factor
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Postprocessors]
  [damage_t_mean]
    type = ElementAverageValue
    variable = damage_t
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [damage_c_mean]
    type = ElementAverageValue
    variable = damage_c
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [stiffness_factor_mean]
    type = ElementAverageValue
    variable = stiffness_factor
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[AuxVariables/kappa_t]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_t]
  type = MaterialRealAux
  variable = kappa_t
  property = cdp_kappa_t
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_t_mean]
  type = ElementAverageValue
  variable = kappa_t
  execute_on = 'INITIAL TIMESTEP_END'
[]

[AuxVariables/kappa_c]
  family = MONOMIAL
  order = CONSTANT
[]
[AuxKernels/kappa_c]
  type = MaterialRealAux
  variable = kappa_c
  property = cdp_kappa_c
  execute_on = 'INITIAL TIMESTEP_END'
[]
[Postprocessors/kappa_c_mean]
  type = ElementAverageValue
  variable = kappa_c
  execute_on = 'INITIAL TIMESTEP_END'
[]
