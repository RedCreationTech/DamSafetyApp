# Interface development fixture, not the aqueduct model or an Abaqus reference.
# Pure end moment: M=1000 Nm, width=.2 m, length=.2 m, h=.1 m, nu=0.
# SPOS expected tip: ry=6e-5 rad, uz=-6e-6 m, ux=3e-6 m.
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = .2
  ymax = .2
  elem_type = QUAD4
[]
[GlobalParams]
  displacements = 'ux uy uz'
  rotations = 'rx ry rz'
[]
[Variables]
  [ux]
  []
  [uy]
  []
  [uz]
  []
  [rx]
  []
  [ry]
  []
  [rz]
  []
[]
[Kernels]
  [ux]
    type = AqueductFlatShell
    variable = ux
    component = 0
  []
  [uy]
    type = AqueductFlatShell
    variable = uy
    component = 1
  []
  [uz]
    type = AqueductFlatShell
    variable = uz
    component = 2
  []
  [rx]
    type = AqueductFlatShell
    variable = rx
    component = 3
  []
  [ry]
    type = AqueductFlatShell
    variable = ry
    component = 4
  []
  [rz]
    type = AqueductFlatShell
    variable = rz
    component = 5
  []
[]
[Materials]
  [shell]
    type = AqueductFlatShellStress
    youngs_modulus = 200e9
    poissons_ratio = 0
    yield_stress = 400e6
    ultimate_stress = 696e6
    plastic_strain_at_ultimate = .1
    thickness = .1
    offset = .5
    shear_correction = .8333333333333333
    drilling_factor = .001
  []
[]
[BCs]
  [ux]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  []
  [uy]
    type = DirichletBC
    variable = uy
    boundary = left
    value = 0
  []
  [uz]
    type = DirichletBC
    variable = uz
    boundary = left
    value = 0
  []
  [rx]
    type = DirichletBC
    variable = rx
    boundary = left
    value = 0
  []
  [ry]
    type = DirichletBC
    variable = ry
    boundary = left
    value = 0
  []
  [rz]
    type = DirichletBC
    variable = rz
    boundary = left
    value = 0
  []
[]
[NodalKernels]
  [moment]
    type = UserForcingFunctionNodalKernel
    variable = ry
    boundary = right
    function = 500
  []
[]
[Postprocessors]
  [tip_ux]
    type = PointValue
    variable = ux
    point = '.2 .1 0'
  []
  [tip_uz]
    type = PointValue
    variable = uz
    point = '.2 .1 0'
  []
  [tip_ry]
    type = PointValue
    variable = ry
    point = '.2 .1 0'
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[Outputs]
  exodus = true
  csv = true
[]
