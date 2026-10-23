#!/usr/bin/env python3
"""Independent cyclic steel and offset shell-mass fixtures, never the full model."""
import argparse
import json
from pathlib import Path

from build_aqueduct_interface_cases import block


FIELDS = ('ux', 'uy', 'uz', 'rx', 'ry', 'rz')


def execution(end, dt, integrator=''):
    return block('Preconditioning', block('smp', 'type = SMP\nfull = true')) + block(
        'Executioner', 'type = Transient\nsolve_type = NEWTON\nend_time = '+str(end)+
        '\ndt = '+str(dt)+ '\nnl_abs_tol = 1e-7\nnl_rel_tol = 1e-11\nline_search = basic\n'
        "petsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'\n"+integrator) + block(
        'Outputs', 'csv = true\nexodus = true')


def truss_cycle():
    text = block('Mesh', 'type = GeneratedMesh\ndim = 1\nnx = 1\nxmax = 1\nallow_renumbering = false')
    text += block('Variables', ''.join(block(v, '') for v in FIELDS[:3]))
    props = ('axial_stress', 'signed_plastic_stretch', 'equivalent_plastic_stretch')
    text += block('AuxVariables', block('area', 'family = MONOMIAL\norder = CONSTANT\ninitial_condition = .001')+
                  ''.join(block(v, 'family = MONOMIAL\norder = CONSTANT') for v in props))
    text += block('AuxKernels', ''.join(block(v, f'type = MaterialRealAux\nvariable = {v}\nproperty = {v}\nexecute_on = TIMESTEP_END') for v in props))
    text += block('Materials', block('steel', "type = AqueductSteelTruss\ndisplacements = 'ux uy uz'\nyoungs_modulus = 200e9\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1"))
    text += block('Kernels', ''.join(block(v, f"type = StressDivergenceTensorsTruss\nvariable = {v}\ncomponent = {j}\ndisplacements = 'ux uy uz'\narea = area\nuse_displaced_mesh = false") for j, v in enumerate(FIELDS[:3])))
    text += block('BCs', ''.join(block(v, f"type = DirichletBC\nvariable = {v}\nboundary = '{'left' if v == 'ux' else 'left right'}'\nvalue = 0") for v in FIELDS[:3]))
    text += block('Functions', block('force', "type = PiecewiseLinear\nx = '0 .5 1 1.5 2'\ny = '0 450000 0 -500000 0'"))
    text += block('NodalKernels', block('force', 'type = UserForcingFunctorNodalKernel\nvariable = ux\nboundary = right\nfunctor = force'))
    text += block('Postprocessors', block('right_u', 'type = NodalVariableValue\nvariable = ux\nnodeid = 1')+
                  ''.join(block(v, f'type = ElementAverageValue\nvariable = {v}') for v in props))
    return text+execution(2, .05)


def shell_mass(component):
    # The first nonzero Newmark displacement is a*dt^2/4. Subsequent samples
    # maintain the same discrete acceleration from zero initial velocity/accel.
    acceleration = .01 if component == 0 else .0001
    text = block('Mesh', 'type = GeneratedMesh\ndim = 2\nnx = 1\nny = 1\nxmax = .2\nymax = .2\nelem_type = QUAD4')
    text += block('GlobalParams', "displacements = 'ux uy uz'\nrotations = 'rx ry rz'")
    text += block('Variables', ''.join(block(v, '') for v in FIELDS))
    text += block('AuxVariables', ''.join(block('mass_'+v, '') for v in FIELDS))
    text += block('Kernels', ''.join(block(v, f'type = AqueductFlatShellInertia\nvariable = {v}\ncomponent = {j}\ndensity = 7850\nthickness = .1\noffset = .5\nsave_in = mass_{v}') for j, v in enumerate(FIELDS)))
    text += block('Functions', block('motion', f"type = ParsedFunction\nexpression = 'if(t<.05,0,{acceleration}*(.5*t*t-.05*t+.0025))'"))
    text += block('BCs', ''.join(block(v, f"type = FunctionDirichletBC\nvariable = {v}\nboundary = 'left right'\nfunction = {'motion' if j == component else '0'}") for j, v in enumerate(FIELDS)))
    text += block('Postprocessors', ''.join(block('mass_'+v, f"type = NodalSum\nvariable = mass_{v}\nboundary = 'left right'") for v in FIELDS))
    return text+execution(.3, .1, block('TimeIntegrator', 'type = NewmarkBeta\nbeta = .25\ngamma = .5'))


def build(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    for name, text in {'steel_truss_cycle': truss_cycle(), 'shell_mass_translation': shell_mass(0),
                       'shell_mass_rotation': shell_mass(4)}.items():
        (output/(name+'.i')).write_text('# Independent interface fixture; no Abaqus reference.\n'+text)
    # Independent force-controlled 1D isotropic hardening, not the return-map implementation.
    H = (696e6-400e6)/.1
    expected = [dict(time=t, axial_stress=s, signed_plastic_stretch=p,
                     equivalent_plastic_stretch=q, right_u=s/200e9+p)
                for t, s, p, q in [(0.5, 450e6, 50e6/H, 50e6/H),
                                   (1, 0, 50e6/H, 50e6/H),
                                   (1.5, -500e6, 0, 100e6/H),
                                   (2, 0, 0, 100e6/H)]]
    (output/'independent-expectations.json').write_text(json.dumps(dict(
        steel_truss_cycle=expected,
        shell_mass_translation=dict(mass_ux=.314, mass_ry=-.0157),
        shell_mass_rotation=dict(mass_ux=-.000157, mass_ry=.000010466666666666667),
        limitations=['Prescribed rigid shell motion tests mass assembly, not dynamic free response.',
                     'Cyclic axial steel fixture is not the full aqueduct or CDP regression.']), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    build(parser.parse_args().output)
