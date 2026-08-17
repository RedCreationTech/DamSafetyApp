#!/usr/bin/env python3
"""Generate a linear-elastic MOOSE static-to-dynamic 2D dam case.

The generated dynamic input initializes displacement and out-of-plane strain
from the static Exodus result, retains gravity and hydrostatic pressure, adds
directional nodal mass, and applies the Abaqus acceleration amplitude for the
full dynamic step.
"""

import argparse
import json
import os
from pathlib import Path


DOF_VARIABLES = {1: 'disp_x', 2: 'disp_y'}


def _number(value):
    return f'{float(value):.16g}'


def _relative(path, base):
    return Path(os.path.relpath(Path(path).resolve(), Path(base).resolve())).as_posix()


def _find_step(report, field):
    matches = [step for step in report['steps'] if step.get(field)]
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one {field} step, found {len(matches)}')
    return matches[0]


def _resolve_name(names, raw_name):
    if raw_name in names:
        return raw_name
    prefix = f'{raw_name}__'
    matches = [name for name in names if name.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f'Cannot uniquely resolve set {raw_name!r}: {matches}')
    return matches[0]


def _material(report):
    used = {block['material'] for block in report['blocks'].values()
            if block['material'] and block['material'].upper() != 'RIGID'}
    if len(used) != 1:
        raise ValueError(f'Expected one physical material, found {sorted(used)}')
    requested = next(iter(used))
    matches = [value for name, value in report['materials'].items()
               if name.casefold() == requested.casefold()]
    if len(matches) != 1:
        raise ValueError(f'Cannot resolve material {requested!r}')
    material = matches[0]
    try:
        youngs, poisson = material['elastic'][0][:2]
        density = material['density'][0][0]
    except (KeyError, IndexError) as error:
        raise ValueError('Material must define Elastic(E, nu) and Density') from error
    damping = material.get('damping', {})
    return {
        'youngs': youngs,
        'poisson': poisson,
        'density': density,
        'mass_damping': damping.get('alpha', 0.0),
        'stiffness_damping': damping.get('beta', 0.0),
    }


def _loads(report, static_step):
    gravity = [load for load in static_step['loads']
               if load['type'].upper() == 'GRAV']
    pressure = [load for load in static_step['loads']
                if load['type'].upper() == 'HP']
    if len(gravity) != 1 or len(pressure) != 1:
        raise ValueError('Expected one GRAV and one HP load in the static step')
    grav = gravity[0]
    direction = grav.get('parameters', [])
    if len(direction) < 2:
        raise ValueError('GRAV must provide at least x/y direction components')
    hp = pressure[0]
    parameters = hp.get('parameters', [])
    if len(parameters) < 1 or parameters[0] == 0:
        raise ValueError('HP must provide a nonzero zero-pressure y coordinate')
    return {
        'gravity_x': grav['value'] * direction[0],
        'gravity_y': grav['value'] * direction[1],
        'pressure_surface': _resolve_name(report['sidesets'], hp['surface']),
        'pressure_value': hp['value'],
        'water_level': parameters[0],
    }


def _fixed_bcs(report, definitions):
    blocks = []
    for index, boundary in enumerate(definitions, start=1):
        if boundary.get('amplitude'):
            continue
        if boundary.get('value', 0.0) != 0.0:
            raise ValueError('Only zero-valued Dirichlet boundaries are supported')
        boundary_name = _resolve_name(report['nodesets'], boundary['set'])
        for dof in range(boundary['dof1'], boundary['dof2'] + 1):
            variable = DOF_VARIABLES.get(dof)
            if not variable:
                continue
            blocks.append((f'fix_{index}_{dof}', variable, boundary_name))
    return blocks


def _bc_text(blocks):
    return '\n'.join(
        f'''  [{name}]
    type = DirichletBC
    variable = {variable}
    boundary = {boundary}
    value = 0
  []'''
        for name, variable, boundary in blocks)


def _load_text(loads):
    gravity = []
    for axis, value in (('x', loads['gravity_x']), ('y', loads['gravity_y'])):
        if value:
            gravity.append(f'''  [gravity_{axis}]
    type = Gravity
    variable = disp_{axis}
    value = {_number(value)}
  []''')
    gravity_text = '\n'.join(gravity)
    pressure_text = f'''  [water_pressure_x]
    type = Pressure
    variable = disp_x
    boundary = {loads['pressure_surface']}
    function = hydrostatic_pressure
  []
  [water_pressure_y]
    type = Pressure
    variable = disp_y
    boundary = {loads['pressure_surface']}
    function = hydrostatic_pressure
  []'''
    return gravity_text, pressure_text


def _static_input(mesh, material, loads, bcs, prefix):
    gravity, pressure = _load_text(loads)
    expression = (f"max(0, {_number(loads['pressure_value'])} * "
                  f"({_number(loads['water_level'])} - y) / "
                  f"{_number(loads['water_level'])})")
    return f'''# Generated linear-elastic reinterpretation of the Abaqus static step.
[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = {mesh}
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [strain_zz]
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [dam]
    planar_formulation = WEAK_PLANE_STRESS
    strain = SMALL
    generate_output = 'stress_xx stress_yy stress_xy vonmises_stress'
  []
[]

[Functions]
  [hydrostatic_pressure]
    type = ParsedFunction
    expression = '{expression}'
  []
[]

[Kernels]
{gravity}
[]

[BCs]
{_bc_text(bcs)}
{pressure}
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = {_number(material['youngs'])}
    poissons_ratio = {_number(material['poisson'])}
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = {_number(material['density'])}
  []
[]

[Postprocessors]
  [max_abs_disp_x]
    type = NodalExtremeValue
    variable = disp_x
    value_type = max_abs
  []
  [max_abs_disp_y]
    type = NodalExtremeValue
    variable = disp_y
    value_type = max_abs
  []
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
[]

[Executioner]
  type = Steady
  solve_type = LINEAR
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  l_tol = 1e-10
  l_max_its = 300
[]

[Outputs]
  exodus = true
  csv = true
  file_base = results/{prefix}_static
[]
'''


def _dynamic_input(static_result, accel_file, added_mass_x, added_mass_y,
                   material, loads, fixed_bcs, accel_bc, static_duration,
                   dynamic_duration, dt, output_interval, prefix):
    gravity, pressure = _load_text(loads)
    expression = (f"max(0, {_number(loads['pressure_value'])} * "
                  f"({_number(loads['water_level'])} - y) / "
                  f"{_number(loads['water_level'])})")
    end_time = static_duration + dynamic_duration
    return f'''# Generated full linear-elastic dynamic step initialized from the static result.
# Known limitation: Abaqus CDP and the preceding Lanczos frequency step are not migrated.
[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = {static_result}
    use_for_exodus_restart = true
  []
[]

[Variables]
  [disp_x]
    initial_from_file_var = disp_x
  []
  [disp_y]
    initial_from_file_var = disp_y
  []
  [strain_zz]
    initial_from_file_var = strain_zz
  []
[]

[AuxVariables]
  [vel_x]
  []
  [accel_x]
  []
  [vel_y]
  []
  [accel_y]
  []
[]

[Physics/SolidMechanics/Dynamic]
  [dam]
    add_variables = false
    velocities = 'vel_x vel_y'
    accelerations = 'accel_x accel_y'
    newmark_beta = 0.25
    newmark_gamma = 0.5
    hht_alpha = 0.0
    mass_damping_coefficient = {_number(material['mass_damping'])}
    stiffness_damping_coefficient = {_number(material['stiffness_damping'])}
    strain = SMALL
    incremental = false
    planar_formulation = WEAK_PLANE_STRESS
    density = density
    generate_output = 'stress_xx stress_yy stress_xy vonmises_stress'
  []
[]

[Functions]
  [hydrostatic_pressure]
    type = ParsedFunction
    expression = '{expression}'
  []
  [base_acceleration]
    type = PiecewiseLinear
    data_file = {accel_file}
    format = columns
  []
[]

[Kernels]
{gravity}
[]

[NodalKernels]
  [added_mass_x]
    type = NodalTranslationalInertia
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.25
    gamma = 0.5
    boundary = POINT_MASS
    nodal_mass_file = {added_mass_x}
  []
  [added_mass_y]
    type = NodalTranslationalInertia
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.25
    gamma = 0.5
    boundary = POINT_MASS
    nodal_mass_file = {added_mass_y}
  []
[]

[BCs]
{_bc_text(fixed_bcs)}
  [base_acceleration_x]
    type = PresetAcceleration
    variable = disp_x
    boundary = {accel_bc}
    function = base_acceleration
    beta = 0.25
    velocity = vel_x
    acceleration = accel_x
  []
{pressure}
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = {_number(material['youngs'])}
    poissons_ratio = {_number(material['poisson'])}
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = {_number(material['density'])}
  []
[]

[Postprocessors]
  [max_abs_disp_x]
    type = NodalExtremeValue
    variable = disp_x
    value_type = max_abs
  []
  [max_abs_disp_y]
    type = NodalExtremeValue
    variable = disp_y
    value_type = max_abs
  []
  [max_abs_accel_x]
    type = NodalExtremeValue
    variable = accel_x
    value_type = max_abs
  []
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
[]

[Executioner]
  type = Transient
  start_time = {_number(static_duration)}
  end_time = {_number(end_time)}
  dt = {_number(dt)}
  solve_type = LINEAR
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  l_tol = 1e-9
  l_max_its = 300
[]

[Outputs]
  exodus = true
  csv = true
  time_step_interval = {output_interval}
  file_base = results/{prefix}_dynamic
[]
'''


def generate_case(report_path, mesh_path, added_mass_x, added_mass_y,
                  output_dir, prefix='dam_static_then_dynamic',
                  output_interval=10):
    report_path = Path(report_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / 'results').mkdir(exist_ok=True)
    report = json.loads(report_path.read_text(encoding='utf-8'))

    static_step = _find_step(report, 'static')
    dynamic_step = _find_step(report, 'dynamic')
    material = _material(report)
    loads = _loads(report, static_step)
    static_bcs = _fixed_bcs(report, report.get('initial_boundaries', []))
    if not static_bcs:
        raise ValueError('No model-level static boundary conditions were converted')

    accel_boundaries = [bc for bc in dynamic_step['boundaries']
                        if bc.get('amplitude')]
    if len(accel_boundaries) != 1:
        raise ValueError('Expected exactly one amplitude-driven dynamic boundary')
    accel_boundary = accel_boundaries[0]
    if str(accel_boundary.get('type', '')).upper() != 'ACCELERATION':
        raise ValueError('Amplitude-driven boundary must be type=ACCELERATION')
    if accel_boundary['dof1'] != 1 or accel_boundary['dof2'] != 1:
        raise ValueError('This generator currently supports x acceleration only')
    accel_name = accel_boundary['amplitude']
    amplitude = report['amplitudes'][accel_name]
    amplitude_time = str(report.get('amplitude_options', {})
                         .get(accel_name, {}).get('time', 'STEP TIME')).upper()

    static_duration = float(static_step['static'][1])
    dt, dynamic_duration = map(float, dynamic_step['dynamic'][:2])
    time_shift = 0.0 if amplitude_time == 'TOTAL TIME' else static_duration
    acceleration_path = output_dir / f'{prefix}_acceleration.csv'
    acceleration_path.write_text(''.join(
        f'{_number(time + time_shift)},{_number(value * accel_boundary["value"])}\n'
        for time, value in amplitude), encoding='utf-8')

    fixed_dynamic = _fixed_bcs(report, dynamic_step['boundaries'])
    accel_set = _resolve_name(report['nodesets'], accel_boundary['set'])
    mesh_relative = _relative(mesh_path, output_dir)
    mass_x_relative = _relative(added_mass_x, output_dir)
    mass_y_relative = _relative(added_mass_y, output_dir)
    acceleration_relative = acceleration_path.name

    static_path = output_dir / f'{prefix}_static.i'
    dynamic_path = output_dir / f'{prefix}_dynamic.i'
    static_path.write_text(_static_input(
        mesh_relative, material, loads, static_bcs, prefix), encoding='utf-8')
    dynamic_path.write_text(_dynamic_input(
        f'results/{prefix}_static.e', acceleration_relative,
        mass_x_relative, mass_y_relative, material, loads, fixed_dynamic,
        accel_set, static_duration, dynamic_duration, dt,
        output_interval, prefix), encoding='utf-8')
    return static_path, dynamic_path, acceleration_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', required=True)
    parser.add_argument('--mesh', required=True)
    parser.add_argument('--added-mass-x', required=True)
    parser.add_argument('--added-mass-y', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--prefix', default='dam_static_then_dynamic')
    parser.add_argument('--output-interval', type=int, default=10)
    args = parser.parse_args()
    paths = generate_case(
        args.report, args.mesh, args.added_mass_x, args.added_mass_y,
        args.output_dir, args.prefix, args.output_interval)
    for path in paths:
        print(path)


if __name__ == '__main__':
    main()
