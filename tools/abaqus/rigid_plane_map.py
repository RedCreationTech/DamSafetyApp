"""Reference-face integration for free RP rotations; no imposed top rotation."""
import csv
import numpy as np
import abaqus2exodus as mesh


def surface_integrals(model, gm, surface, reference):
    weights = {}
    gauss, rule_weights = np.polynomial.legendre.leggauss(3)
    signs = np.asarray(((-1, -1), (1, -1), (1, 1), (-1, 1)))
    reference = np.asarray(reference)
    count = 0
    for elset, face in model.surfaces[surface]:
        for instance, ids in model.asm_elsets[elset].items():
            part = model.parts[mesh.model_inst_part(model, instance)]
            for eid in ids:
                nodes = [gm.node_map[(instance, part.elems[eid][i - 1])]
                         for i in mesh.C3D8_FACES[face]]
                points = np.asarray([gm.coords[node - 1] for node in nodes])
                for node in nodes:
                    weights.setdefault(node, np.zeros(4))
                for xi, wx in zip(gauss, rule_weights):
                    for eta, wy in zip(gauss, rule_weights):
                        shape = (1 + signs[:, 0] * xi) * (1 + signs[:, 1] * eta) / 4
                        dx = signs[:, 0] * (1 + signs[:, 1] * eta) / 4
                        dy = signs[:, 1] * (1 + signs[:, 0] * xi) / 4
                        measure = np.linalg.norm(np.cross(dx @ points, dy @ points)) * wx * wy
                        position = shape @ points - reference
                        for i, node in enumerate(nodes):
                            weights[node] += measure * shape[i] * np.r_[1, position]
                count += 1
    return weights, count


def prepare_rigid_plane(model, gm, output):
    if len(model.couplings) != 1:
        raise ValueError('Expected one original RP coupling')
    coupling = model.couplings[0]
    rp_ids = model.asm_nsets[coupling['ref node']]
    if len(rp_ids) != 1:
        raise ValueError('Expected one reference point')
    position = gm.coords[gm.asm_node_map[rp_ids[0]] - 1]
    weights, count = surface_integrals(model, gm, coupling['surface'], position)
    area = sum(row[0] for row in weights.values())
    moments = sum((row[1:] for row in weights.values()), np.zeros(3))
    if abs(area - 1.056) > 1e-10 or np.max(np.abs(moments)) > 1e-10:
        raise ValueError('Unexpected top geometry/RP reference position')
    with (output / 'rp-surface-integrals.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('moose_node_id', 'x', 'y', 'z', 'area_weight',
                         'moment_x', 'moment_y', 'moment_z'))
        for node, row in sorted(weights.items()):
            writer.writerow((node - 1, *gm.coords[node - 1], *row))
    # Included by the main input once the displacement variables/mesh exist.
    # Top-face displacement DOFs are shared by the HEX8 and lower face block.
    lines = ['# RP fragment; requires surface block 500 and disp_x/disp_y/disp_z.',
             '[Variables]']
    for axis in ('x', 'y', 'z'):
        lines += [f'  [rp_traction_{axis}]', '    family = LAGRANGE',
                  '    order = FIRST', '    block = 500', '  []']
    # Locked SystemBase::applyScalingFactors applies field factors first,
    # then scalar factors. Register all field variables before this scalar.
    lines += ['  [rp_rotations]', '    family = SCALAR', '    order = THIRD', '  []',
              '[]', '[Functions]', '  [rp_y_loading]', '    type = PiecewiseLinear',
              "    x = '0 1'", "    y = '0 -0.020'", '  []', '[]', '[Kernels]']
    for component, axis in enumerate(('x', 'y', 'z')):
        lines += [f'  [rp_kinematics_{axis}]', '    type = RigidPlaneConstraint',
                  f'    variable = rp_traction_{axis}', f'    displacement = disp_{axis}',
                  '    rotations = rp_rotations', f'    component = {component}',
                  "    reference_point = '" + ' '.join(format(v, '.17g') for v in position) + "'",
                  '    block = 500',
                  f"    translation = {'rp_y_loading' if axis == 'y' else '0'}", '  []',
                  f'  [rp_force_{axis}]', '    type = CoupledForce', f'    variable = disp_{axis}',
                  f'    v = rp_traction_{axis}', '    coef = -1', '    block = 500',
                  '    use_displaced_mesh = false', '  []']
    lines += ['[]', '[ScalarKernels]', '  [rp_free_moments]',
              '    type = RigidPlaneMomentConstraint', '    variable = rp_rotations',
              "    tractions = 'rp_traction_x rp_traction_y rp_traction_z'",
              f"    boundary = {mesh.sanitize('SURF_' + coupling['surface'])}",
              '    moment_file = rp-surface-integrals.csv', '  []', '[]', '']
    (output / 'rp-constraint-fragment.i').write_text('\n'.join(lines))
    return {'rp_position_m': position, 'surface': coupling['surface'], 'face_count': count,
            'node_count': len(weights), 'area_m2': area, 'first_moments_m3': list(moments),
            'prescribed_translation_m': [0, -0.020, 0], 'free_rotations': ['x', 'y', 'z'],
            'formulation': 'surface LM kinematics plus zero generalized RP moments',
            'runtime_verification': 'pending; no top translations replaced by fixed values',
            'fragment': 'rp-constraint-fragment.i', 'required_lower_face_block': 500}
