"""Strict HEX8 host containment for an audited embedded-element region.

Coordinates remain unchanged. No nearest-node fallback or snapping is used.
The map is preparation evidence; MOOSE evaluates the host shape functions in
EqualValueEmbeddedConstraint rather than reading these weights as penalties.
"""
import collections
import csv

import numpy as np
from scipy.spatial import cKDTree

import abaqus2exodus as mesh


def inverse_hex(point, corners, physical_tolerance=1e-10):
    """Return trilinear weights only for a genuinely contained point (SI units)."""
    natural = np.zeros(3)
    for _ in range(30):
        mapped, jacobian_rows = mesh._hex_map(corners, *natural)
        residual = point - np.asarray(mapped)
        if np.max(np.abs(residual)) <= physical_tolerance:
            if np.max(np.abs(natural)) > 1 + 1e-8:
                return None
            weights = np.asarray(mesh._hex_shape(*natural))
            if abs(weights.sum() - 1) > 1e-12:
                raise ValueError('HEX8 interpolation does not reproduce a constant field')
            if np.max(np.abs(weights @ corners - point)) > physical_tolerance:
                raise ValueError('HEX8 interpolation does not reproduce an affine field')
            return weights
        try:
            # _hex_map stores derivatives in ROWS; the physical mapping
            # Jacobian used by Newton has derivatives in COLUMNS.
            increment = np.linalg.solve(np.asarray(jacobian_rows).T, residual)
        except np.linalg.LinAlgError:
            return None
        natural += increment
        if not np.all(np.isfinite(natural)) or np.max(np.abs(natural)) > 10:
            return None
    return None


def prepare_embedded_hosts(model, gm, blocks, types, output):
    if len(model.embedded) != 1:
        raise ValueError('Expected one audited embedded region')
    region = model.embedded[0]
    host_origins = {(inst, eid) for inst, eids in
                    model.asm_elsets[region['host elset']].items() for eid in eids}
    embedded_origins = {(inst, eid) for name in region['embedded']
                        for inst, eids in model.asm_elsets[name].items() for eid in eids}
    hosts = []
    nodes_to_blocks = collections.defaultdict(set)
    for block, elems in blocks.items():
        for eid, conn in elems:
            origin = gm.elem_origin[eid]
            if origin in host_origins:
                if types[block] != 'HEX8':
                    raise ValueError('Embedded host region contains a non-HEX8 element')
                corners = np.asarray([gm.coords[n - 1] for n in conn])
                hosts.append((eid, block, conn, corners))
            if origin in embedded_origins:
                if types[block] != 'TRUSS':
                    raise ValueError('Embedded region contains a non-TRUSS element')
                for node in conn:
                    nodes_to_blocks[node].add(block)
    if not hosts or not nodes_to_blocks:
        raise ValueError('Empty embedded/host region')
    lower = np.asarray([corners.min(axis=0) for _, _, _, corners in hosts])
    upper = np.asarray([corners.max(axis=0) for _, _, _, corners in hosts])
    centers = (lower + upper) / 2
    search_radius = np.max(np.linalg.norm(upper - lower, axis=1)) / 2 + 1e-10
    tree = cKDTree(centers)
    pairs = set()
    failures = []
    cross_instance_ambiguities = []
    max_affine_error = 0.0
    multiplicities = collections.Counter()
    rows = []
    for node, secondary_blocks in sorted(nodes_to_blocks.items()):
        point = np.asarray(gm.coords[node - 1])
        contained = []
        for index in tree.query_ball_point(point, search_radius):
            if np.any(point < lower[index] - 1e-10) or np.any(point > upper[index] + 1e-10):
                continue
            eid, primary_block, conn, corners = hosts[index]
            weights = inverse_hex(point, corners)
            if weights is not None:
                contained.append((eid, primary_block, conn, weights))
                max_affine_error = max(max_affine_error,
                                       float(np.max(np.abs(weights @ corners - point))))
        if not contained:
            failures.append(node)
            continue
        multiplicities[len(contained)] += 1
        instances = {gm.elem_origin[eid][0] for eid, _, _, _ in contained}
        if len(instances) > 1:
            cross_instance_ambiguities.append({'node': node, 'instances': sorted(instances)})
        for _, primary_block, _, _ in contained:
            pairs.update((primary_block, secondary_block) for secondary_block in secondary_blocks)
        # Shared faces inside one instance have identical interpolated values;
        # keep one deterministic representative for the preparation CSV.
        eid, primary_block, conn, weights = min(contained, key=lambda row: row[0])
        rows.append((node, eid, primary_block, *conn, *weights))
    if failures or cross_instance_ambiguities:
        raise ValueError(f'Embedded host mapping unresolved: {len(failures)} outside nodes '
                         f'{failures[:8]}, {len(cross_instance_ambiguities)} cross-instance '
                         f'ambiguities {cross_instance_ambiguities[:3]}')
    with (output / 'embedded-host-map.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('secondary_exodus_node_number', 'host_exodus_element_label', 'host_block',
                         *(f'host_node_{i}' for i in range(1, 9)),
                         *(f'shape_weight_{i}' for i in range(1, 9))))
        writer.writerows(rows)
    # A pair constrains only secondary nodes found inside that primary block.
    # The complete, strict pre-audit above prevents silently uncovered steel.
    lines = ['# Embedded fragment only. Requires the full displacement problem.',
             '# Kinematic formulation transfers the original steel residual to the host.',
             '# penalty is a numerical constraint scale in N/m, not a bond stiffness.',
             '[Constraints]']
    for index, (primary, secondary) in enumerate(sorted(pairs)):
        for axis in 'xyz':
            lines.extend([f'  [embedded_{index}_{axis}]',
                          '    type = EqualValueEmbeddedConstraint',
                          f'    variable = disp_{axis}', f'    primary_variable = disp_{axis}',
                          f"    primary = '{primary}'", f"    secondary = '{secondary}'",
                          '    formulation = kinematic', '    penalty = 1e10', '  []'])
    lines.extend(['[]', ''])
    (output / 'embedded-constraint-fragment.i').write_text('\n'.join(lines))
    return {'embedded_element_count': len(embedded_origins), 'host_element_count': len(hosts),
            'embedded_node_count': len(nodes_to_blocks), 'outside_nodes': [],
            'cross_instance_ambiguities': [], 'nearest_node_fallbacks': 0,
            'max_affine_reconstruction_error_m': max_affine_error,
            'host_multiplicity_histogram': dict(sorted(multiplicities.items())),
            'constraint_pairs': sorted(pairs), 'constraint_count': len(pairs) * 3,
            'runtime_verification': 'pending; fragment is not a standalone input'}
