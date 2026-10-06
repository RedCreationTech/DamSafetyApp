#!/usr/bin/env python3
"""Prepare audited geometry/material data, without claiming a runnable wall input.

This case-specific tool preserves all instance identities and only accepts the
reviewed mm-N-s input. It does not generate a submission manifest or run MOOSE.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

import abaqus2exodus as mesh
import abaqus_cdp as cdp
from embedded_host_map import prepare_embedded_hosts


SOURCE_SHA = '21299c9694ecb4597c84e1259050d645b5299ddf72d63e9647782b9c37327d25'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def prepare(source, output):
    if sha(source) != SOURCE_SHA:
        raise ValueError('Source hash changed; review the new INP before conversion')
    # Prevent accidental replacement of an earlier input/evidence snapshot.
    output.mkdir(parents=True, exist_ok=False)
    model = mesh.parse_inp(source)
    if model.ties or model.mpcs:
        raise ValueError('Unexpected Tie/MPC in the reviewed input')
    gm, blocks, types, meta, nodesets, sidesets = mesh.build_global_mesh(
        model, 1e-9, preserve_instance_nodes=True)
    counts = collections.Counter()
    for block, elems in blocks.items():
        counts[types[block]] += len(elems)
    if len(gm.coords) != 24705 or dict(counts) != {'HEX8': 14900, 'TRUSS': 3646}:
        raise ValueError(f'Unexpected mesh inventory: {len(gm.coords)} / {counts}')
    if gm.merged_count != 0 or len(set(gm.node_map.values())) != len(gm.coords):
        raise ValueError('Node identities were merged')
    gm.coords = [tuple(value * 1e-3 for value in point) for point in gm.coords]
    if mesh.find_inverted_hexes(gm, blocks, types):
        raise ValueError('Non-positive HEX8 corner Jacobian')
    embedded_audit = prepare_embedded_hosts(model, gm, blocks, types, output)
    write_json(output / 'embedded-host-audit.json', embedded_audit)
    mesh.write_exodus(output / 'wall_mesh.e', gm, blocks, types, meta,
                       nodesets, sidesets, 'luwen-wxz geometry in metres; formulation pending')
    mesh.write_node_map_csv(gm, output / 'node-map.csv')
    block_names = sorted(blocks)
    with (output / 'element-map.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('exodus_internal_element_number', 'exodus_element_label',
                         'instance', 'abaqus_element_label', 'block', 'original_type'))
        index = 0
        for block in block_names:
            for eid, conn in blocks[block]:
                index += 1
                instance, original = gm.elem_origin[eid]
                part = model.parts[mesh.model_inst_part(model, instance)]
                writer.writerow((index, eid, instance, original, block, part.elem_types[original]))
    active_materials = {value[1] for value in meta.values()}
    material_data = {}
    for name in sorted(active_materials):
        data = model.materials[name]
        material_data[name] = {
            'youngs_modulus_pa': data['elastic'][0][0] * 1e6,
            'poissons_ratio': data['elastic'][0][1],
            'density_kg_m3': data['density'][0][0] * 1e12,
            'plastic_stress_pa_strain': [[row[0] * 1e6, row[1]]
                                        for row in data.get('plastic', [])]}
    aac = cdp.validate_material(cdp.parse_materials(source)['aac201-lunwen'])
    headers = {
        'compression_hardening': ('stress_pa', 'inelastic_strain'),
        'compression_damage': ('damage_c', 'inelastic_strain'),
        'tension_stiffening': ('stress_pa', 'cracking_strain'),
        'tension_damage': ('damage_t', 'cracking_strain')}
    for key, table in aac['tables'].items():
        scale = 1 if key.endswith('damage') else 1e6
        with (output / f'{key}.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(headers[key])
            writer.writerows((format(y * scale, '.17g'), format(x, '.17g')) for y, x in table)
    material_data['aac201-lunwen'].update({
        'cdp': aac['plasticity'], 'recovery': aac['recovery'],
        'table_rows': {key: len(table) for key, table in aac['tables'].items()}})
    write_json(output / 'materials-si.json', material_data)
    aac_blocks = ' '.join(name for name in block_names if meta[name][1] == 'aac201-lunwen')
    p = aac['plasticity']
    params = {
        'youngs_modulus': material_data['aac201-lunwen']['youngs_modulus_pa'],
        'poissons_ratio': material_data['aac201-lunwen']['poissons_ratio'],
        'dilation_angle': p['dilation_angle_degrees'], 'eccentricity': p['eccentricity'],
        'biaxial_to_uniaxial_compression_ratio': p['fb0_fc0'],
        'tensile_meridian_ratio': p['kc'], 'viscosity': p['viscosity_seconds'],
        **aac['recovery']}
    lines = ['# Material fragment only; not a standalone executable input.',
             '[Materials]', '  [aac_cdp]', '    type = AbaqusCDPStressUpdate',
             f"    block = '{aac_blocks}'"]
    lines += [f'    {key}_file = {key}.csv' for key in headers]
    lines += [f'    {key} = {value:.17g}' for key, value in params.items()]
    lines += ['  []', '[]', '']
    (output / 'aac-material-fragment.i').write_text('\n'.join(lines))
    # Retain exact keyword/options/data and line provenance for the missing
    # contact implementation; these are source parameters, not MOOSE objects.
    cards = []
    for number, line in enumerate(source.read_bytes().decode('latin1').splitlines(), 1):
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            keyword, options = cdp._parse_keyword(line)
            cards.append({'keyword': keyword, 'options': options, 'line': number, 'rows': []})
        else:
            cards[-1]['rows'].append(line)
    contact_keywords = {'surface interaction', 'friction', 'surface behavior',
                        'cohesive behavior', 'damage initiation', 'damage evolution',
                        'damage stabilization', 'contact pair'}
    contact_cards = [card for card in cards if card['keyword'] in contact_keywords]
    pairs = [card for card in contact_cards if card['keyword'] == 'contact pair']
    if len(pairs) != 7:
        raise ValueError('Expected seven contact pairs')
    for pair in pairs:
        pair['mapped_surfaces'] = []
        for name in pair['rows'][0].split(','):
            name = name.strip()
            converted = mesh.sanitize(name)
            entries = sidesets[converted]
            if not entries:
                raise ValueError(f'Empty contact surface {name}')
            pair['mapped_surfaces'].append({'abaqus': name, 'exodus_sideset': converted,
                                           'face_count': len(entries)})
        a, b = [mesh.sanitize('SURF_' + name.strip()) for name in pair['rows'][0].split(',')]
        if set(nodesets[a]) & set(nodesets[b]):
            raise ValueError('Contact surfaces share node identities')
    write_json(output / 'contact-source-and-map.json', {
        'original_cards': contact_cards,
        'si_parameters': {'friction_coefficient': 0.6, 'slip_tolerance': 0.005,
                          'normal_pressure_overclosure': 'HARD',
                          'cohesive_stiffness_pa_per_m': [2.559e12, 1.093e12, 1.093e12],
                          'damage_initiation': 'MAXS', 'strength_pa': [5.63e6] * 3,
                          'damage_evolution': 'DISPLACEMENT',
                          'evolution_displacement_m': 0.000241,
                          'damage_stabilization_s': 0.0005},
        'implementation_status': 'pending; do not substitute BK/POWER_LAW or a scalar stiffness'})
    write_json(output / 'geometry-constraint-map.json', {
        'couplings_source': model.couplings, 'embedded_source': model.embedded,
        'initial_boundaries_source_units': model.initial_boundaries,
        'steps_source_units': model.steps, 'amplitudes_source': model.amplitudes,
        'rp_nodes': gm.asm_node_map,
        'sections_original_mm2': {part.name: part.sections for part in model.parts.values()},
        'blocks': {b: {'id': i, 'material': meta[b][1], 'count': len(blocks[b]), 'type': types[b]}
                   for i, b in enumerate(block_names, 1)},
        'nodesets': {name: len(ids) for name, ids in sorted(nodesets.items())},
        'sidesets': {name: len(entries) for name, entries in sorted(sidesets.items())},
        'constraint_implementation_status': 'embedded-fragment-prepared-RP-runtime-pending',
        'embedded_host_audit': embedded_audit})
    files = [{'path': path.name, 'sha256': sha(path), 'size': path.stat().st_size}
             for path in sorted(output.iterdir())]
    manifest = {'status': 'geometry-and-material-data-prepared-not-runnable',
                'source': str(source), 'source_sha256': SOURCE_SHA,
                'generator_sha256': sha(__file__), 'mesh_converter_sha256': sha(mesh.__file__),
                'coordinate_units': 'm', 'material_units': 'Pa-kg-m-s',
                'nodes': len(gm.coords), 'elements': dict(counts),
                'merged_nodes': gm.merged_count, 'contact_pairs': len(pairs),
                'pending': ['C3D8R reduced integration/hourglass formulation',
                            'original MAXS/displacement cohesive contact implementation',
                            'embedded runtime verification and RP coupling',
                            'full input/output/executioner and LIMS submission'],
                'files': files}
    write_json(output / 'preparation-manifest.json', manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inp', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = prepare(args.inp.resolve(), args.output.resolve())
    print(json.dumps({key: result[key] for key in ('status', 'nodes', 'elements', 'contact_pairs')}))
