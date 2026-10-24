#!/usr/bin/env python3
"""Prepare the approved AQ1023 full assembly without dropping interfaces.

This creates geometry and a source-semantic contract, not a runnable .i or Job.
All instances retain separate nodes, even at coincident coordinates. Contact,
Tie, embedded, shell-solid and connector cards remain explicit pending physics.
"""
import argparse
import copy
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

import abaqus2exodus as mesh
import abaqus_cdp as cdp


SOURCE_SHA = '4d6ce44e5f66ddc6abce3a3cae506402653bef660d43aa94738e78d273f2a95a'
KNOWN = set('''heading|preprint|part|node|element|nset|elset|solid section|
end part|shell section|assembly|instance|end instance|connector section|
surface|tie|embedded element|shell to solid coupling|end assembly|
connector behavior|connector elasticity|amplitude|baseline correction|
material|density|elastic|concrete damaged plasticity|concrete compression hardening|
concrete tension stiffening|concrete compression damage|concrete tension damage|
plastic|acoustic medium|surface interaction|friction|surface behavior|
impedance property|boundary|contact pair|step|static|restart|output|
node output|element output|end step|dynamic'''.replace('\n', '').split('|'))
PENDING = [
    'S4R shell thickness offset plasticity and dynamic rotation',
    'shell-to-solid force and moment transfer',
    'Cartesian connector translation and Align rotation constraints',
    'host-restricted embedded interpolation without nearest-node fallback',
    'nonmatching Tie interpolation and adjust semantics',
    'hard contact friction initial gap and slip regularization',
    'steel solid and truss plasticity and section inertia',
    'baseline-corrected three-axis motion and time-integration initial state',
    'static stabilization and complete static-to-dynamic history',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_confirmation(approval, source_sha):
    if (source_sha != SOURCE_SHA or approval.get('source_sha256') != SOURCE_SHA
            or approval.get('case_id') != 'AQ1023-FULL-ELC'):
        raise ValueError('source does not match the expert-approved AQ1023 input')
    a1, a2 = approval.get('A01', {}), approval.get('A02', {})
    units = dict(length='mm', force='N', time='s', stress='MPa', mass='tonne')
    if (a1.get('status') != 'confirmed' or a2.get('status') != 'confirmed'
            or a1.get('authoritative_complete_input') is not True
            or a1.get('external_overrides') is not False
            or a1.get('add_gravity_or_water') is not False
            or a2.get('source_units') != units
            or a2.get('preserve_materials_sections_and_connector_parameters') is not True
            or a2.get('cdp_viscosity_seconds') != 0.05 or a2.get('output_units') != 'SI'):
        raise ValueError('confirmation must preserve the approved complete source and mm-N-s-MPa-tonne units')


def cards_from_text(text):
    """Preserve every keyword/data block with its original line and scope."""
    cards = []
    part = step = material = None
    for number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            name = line[1:].split(',')[0].strip().lower()
            if name not in KNOWN:
                raise ValueError(f'line {number}: unknown keyword *{name}')
            options = mesh._parse_kv(line[1 + len(line[1:].split(',')[0]):])
            if name == 'part':
                part = options['name']
            if name == 'step':
                step = options['name']
            if name == 'material':
                material = options['name']
            cards.append(dict(keyword=name, options=options, line=number,
                              part=part, step=step, material=material, data=[]))
            if name == 'end part':
                part = None
            if name == 'end step':
                step = None
        else:
            if not cards:
                raise ValueError(f'line {number}: data before keyword')
            cards[-1]['data'].append(line)
    return cards


def numbers(row):
    values = [float(x.strip()) for x in row.rstrip(',').split(',') if x.strip()]
    if not all(math.isfinite(x) for x in values):
        raise ValueError('nonfinite numeric data')
    return values


def validate_rows(cards):
    """The legacy parser tolerates malformed rows; this entry must not."""
    seen = set()
    for c in cards:
        kind, part = c['keyword'], c['part']
        if kind in ('node', 'element') and part:
            width = 4 if kind == 'node' else {
                'C3D8R': 9, 'C3D8': 9, 'S4R': 5, 'T3D2': 3,
            }.get(c['options'].get('type', '').upper())
            if width is None:
                raise ValueError(f'unsupported Part element at line {c["line"]}')
            for row in c['data']:
                values = numbers(row)
                if len(values) != width or values[0] != int(values[0]):
                    raise ValueError(f'malformed {kind} row near line {c["line"]}')
                key = (part.upper(), kind, int(values[0]))
                if key in seen:
                    raise ValueError(f'duplicate identity {key}')
                seen.add(key)
                if kind == 'element' and any(x != int(x) for x in values):
                    raise ValueError('noninteger connectivity')
        elif kind == 'amplitude':
            values = [v for row in c['data'] for v in numbers(row)]
            if len(values) % 2 or not values:
                raise ValueError('amplitude requires complete time/value pairs')
            times = values[::2]
            if any(b <= a for a, b in zip(times, times[1:])):
                raise ValueError('amplitude time is not strictly increasing')


def resolve_elset(model, name):
    """Expand source elsets to (instance, element) identities, never nodes."""
    groups = model.asm_elsets.get(name)
    if not isinstance(groups, dict) or not groups:
        raise ValueError(f'missing assembly element set: {name}')
    instances = {i.name: i.part for i in model.instances}
    result = []
    for inst, labels in groups.items():
        part = model.parts[instances[inst]]
        for label in labels:
            if label not in part.elems:
                raise ValueError(f'unknown element {inst}.{label} in {name}')
            result.append((inst, label))
    return sorted(set(result))


def resolve_surfaces(model, gm):
    faces = dict(mesh.C3D8_FACES)
    # Legacy mesh.C3D8_FACES is an unordered node-set map. Numerical QUAD4
    # integration requires a cyclic S6 order; (1,4,5,8) is a bow tie.
    faces['S6'] = (4,1,5,8)
    shell = {'SPOS': (1, 2, 3, 4), 'SNEG': (1, 2, 3, 4),
             'E1': (1, 2), 'E2': (2, 3), 'E3': (3, 4), 'E4': (4, 1)}
    instances = {i.name: i.part for i in model.instances}
    result = {}
    for name, entries in model.surfaces.items():
        facets = []
        for es, face in entries:
            for inst, label in resolve_elset(model, es):
                p = model.parts[instances[inst]]
                typ = p.elem_types[label]
                local = (shell if typ == 'S4R' else faces).get(face)
                if local is None:
                    raise ValueError(f'unsupported face {typ}/{face}')
                facets.append(dict(instance=inst, element_label=label,
                                   face=face, element_type=typ,
                                   nodes=[gm.node_map[(inst, p.elems[label][n-1])]
                                          for n in local]))
        if not facets:
            raise ValueError(f'empty surface: {name}')
        result[name] = facets
    return result


def section_contract(model, cards):
    result = []
    for c in cards:
        if c['keyword'] not in ('solid section', 'shell section'):
            continue
        o = c['options']
        p = model.parts[c['part']]
        labels = sorted(set(p.elsets.get(o['elset'], [])))
        if not labels:
            raise ValueError(f'empty section elset {o["elset"]}')
        vals = numbers(c['data'][0]) if c['data'] else []
        entry = dict(part=c['part'], elset=o['elset'], material=o['material'],
                     keyword=c['keyword'], source_line=c['line'],
                     element_labels=labels)
        if c['keyword'] == 'shell section':
            if len(vals) != 2 or vals[0] <= 0 or vals[1] != int(vals[1]):
                raise ValueError('shell thickness/integration definition required')
            entry.update(thickness_m=vals[0] * 1e-3,
                         thickness_integration_points=int(vals[1]),
                         offset=o.get('offset', '0'))
            p.sections.append((o['elset'], o['material'], None))
        elif any(p.elem_types[e] == 'T3D2' for e in labels):
            if len(vals) != 1 or vals[0] <= 0:
                raise ValueError('truss cross-sectional area required')
            entry['area_m2'] = vals[0] * 1e-6
        result.append(entry)
    # Every element must have exactly one section; overlap can hide material loss.
    covered = Counter((s['part'], e) for s in result for e in s['element_labels'])
    for name, p in model.parts.items():
        for eid in p.elems:
            if covered[(name, eid)] != 1:
                raise ValueError(f'element needs one section: {name}/{eid}')
    return result


def interface_contract(model, cards, gm, surfaces):
    result = {'ties': [], 'embedded': [], 'shell_solid': [],
              'connectors': [], 'connector_sections': [], 'contacts': [],
              'source_cards': [], 'connector_behaviors_si': {}}
    group = {'tie': 'ties', 'embedded element': 'embedded',
             'shell to solid coupling': 'shell_solid', 'contact pair': 'contacts',
             'connector section': 'connector_sections'}
    retained = set(group) | {'connector behavior', 'connector elasticity',
                            'surface interaction', 'friction', 'surface behavior',
                            'impedance property', 'baseline correction',
                            'static', 'dynamic', 'boundary', 'output',
                            'node output', 'element output', 'acoustic medium'}
    connector_ids = set()
    behavior = None
    for c in cards:
        key = c['keyword']
        if key in retained:
            result['source_cards'].append(c)
        if key == 'connector behavior':
            behavior = c['options']['name']
            if behavior in result['connector_behaviors_si']:
                raise ValueError(f'duplicate connector behavior {behavior}')
            result['connector_behaviors_si'][behavior] = {}
        elif key == 'connector elasticity':
            component = int(c['options']['component'])
            if behavior is None or component not in (1, 2, 3):
                raise ValueError('unsupported connector elasticity component')
            if set(c['options']) != {'component'} or len(c['data']) != 1:
                raise ValueError('only uncoupled linear connector elasticity is supported')
            values = numbers(c['data'][0])
            if len(values) != 1 or values[0] <= 0:
                raise ValueError('positive connector stiffness required')
            stiffness = result['connector_behaviors_si'][behavior]
            if str(component) in stiffness:
                raise ValueError('duplicate connector stiffness component')
            stiffness[str(component)] = values[0] * 1000
        if key == 'element' and not c['part']:
            if c['options'].get('type') != 'CONN3D2':
                raise ValueError('unsupported assembly element')
            for row in c['data']:
                tokens = [t.strip() for t in row.split(',')]
                if len(tokens) != 3:
                    raise ValueError('connector requires two qualified endpoints')
                label = int(tokens[0])
                if label in connector_ids:
                    raise ValueError(f'duplicate connector {label}')
                connector_ids.add(label)
                ends = []
                for tok in tokens[1:]:
                    instance, node = tok.rsplit('.', 1)
                    if (instance, int(node)) not in gm.node_map:
                        raise ValueError(f'unknown connector endpoint: {tok}')
                    gid = gm.node_map[(instance, int(node))]
                    ends.append(dict(instance=instance, node_label=int(node), node=gid))
                result['connectors'].append(dict(label=label, endpoints=ends,
                                                source_line=c['line']))
        if key in ('tie', 'shell to solid coupling', 'contact pair'):
            if len(c['data']) != 1:
                raise ValueError(f'{key} needs a surface pair')
            pair = [x.strip() for x in c['data'][0].split(',')]
            if len(pair) != 2 or any(x not in surfaces for x in pair):
                raise ValueError(f'unresolved surface pair: {c}')
            result[group[key]].append(dict(options=c['options'],
                                           surfaces=pair, source_line=c['line']))
        elif key == 'embedded element':
            host = resolve_elset(model, c['options']['host elset'])
            children = [t.strip() for row in c['data'] for t in row.split(',') if t.strip()]
            entries = [e for name in children for e in resolve_elset(model, name)]
            result['embedded'].append(dict(options=c['options'], host=host,
                                           embedded=entries, source_line=c['line']))
        elif key == 'connector section':
            labels = model.asm_elsets.get(c['options']['elset'])
            if not isinstance(labels, list) or not labels:
                raise ValueError('connector section must resolve assembly labels')
            if any(label not in connector_ids for label in labels):
                raise ValueError('connector section references unknown connector')
            result['connector_sections'].append(dict(options=c['options'], labels=labels,
                                                     data=c['data'], source_line=c['line']))
    assigned = Counter()
    for section in result['connector_sections']:
        name = section['options']['behavior']
        if name not in result['connector_behaviors_si']:
            raise ValueError(f'unknown connector behavior {name}')
        assigned.update(section['labels'])
    if set(assigned) != connector_ids or any(n != 1 for n in assigned.values()):
        raise ValueError('every connector needs exactly one section')
    interactions = {c['options']['name'] for c in cards
                    if c['keyword'] == 'surface interaction'}
    if any(c['options']['interaction'] not in interactions for c in result['contacts']):
        raise ValueError('contact references missing interaction')
    return result


def material_si(model):
    result = {}
    stresses = {'plastic', 'cdp_compression_hardening', 'cdp_tension_stiffening'}
    for name, fields in model.materials.items():
        dst = json.loads(json.dumps(fields))
        for key, rows in dst.items():
            if key == 'density':
                for row in rows:
                    row[0] *= 1e12
            elif key == 'elastic' or key in stresses:
                for row in rows:
                    row[0] *= 1e6
        result[name] = dst
    return result


def prepare(source, output, *, encoding='gbk'):
    source, output = Path(source), Path(output)
    text = source.read_bytes().decode(encoding, errors='strict')
    cards = cards_from_text(text)
    validate_rows(cards)
    if output.exists():
        raise ValueError('output already exists; use a new directory')
    with tempfile.TemporaryDirectory() as tmp:
        normalized = Path(tmp) / 'source.inp'
        normalized.write_text(text, encoding='utf-8')
        model = mesh.parse_inp(normalized)
        sections = section_contract(model, cards)
        materials = cdp.parse_materials(normalized)
        for material in materials.values():
            cdp.validate_material(material)
    # Do not mutate the global supported types or enable S4R in the legacy CLI.
    type_map = dict(mesh.ELEMENT_TYPE_MAP, S4R=('SHELL4', 4))
    # Assembly-only connector elsets have no Part connectivity. Retain them in
    # the interface contract, not in the legacy Part-elset-to-nodeset loop.
    geometry_model = copy.copy(model)
    geometry_model.asm_elsets = {k: v for k, v in model.asm_elsets.items()
                                if isinstance(v, dict)}
    gm, blocks, types, meta, nodesets, sidesets = mesh.build_global_mesh(
        geometry_model, 1e-9, merge_nodes=False, element_type_map=type_map)
    gm.coords = [tuple(v * 1e-3 for v in xyz) for xyz in gm.coords]
    surfaces = resolve_surfaces(model, gm)
    interfaces = interface_contract(model, cards, gm, surfaces)
    instances = {i.name: i.part for i in model.instances}
    counts = Counter(model.parts[instances[i]].elem_types[e] for i, e in gm.elem_origin.values())
    for name, (_, material) in meta.items():
        if material == 'UNASSIGNED':
            raise ValueError(f'geometry block lacks material: {name}')
    expected_nodes = sum(len(model.parts[i.part].nodes) for i in model.instances)
    if len(gm.coords) != expected_nodes or gm.merged_count:
        raise ValueError('instance nodes were lost or merged')
    report = dict(schema_version='aqueduct-assembly-preparation/1',
                  case_id='AQ1023-FULL-ELC', status='geometry-and-semantic-preparation',
                  ready_for_solve=False, pending_physics=PENDING,
                  source=dict(filename=source.name, sha256=sha(source), encoding=encoding),
                  units=dict(source='mm-N-s-MPa-tonne', geometry='m',
                             materials_si='Pa-kg-m-s', sections='m-m2',
                             connector_behaviors_si='N/m', source_cards='source units'),
                  parts=len(model.parts), instances=len(model.instances), nodes=len(gm.coords),
                  mesh_element_counts=dict(counts),
                  mesh_elements=sum(counts.values()),
                  interface_connectors=len(interfaces['connectors']),
                  merged_nodes=gm.merged_count, sections=sections,
                  materials_si=material_si(model), surfaces=surfaces, interfaces=interfaces,
                  amplitude_options_raw=model.amplitude_options,
                  amplitudes_raw=model.amplitudes,
                  raw_amplitudes_are_effective_motion=False,
                  initial_boundaries_raw=model.initial_boundaries, steps_raw=model.steps)
    output.mkdir(parents=True)
    (output / 'source-utf8.inp').write_text(text, encoding='utf-8')
    mesh.write_exodus(output / 'assembly-geometry.e', gm, blocks, types, meta,
                      nodesets, sidesets, 'AQ1023 geometry only; no interface physics')
    with (output / 'node-map.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['instance', 'abaqus_node_label', 'exodus_node', 'x_m', 'y_m', 'z_m'])
        for (inst, label), gid in gm.node_map.items():
            w.writerow([inst, label, gid, *gm.coords[gid-1]])
    with (output / 'element-map.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['instance', 'abaqus_element_label', 'exodus_element', 'block'])
        for eid, (inst, label) in gm.elem_origin.items():
            w.writerow([inst, label, eid, gm.elem_block[eid]])
    report['artifacts'] = {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output / 'assembly-contract.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    summary = {k: report[k] for k in ('schema_version', 'case_id', 'status', 'ready_for_solve',
                                     'source', 'parts', 'instances', 'nodes', 'mesh_element_counts',
                                     'mesh_elements', 'interface_connectors', 'merged_nodes', 'pending_physics')}
    summary['contract_sha256'] = sha(output / 'assembly-contract.json')
    summary['artifacts'] = report['artifacts']
    (output / 'preparation-summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inp', required=True, type=Path)
    ap.add_argument('--confirmation', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    approval = json.loads(args.confirmation.read_text())
    try:
        validate_confirmation(approval, sha(args.inp))
    except ValueError as e:
        ap.error(str(e))
    summary = prepare(args.inp, args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
