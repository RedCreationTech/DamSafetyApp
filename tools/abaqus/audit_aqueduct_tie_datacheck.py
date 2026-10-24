#!/usr/bin/env python3
"""Audit printed Abaqus Tie membership without guessing suppressed warning counts.

DAT global labels are independent of converter/Exodus labels. Keep printed main
node lists and DOF strings verbatim; they are not interpolation coefficients.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_dat(text):
    lines = text.splitlines()
    nodes = {}
    ties = {}
    adjustments = {}
    contacts = []
    node_map = tie_table = contact_table = False
    current_tie = previous = None
    for line_number, line in enumerate(lines, 1):
        if 'LOCAL TO GLOBAL NODE AND ELEMENT MAPS' in line:
            node_map = True
        if node_map and re.search(r'element\s+element', line):
            node_map = False
        if node_map:
            match = re.match(r'^\s*(\S+)\s+(\d+)\s+(\d+)\s*$', line)
            if match:
                instance, local, global_label = match.groups()
                label = int(global_label)
                if label in nodes:
                    raise ValueError('duplicate DAT global node label')
                nodes[label] = (instance.upper(), int(local))
        if 'T I E  C O N S T R A I N T S' in line:
            tie_table = True
        if tie_table and 'E M B E D D E D' in line:
            tie_table = False
        if tie_table:
            match = re.match(r'^\s*(\d+)\s+(TRANSLATION|NBEAM-HINGE)\s+([\d\s]+)$', line)
            if match:
                label, dof, main = match.groups()
                previous = int(label)
                if previous in ties:
                    raise ValueError('duplicate secondary Tie label')
                ties[previous] = {'dof_printed': dof, 'main_global_labels': list(map(int, main.split())), 'dat_line': line_number}
            elif previous is not None and line[:28].strip() == '' and re.fullmatch(r'\s*\d[\d\s]*', line):
                ties[previous]['main_global_labels'].extend(map(int, line.split()))
        match = re.search(r'\*tie,\s*name=([^,\s]+)', line, re.I)
        if match:
            current_tie = match.group(1).upper()
        match = re.match(r'^\s*(\d+)\s+([\d.Ee+\-]+)\s+ADJUSTED TO\s*\(([^)]+)\)', line)
        if match and current_tie:
            label, distance, coordinates = match.groups()
            rows = adjustments.setdefault(current_tie, [])
            rows.append({'global_label': int(label), 'distance_mm_printed': float(distance),
                         'adjusted_xyz_mm_printed': [float(x) for x in coordinates.split(',')], 'dat_line': line_number})
        if 'C O N T A C T   P A I R' in line:
            contact_table = True
        if contact_table:
            match = re.match(r'^\s*-(\d+)-\s+-(\d+)-\s+-(\d+)-\s+([\d.Ee+\-]+)\s+FINITE SLIDING', line)
            if match:
                secondary, main, interaction, length = match.groups()
                contacts.append({'secondary_print_label': int(secondary), 'main_print_label': int(main),
                                 'interaction_print_label': int(interaction), 'characteristic_length_mm_printed': float(length),
                                 'dat_line': line_number})
    if not nodes or not ties:
        raise ValueError('required complete node map or Tie table missing')
    return nodes, ties, adjustments, contacts


def audit(dat, contract, node_csv, out):
    out = Path(out)
    if out.exists():
        raise ValueError('output exists; preserve previous audit')
    model = json.loads(Path(contract).read_text())
    nodes, table, adjustments, contacts = parse_dat(Path(dat).read_text(encoding='gbk'))
    converted = {}
    with open(node_csv, newline='') as stream:
        for row in csv.DictReader(stream):
            key = (row['instance'].upper(), int(row['abaqus_node_label']))
            if key in converted:
                raise ValueError('duplicate source instance/local node')
            converted[key] = int(row['exodus_node'])
    if set(nodes.values()) != set(converted):
        raise ValueError('DAT instance/local node identity differs from source mesh')
    global_to_exodus = {g: converted[key] for g, key in nodes.items()}
    exodus_to_global = {v: k for k, v in global_to_exodus.items()}
    if len(exodus_to_global) != len(nodes):
        raise ValueError('node map is not bijective')
    summaries = []
    assigned = set()
    rows = []
    for definition in model['interfaces']['ties']:
        name = definition['options']['name']
        secondary_name, main_name = definition['surfaces']
        secondary = {n for face in model['surfaces'][secondary_name] for n in face['nodes']}
        main = {n for face in model['surfaces'][main_name] for n in face['nodes']}
        bound = []
        unbound = []
        for node in sorted(secondary):
            label = exodus_to_global[node]
            entry = table.get(label)
            if entry:
                if label in assigned:
                    raise ValueError('one DAT secondary assigned to multiple source Tie surfaces')
                assigned.add(label)
                main_exodus = [global_to_exodus[g] for g in entry['main_global_labels']]
                if not set(main_exodus) <= main:
                    raise ValueError('printed Tie main node lies outside source main surface')
                bound.append(node)
            else:
                main_exodus = []
                unbound.append(node)
            rows.append({'tie': name, 'exodus_node_1based': node, 'dat_global_node': label,
                         'instance': nodes[label][0], 'abaqus_local_node': nodes[label][1],
                         'bound_in_printed_table': bool(entry), 'constraint': entry,
                         'main_exodus_nodes_1based': main_exodus})
        adjusted = adjustments.get('ASSEMBLY_' + name.upper(), [])
        if not {r['global_label'] for r in adjusted} <= {exodus_to_global[n] for n in bound}:
            raise ValueError('adjustment record not in actual bound membership')
        sizes = [len(table[exodus_to_global[n]]['main_global_labels']) for n in bound]
        summaries.append({'tie': name, 'secondary_surface': secondary_name, 'main_surface': main_name,
                          'candidate_nodes': len(secondary), 'bound_nodes': len(bound), 'unbound_nodes': len(unbound),
                          'bound_exodus_nodes_1based': bound, 'unbound_exodus_nodes_1based': unbound,
                          'main_nodes_per_bound_min': min(sizes, default=0), 'main_nodes_per_bound_max': max(sizes, default=0),
                          'printed_adjustment_count': len(adjusted),
                          'max_printed_adjustment_mm': max((r['distance_mm_printed'] for r in adjusted), default=None),
                          'adjustments_printed': adjusted})
    if assigned != set(table):
        raise ValueError('unassigned printed Tie secondary labels')
    result = {'schema_version': 1, 'status': 'printed-membership-audited',
              'inputs': {str(Path(p).resolve()): sha(p) for p in (dat, contract, node_csv)},
              'node_identity_match': {'nodes': len(nodes), 'instances': len({x[0] for x in nodes.values()}),
                                      'method': 'exact instance/local node identity set; not an INP-to-ODB cryptographic binding'},
              'printed_tie_rows': len(table), 'ties': summaries, 'contact_characteristic_lengths': contacts,
              'limitations': ['No Tie interpolation coefficients printed; embedded-element weight table is a different table.',
                              'Adjusted XYZ values are rounded DAT text, not a precise adjustment vector.',
                              'ODB was received and hashed but cannot be opened without Abaqus odbAccess.',
                              'Membership audit does not validate a MOOSE surface-to-surface Tie operator.']}
    out.mkdir(parents=True)
    (out/'tie-audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    (out/'tie-membership.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2)+'\n')
    with (out/'global-node-map.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['dat_global_node', 'instance', 'abaqus_local_node', 'exodus_node_1based'])
        writer.writerows((g, *nodes[g], global_to_exodus[g]) for g in sorted(nodes))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('dat', 'contract', 'node-csv', 'out'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args()
    result = audit(args.dat, args.contract, args.node_csv, args.out)
    print(json.dumps({x['tie']: [x['candidate_nodes'], x['bound_nodes'], x['unbound_nodes']] for x in result['ties']}))
