#!/usr/bin/env python3
"""Freeze original OFF prefix and passive observation input; never submit or run a solver."""
import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

SOURCE_SHA = '9813ae91f1433b7798fec4233be0efdb23739009e544405247f04cac9a9da804'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(source, nodes_csv, elements_csv, output, *, labels=(1881, 767),
            source_sha=SOURCE_SHA, preserve_source=False, observation_start=2.698,
            observation_end=2.720, byte_limit=65536000):
    assert sha(source) == source_sha, 'Frozen source input changed'
    labels = list(labels)
    assert labels and len(set(labels)) == len(labels) and all(x > 0 for x in labels)
    assert math.isfinite(observation_start) and math.isfinite(observation_end)
    assert observation_start <= observation_end and byte_limit >= 0
    with nodes_csv.open(encoding='utf-8-sig') as f:
        nodes = {int(r['node_label']): [float(r['X_mm'])/1000, float(r['Y_mm'])/1000] for r in csv.DictReader(f)}
    with elements_csv.open(encoding='utf-8-sig') as f:
        elements = {int(r['element_label']): json.loads(r['connectivity']) for r in csv.DictReader(f)}
    assert all(len(elements[k]) == 4 and len(set(elements[k])) == 4 for k in labels)
    xy = [v for label in labels for n in elements[label] for v in nodes[n]]
    assert all(math.isfinite(v) for v in xy)
    original = source.read_text()
    if preserve_source:
        return prepare_existing(source, original, nodes_csv, elements_csv, output,
                                labels, xy, elements, observation_start, observation_end, byte_limit)
    cut = original.replace('  end_time = 21.0\n', '  end_time = 2.720000\n')
    assert cut != original
    edits = []
    def trim(match):
        prefix, seq = match.groups()
        values = seq.split()
        keep = [x for x in values if float(x) <= 2.7200000001]
        assert keep and abs(float(keep[-1])-2.720) < 1e-12
        edits.append({'parameter': prefix.strip(), 'original_count': len(values), 'prefix_count': len(keep)})
        return prefix + "'" + ' '.join(keep) + "'"
    cut = re.sub(r"(^\s*(?:time_sequence|sync_times)\s*=\s*)'([^']*)'", trim, cut, flags=re.M)
    assert len(edits) == 3, edits
    assert 'enable_finite_newmark_damping = false' in cut
    control = cut.replace('    type = TimeSequenceStepper\n', '    type = CDPObservationTimeStepper\n')
    assert control != cut
    extra = '\n[Problem]\n  type = CDPObservationProblem\n  observation_enabled = true\n  observation_start = 2.698\n  observation_end = 2.720\n  observation_labels = \'1881 767\'\n  observation_vertices = \'' + ' '.join(format(v, '.17g') for v in xy) + "'\n[]\n\n[Materials/passive_observer]\n  type = CDPObservationMaterial\n[]\n"
    extra = extra.replace("'1881 767'", "'" + ' '.join(map(str, labels)) + "'")
    extra = extra.replace('observation_start = 2.698', f'observation_start = {observation_start:.3f}')
    extra = extra.replace('observation_end = 2.720', f'observation_end = {observation_end:.3f}')
    if byte_limit != 65536000:
        extra = extra.replace('  observation_enabled = true\n',
                              f'  observation_enabled = true\n  observation_byte_limit = {byte_limit}\n')
    observed = control + extra
    assert observed.removesuffix(extra).replace('    type = CDPObservationTimeStepper\n','    type = TimeSequenceStepper\n') == cut
    output.mkdir(parents=True, exist_ok=False)
    for p in source.parent.iterdir():
        if p.suffix in {'.csv', '.e'}:
            shutil.copyfile(p, output/p.name)
    (output/'ds04g_observe_eq1p720.i').write_text(observed)
    # Same candidate / same diagnostic objects; this arm exercises the disabled default.
    (output/'ds04g_observe_disabled.i').write_text(observed.replace('  observation_enabled = true\n','  observation_enabled = false\n'))
    (output/'ds04g_prefix_control.i').write_text(cut)
    receipt = {'source_input_sha256': SOURCE_SHA, 'no_job_submitted': True,
               'eq_range_s': [0, 1.720], 'cumulative_end_time_s': 2.720,
               'unchanged_physics_and_prefix_restoration_exact': True,
               'sequence_truncations': edits, 'source_connectivity': {k: elements[k] for k in labels},
               'source_geometry_sha256': {str(nodes_csv): sha(nodes_csv), str(elements_csv): sha(elements_csv)},
               'files_sha256': {p.name: sha(p) for p in sorted(output.iterdir())}}
    (output/'input-audit.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n')
    return receipt


def prepare_existing(source, original, nodes_csv, elements_csv, output, labels, xy,
                     elements, start, end, byte_limit):
    """Only replace the already-present observation block; never regenerate solver times."""
    pattern = r'^\[Problem\]\n.*?^\[\]\n'
    blocks = re.findall(pattern, original, flags=re.M | re.S)
    assert len(blocks) == 1 and 'type = CDPObservationProblem' in blocks[0]
    assert 'enable_finite_newmark_damping = false' in original
    assert 'type = CDPObservationTimeStepper' in original
    block = blocks[0]
    values = {'observation_enabled': 'true', 'observation_start': format(start, '.17g'),
              'observation_end': format(end, '.17g'),
              'observation_labels': "'" + ' '.join(map(str, labels)) + "'",
              'observation_vertices': "'" + ' '.join(format(v, '.17g') for v in xy) + "'"}
    for key, value in values.items():
        block, count = re.subn(r'^(  '+key+r'\s*=\s*)[^\n]*',
                               lambda m: m[1] + value, block, flags=re.M)
        assert count == 1, 'Missing or duplicate observation parameter: ' + key
    if re.search(r'^  observation_byte_limit\s*=', block, flags=re.M):
        block = re.sub(r'^(  observation_byte_limit\s*=\s*)[^\n]*',
                       lambda m: m[1] + str(byte_limit), block, flags=re.M)
    else:
        block = block[:-3] + '  observation_byte_limit = ' + str(byte_limit) + '\n[]\n'
    observed = original.replace(blocks[0], block, 1)
    disabled = observed.replace('  observation_enabled = true\n',
                                '  observation_enabled = false\n', 1)
    assert observed.replace(block, blocks[0], 1) == original
    output.mkdir(parents=True, exist_ok=False)
    for p in source.parent.iterdir():
        if p.suffix in {'.csv', '.e'}:
            shutil.copyfile(p, output/p.name)
    for name, text in [('ds04g_prefix_control.i', original),
                       ('ds04g_observe_disabled.i', disabled),
                       ('ds04g_observe_eq1p720.i', observed)]:
        (output/name).write_text(text)
    receipt = {'source_input_sha256': sha(source), 'no_job_submitted': True,
               'mode': 'preserve_source_observation_block_only', 'labels': labels,
               'expected_native_QPs': 4*len(labels), 'observation_window_total_s': [start, end],
               'observation_byte_limit_per_rank': byte_limit,
               'outside_observation_block_byte_identical': True,
               'source_connectivity': {k: elements[k] for k in labels},
               'source_geometry_sha256': {str(p): sha(p) for p in (nodes_csv, elements_csv)},
               'files_sha256': {p.name: sha(p) for p in sorted(output.iterdir())}}
    (output/'input-audit.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n')
    return receipt

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['source','nodes_csv','elements_csv','output']: parser.add_argument(name, type=Path)
    parser.add_argument('--labels', nargs='+', type=int, default=[1881, 767])
    parser.add_argument('--source-sha', default=SOURCE_SHA)
    parser.add_argument('--preserve-source', action='store_true')
    parser.add_argument('--observation-start', type=float, default=2.698)
    parser.add_argument('--observation-end', type=float, default=2.720)
    parser.add_argument('--byte-limit', type=int, default=65536000)
    args=parser.parse_args()
    print(json.dumps(prepare(args.source,args.nodes_csv,args.elements_csv,args.output,
        labels=args.labels, source_sha=args.source_sha, preserve_source=args.preserve_source,
        observation_start=args.observation_start, observation_end=args.observation_end,
        byte_limit=args.byte_limit),ensure_ascii=False,indent=2))
