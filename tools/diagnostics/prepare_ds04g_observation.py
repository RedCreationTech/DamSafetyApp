#!/usr/bin/env python3
"""Freeze original OFF prefix and passive observation input; never submit or run a solver."""
import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path

SOURCE_SHA = '9813ae91f1433b7798fec4233be0efdb23739009e544405247f04cac9a9da804'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(source, nodes_csv, elements_csv, output):
    assert sha(source) == SOURCE_SHA, 'Original OFF20 input changed'
    nodes = {int(r['node_label']): [float(r['X_mm'])/1000, float(r['Y_mm'])/1000] for r in csv.DictReader(nodes_csv.open(encoding='utf-8-sig'))}
    elements = {int(r['element_label']): json.loads(r['connectivity']) for r in csv.DictReader(elements_csv.open(encoding='utf-8-sig'))}
    labels = [1881, 767]
    xy = [v for label in labels for n in elements[label] for v in nodes[n]]
    original = source.read_text()
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

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['source','nodes_csv','elements_csv','output']: parser.add_argument(name, type=Path)
    args=parser.parse_args()
    print(json.dumps(prepare(args.source,args.nodes_csv,args.elements_csv,args.output),ensure_ascii=False,indent=2))
