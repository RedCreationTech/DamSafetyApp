#!/usr/bin/env python3
"""Prepare one immutable boundary-line-search comparison; never submit a Job."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

CASES = {
    'B00-single-off': ('D04c-early-window', False),
    'B01-single-on': ('D04c-early-window', True),
    'B02-single-refined-on': ('D04g-early-refined-grid', True),
    'B03-fullmesh-off': ('D05b-fullmesh-short-observe', False),
    'B04-fullmesh-on': ('D05b-fullmesh-short-observe', True),
    'B05-fullmesh-refined-on': ('D05c-fullmesh-early-refined', True),
    'B06-fullmesh-baseline-grid-on': ('D05b-fullmesh-short-observe', True),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(parent_root, root, label, audit_path, grid_job_root=None):
    parent_name, enabled = CASES[label]
    parent = parent_root / parent_name
    old = json.loads((parent / 'submission.json').read_text())
    audit = json.loads(audit_path.read_text())
    source = parent / 'input' / old['input_file']
    assert sha(source) == old['input_sha256']
    text = source.read_text()
    marker = '    type = AbaqusCDPStressUpdate\n'
    assert text.count(marker) == 1
    assert 'use_bound_feasible_line_search' not in text
    # Timing comparison keeps previous material diagnostics/output unchanged.
    assert 'enable_path_diagnostics = true' not in text
    text = text.replace(marker, marker + '    use_bound_feasible_line_search = '
                        + str(enabled).lower() + '\n')
    grid_source = None
    grid_manifest = None
    if label == 'B06-fullmesh-baseline-grid-on':
        assert grid_job_root is not None, 'A completed control Job is required'
        grid_manifest = json.loads((grid_job_root / 'manifest.json').read_text())
        assert grid_manifest['state'] == 'succeeded'
        grid_source = grid_job_root / 'results/ip_history_ip_states_time.csv'
        rows = list(csv.DictReader(grid_source.open()))
        times = [float(row['time']) for row in rows]
        assert times[0] == 0 and abs(times[-1] - 0.1) < 1e-14
        assert all(b > a for a, b in zip(times, times[1:]))
        sequence = ' '.join(row['time'] for row in rows[1:])
        text, count = re.subn(r'  \[TimeStepper\]\n.*?\n  \[\]',
                             "  [TimeStepper]\n    type = TimeSequenceStepper\n"
                             "    time_sequence = '" + sequence + "'\n  []", text, flags=re.S)
        assert count == 1
    out = root / label
    out.mkdir(parents=True, exist_ok=False)
    inp = out / 'input'
    inp.mkdir()
    case_name = 'tpl-cdpc-uc-' + label.lower() + '-v1'
    filename = case_name + '.i'
    (inp / filename).write_text(text)
    for item in old['mesh_files'] + old['extra_files']:
        if item['name'] in {'diagnostic-manifest.json', 'release-audit.json'}:
            continue
        path = parent / 'input' / item['name']
        assert sha(path) == item['sha256']
        (inp / item['name']).write_bytes(path.read_bytes())
    (inp / 'release-audit.json').write_bytes(audit_path.read_bytes())
    if grid_source:
        (inp / 'baseline-accepted-time-grid.csv').write_bytes(grid_source.read_bytes())
    diagnostic = {
        'test_id': 'TEST-CDPC-BOUND-001/' + label,
        'test_variant': json.loads((parent / 'input' / 'diagnostic-manifest.json').read_text())['test_variant'],
        'parent_job_id': json.loads((parent / 'manifest.json').read_text())['job_id'],
        'parent_input_sha256': old['input_sha256'],
        'solver_sha': audit['solver_sha'],
        'use_bound_feasible_line_search': enabled,
        'change': 'Only solver release and explicit boundary-line-search switch; tables, loading, numerical tolerances, time grid and outputs unchanged from parent',
        'status': 'prototype-candidate; not engineering acceptance',
    }
    if grid_source:
        diagnostic.update(grid_reference_job=grid_manifest['job_id'],
                          grid_source_sha256=sha(grid_source),
                          change='Boundary candidate with the completed control Job accepted time grid; retains failed-step cutback capability, original material, tolerances, loading and output times')
    (inp / 'diagnostic-manifest.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
    mesh = {item['name'] for item in old['mesh_files']}
    submission = dict(old, case_name=case_name, input_file=filename,
                      input_sha256=sha(inp / filename), command='DamSafetyApp-opt -i ' + filename)
    submission['extra_files'] = [dict(name=p.name, sha256=sha(p)) for p in sorted(inp.iterdir())
                                 if p.name != filename and p.name not in mesh]
    (out / 'submission.json').write_text(json.dumps(submission, indent=2) + '\n')
    return {'case': label, 'root': str(out), 'input_sha256': submission['input_sha256']}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parents', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--case', choices=CASES, required=True)
    p.add_argument('--release-audit', type=Path, required=True)
    p.add_argument('--grid-job-root', type=Path)
    a = p.parse_args()
    print(json.dumps(prepare(a.parents, a.root, a.case, a.release_audit, a.grid_job_root)))
