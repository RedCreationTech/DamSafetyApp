#!/usr/bin/env python3
"""Prepare one immutable boundary-line-search comparison; never submit a Job."""
import argparse
import hashlib
import json
from pathlib import Path

CASES = {
    'B00-single-off': ('D04c-early-window', False),
    'B01-single-on': ('D04c-early-window', True),
    'B02-single-refined-on': ('D04g-early-refined-grid', True),
    'B03-fullmesh-off': ('D05b-fullmesh-short-observe', False),
    'B04-fullmesh-on': ('D05b-fullmesh-short-observe', True),
    'B05-fullmesh-refined-on': ('D05c-fullmesh-early-refined', True),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(parent_root, root, label, audit_path):
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
    a = p.parse_args()
    print(json.dumps(prepare(a.parents, a.root, a.case, a.release_audit)))
