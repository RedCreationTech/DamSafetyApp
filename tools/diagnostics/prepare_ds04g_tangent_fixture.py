#!/usr/bin/env python3
"""Freeze the seven already-selected observation inputs as numeric unit fixtures."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def values(row):
    return row['values'].split()


def unique(rows, **matches):
    matches = [r for r in rows if all(str(r[k]) == str(v) for k, v in matches.items())]
    assert len(matches) == 1, len(matches)
    return matches[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selected', type=Path, required=True)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    states = json.loads(args.selected.read_text())
    assert [s['id'] for s in states] == ['F1', 'F2', 'F3', 'F4', 'L1', 'L2', 'L3']
    manifest = dict(parent_job_id='job_20261020_223540_ds33p4',
                    source_solver_sha='8a0114cff87eac51b5db1a796434906301388d0e',
                    selected_sha256=hashlib.sha256(args.selected.read_bytes()).hexdigest(), files={}, points=[])
    for s in states:
        material = [r for r in s['raw']['qp-updates'] if r['event'] == 'material_enter']
        get = lambda field: values(unique(material, field=field))
        scalars = get('old_kappa_DT_DC_dt_full_tangent')
        old = get('old_backbone_plastic') + scalars[:2] + get('old_viscous_plastic') + scalars[2:4]
        tangent = unique(s['raw']['algorithmic-tangent'], field='column_major_6x6')
        final = [r for r in s['raw']['substep-events'] if r['event'] == 'state_return' and
                 r['partition'] == tangent['partition'] and r['substep'] == tangent['N']]
        out = lambda field: values(unique(final, field=field))
        ending = out('kappa_DT_DC_backboneDT_backboneDC_dt_mu_stiffness')
        new = out('backbone_plastic') + ending[:2] + out('viscous_plastic') + ending[2:4]
        lines = [s['id'], ' '.join(get('old_total_strain')), ' '.join(get('new_total_strain')),
                 ' '.join(old), scalars[4] + ' ' + tangent['N'],
                 ' '.join(out('cauchy_stress')), ' '.join(new), ' '.join(values(tangent)),
                 '1' if s['local_solve'] is not None else '0']
        if s['local_solve'] is not None:
            selected = [r for r in s['raw']['local-newton'] if int(r['local_solve']) == s['local_solve']]
            enter = [r for r in selected if r['event'] == 'enter']
            trial = unique(selected, event='trial', field='scales_yield_elastic')
            iterate = [r for r in selected if r['event'] == 'iterate' and int(r['local_iteration']) == s['local_iteration']]
            lines += [' '.join(values(unique(enter, field='total_strain'))),
                      ' '.join(values(unique(enter, field='old_plastic')) + values(unique(enter, field='old_kappa'))),
                      ' '.join(values(trial)[:2])]
            for field in ['unknown', 'residual']:
                lines.append(' '.join(values(unique(iterate, field=field))))
            for event in ['AD_Newton', 'FD_fallback_Newton']:
                rows = [r for r in s['raw']['local-jacobian'] if int(r['local_solve']) == s['local_solve'] and
                        int(r['local_iteration']) == s['local_iteration']]
                lines.append(' '.join(values(unique(rows, event=event, field='row_major_9x9'))))
        assert len(old) == len(new) == 16
        (args.output / (s['id'] + '.txt')).write_text('\n'.join(lines) + '\n', encoding='ascii')
        manifest['points'].append(dict(id=s['id'], identity=s['key'], summary=s['summary']))
    for name in ['compression_hardening.csv', 'compression_damage.csv', 'tension_stiffening.csv', 'tension_damage.csv']:
        shutil.copyfile(args.input / name, args.output / name)
    for path in sorted(args.output.iterdir()):
        if path.suffix in ('.txt', '.csv'):
            manifest['files'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.output / 'source-identity.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Frozen 7 numeric fixtures; no integrator evaluation')


if __name__ == '__main__':
    main()
