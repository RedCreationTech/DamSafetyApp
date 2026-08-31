#!/usr/bin/env python3
"""Audit accepted raw-IP probe CSVs (requires numpy); no solver or API calls."""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np


def tensor(row, prefix):
    a = np.zeros((3, 3))
    for i, j, label in [(0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')]:
        a[i,j] = a[j,i] = row['probe_'+prefix+'_'+label]
    return a


def elastic(strain):
    e, nu = 29791500000., .2
    return e/(1+nu)*strain+e*nu/((1+nu)*(1-2*nu))*np.trace(strain)*np.eye(3)


def invariants(stress):
    values = np.linalg.eigvalsh(stress)
    r = max(0.,sum(np.maximum(values,0)))/sum(abs(values)) if any(values) else 0.
    return values, r


def analyze(root):
    manifest = json.loads((root/'manifest.json').read_text())
    assert manifest['state'] == 'succeeded'
    log = re.sub(r'\x1b\[[0-9;]*m','',(root/'logs/solve.log').read_text())
    history_path=root/'results/compression_isolation.csv'
    single_element=history_path.exists()
    if not single_element: history_path=root/'results/uniaxial_compression.csv'
    history = list(csv.DictReader(history_path.open()))
    v = lambda key: [float(r[key]) for r in history]
    times = {int(r['timestep']):float(r['time']) for r in csv.DictReader((root/'results/ip_history_ip_states_time.csv').open())}
    rows, previous = [], {}
    for path in sorted((root/'results').glob('ip_history_ip_states_*.csv')):
        if path.stem.endswith('_time'): continue
        step = int(path.stem.rsplit('_',1)[1]); time = times[step]
        for raw in csv.DictReader(path.open()):
            row = {k:float(v) for k,v in raw.items()}
            key = (int(row['elem_id']),int(row['qp_id']))
            sigma, eps, pb, pv = [tensor(row,p) for p in ['sigma','eps','pb','pv']]
            bone, visc = elastic(eps-pb), elastic(eps-pv)
            eig_b,r_b = invariants(bone); eig_v,r_v = invariants(visc)
            predicted = row['cdp_stiffness_factor']*visc
            prior = previous.get(key)
            result = dict(time=time, step=step, elem_id=key[0], qp_id=key[1],
                damage_t=row['DamageT'],damage_c=row['DamageC'],kappa_t=row['cdp_kappa_t'],
                kappa_c=row['cdp_kappa_c'],r_backbone=r_b,r_viscous=r_v,
                backbone_max_principal=float(eig_b[-1]),viscous_max_principal=float(eig_v[-1]),
                lateral_abs=max(abs(sigma[0,0]),abs(sigma[1,1])),
                shear_abs=max(abs(sigma[0,1]),abs(sigma[0,2]),abs(sigma[1,2])),
                axial_abs=abs(sigma[2,2]),
                reconstructed_cauchy_max_error=float(np.max(abs(predicted-sigma))),
                delta_kappa_t=row['cdp_kappa_t']-(prior['kappa_t'] if prior else 0.),
                accepted_substeps=row['cdp_accepted_substeps'],
                ad_evaluations=row['cdp_automatic_jacobian_evaluations'],
                fd_evaluations=row['cdp_finite_difference_jacobian_evaluations'])
            previous[key] = result; rows.append(result)
    element_ids=sorted({r['elem_id'] for r in rows})
    assert len(rows) == len(times)*8*len(element_ids), 'Requires 8 IPs on every sampled HEX8 at every accepted step'
    peak = max(r['axial_abs'] for r in rows)
    max_value = lambda k: max(r[k] for r in rows)
    first = next((r for r in rows if r['kappa_t']>1e-12),None)
    reasons = Counter(re.findall(r'Nonlinear solve converged due to (\w+)',log))
    residuals = []
    for block in re.split(r'\nTime Step ',log)[1:]:
        if 'Solve Converged!' in block:
            values = re.findall(r'\d+ Nonlinear \|R\| = ([0-9.eE+-]+)',block)
            if values: residuals.append(float(values[-1]))
    result = dict(job_id=manifest['job_id'],state=manifest['state'],
        finished_executing_seconds=(float(m.group(1)) if (m := re.search(r'Finished Executing\s*\[\s*([0-9.]+) s',log)) else None),
        api_solver_elapsed_seconds=json.loads((root/'execution-status.json').read_text())['timings']['solver_elapsed_seconds'],
        accepted_steps=manifest['progress']['step_current'],
        attempted_steps=len(re.findall(r'^Time Step [1-9]\d*, time =',log,re.M)),
        cutbacks=len(re.findall(r'^Time Step [1-9]\d*, time =',log,re.M))-manifest['progress']['step_current'],
        convergence_reasons=dict(reasons),max_accepted_logged_residual=max(residuals),
        ip_frames=len(times),ip_rows=len(rows),sampled_element_ids=element_ids,max_ip_damage_t=max_value('damage_t'),
        max_ip_kappa_t=max_value('kappa_t'),max_ip_lateral_pa=max_value('lateral_abs'),
        max_ip_shear_pa=max_value('shear_abs'),ip_lateral_to_axial_peak=max_value('lateral_abs')/peak,
        max_ip_r_backbone=max_value('r_backbone'),max_ip_r_viscous=max_value('r_viscous'),
        max_cauchy_reconstruction_error_pa=max_value('reconstructed_cauchy_max_error'),
        first_kappa_t_above_1e_12=first,
        mean_curve_max_damage_t=max(v('average_damage_t')) if single_element else None,
        mean_curve_peak_compression_pa=max(abs(x) for x in v('average_stress_zz')) if single_element else None,
        final_mean_damage_c=v('average_damage_c')[-1] if single_element else None,
        final_mean_compression_pa=-v('average_stress_zz')[-1] if single_element else None,
        max_ip_accepted_substeps=max_value('accepted_substeps'),
        max_ip_fd_evaluations=max_value('fd_evaluations'),
        failure_branches=dict(Counter(re.findall(r'branch=(\w+)',log))),
        caveat='Accepted end states only. Reconstructed backbone stress uses raw strain/plastic tensors and locked E/nu; r is not the unobserved substep history. Timing is one diagnostic run, not a benchmark.')
    result['strict_pure_compression_screen_passed'] = bool(
        result['ip_lateral_to_axial_peak'] <= 1e-6 and result['max_ip_shear_pa']/peak <= 1e-6
        and result['max_ip_damage_t'] <= 1e-8 and result['max_ip_kappa_t'] <= 1e-8)
    variant=json.loads((root/'input/diagnostic-manifest.json').read_text())['test_variant']
    result['pure_compression_screen_applicable']=single_element and variant not in ('top-clamp','replay-path')
    if not result['pure_compression_screen_applicable']:
        result['strict_pure_compression_screen_passed']=None
    checks = [{'name':manifest['input_snapshot']['input_file'],'sha256':manifest['input_snapshot']['input_sha256']}]+manifest['input_snapshot']['mesh_files']+manifest['input_snapshot']['extra_files']
    assert all(hashlib.sha256((root/'input'/r['name']).read_bytes()).hexdigest()==r['sha256'] for r in checks)
    result['input_hashes_verified'] = len(checks)
    with (root/'accepted-ip-derived.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (root/'probe-summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args()
    result=analyze(a.root)
    print(json.dumps({k:v for k,v in result.items() if k not in ['first_kappa_t_above_1e_12','caveat']},ensure_ascii=False))
