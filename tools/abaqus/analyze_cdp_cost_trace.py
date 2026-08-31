#!/usr/bin/env python3
"""Summarize bounded CDP diagnostic traces and complete per-rank counters.

Requires numpy. This does not label a nonsmooth tangent failure as a solver bug,
and does not sum nested inclusive timers into a wall-clock estimate.
"""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def newton_boundary_audit(matrix, residual, unknown):
    """Reconstruct a direction; do not alter state or claim nonlinear convergence."""
    matrix = np.asarray(matrix, dtype=float)
    residual = np.asarray(residual, dtype=float)
    unknown = np.asarray(unknown, dtype=float)
    try:
        direction = np.linalg.solve(matrix, -residual)
    except np.linalg.LinAlgError:
        return {"singular": True}
    decreasing = [i for i in range(6, 9) if direction[i] < 0.0]
    bounds = [-unknown[i] / direction[i] for i in decreasing]
    return {
        "singular": False,
        "scaled_unknowns_lambda_delta_kappa_t_delta_kappa_c": unknown[6:].tolist(),
        "scaled_newton_direction": direction[6:].tolist(),
        "maximum_nonnegative_step_capped_at_one": float(min([1.0] + bounds)),
        "blocked_zero_indices": [i for i in decreasing if unknown[i] == 0.0],
        "relative_linear_solve_residual": float(
            np.linalg.norm(matrix @ direction + residual)
            / max(np.linalg.norm(residual), 1e-30)
        ),
    }


def analyze(root):
    manifest=json.loads((root/'manifest.json').read_text())
    assert manifest['state']=='succeeded'
    ranks=defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
    for path in (root/'results').glob('cdp_cost_*.csv'):
        rank=int(re.search(r'_rank(\d+)_',path.name)[1])
        for row in csv.DictReader(path.open()):
            for key in ['calls','failed_calls','inclusive_microseconds','failed_inclusive_microseconds']:
                ranks[rank][row['category']][key]+=float(row[key])
    totals=defaultdict(lambda:defaultdict(float))
    for rank in ranks.values():
        for category,values in rank.items():
            for key,value in values.items():totals[category][key]+=value
    events=[]
    for path in (root/'logs').glob('cdp_trace_*.jsonl'):
        rank=int(re.search(r'_rank(\d+)_',path.name)[1])
        for line in path.read_text().splitlines():
            event=json.loads(line);event['rank']=rank;events.append(event)
    failures=[];tangents=[]
    for event in events:
        identity={k:event[k] for k in ['rank','time','dt','step','element','qp','call','partition','substep']}
        if event['event']=='local_failure_jacobian':
            if event.get('diagnostic_error'):
                failures.append(dict(identity,diagnostic_error=True));continue
            ad=np.array(event['ad'],float); sweeps=[]
            for sample in event['difference_sweep']:
                norms={}
                for key in ['centered','forward']:
                    fd=np.array(sample[key],float)
                    norms[key+'_relative_error']=float(np.linalg.norm(ad-fd)/max(np.linalg.norm(fd),1e-30))
                sweeps.append(dict(h=sample['h'],crosses_nonnegative_unknown_boundary=bool(min(event['unknown'][6:])<sample['h']),
                    centered_newton_boundary=newton_boundary_audit(sample['centered'], event['residual'], event['unknown']),**norms))
            s=np.array(event['unknown'][:6])*event['stress_scale']
            eig=np.linalg.eigvalsh([[s[0],s[3],s[5]],[s[3],s[1],s[4]],[s[5],s[4],s[2]]])
            failures.append(dict(identity,residual_inf=float(max(abs(np.array(event['residual'])))),
                principal_stress_pa=eig.tolist(),minimum_relative_principal_gap=float(min(np.diff(eig))/max(np.linalg.norm(eig),1.)),
                minimum_scaled_history_unknown=min(event['unknown'][6:]),
                ad_newton_boundary=newton_boundary_audit(ad, event['residual'], event['unknown']),
                difference_sweep=sweeps))
        if event['event']=='material_tangent':
            groups=[]
            for h in sorted({s['h'] for s in event['samples']},reverse=True):
                samples=[s for s in event['samples'] if s['h']==h]
                valid=[s for s in samples if s['ok'] and s['same_final_branch'] and
                       s['plus_partition']==s['minus_partition']==event['base_partition']]
                groups.append(dict(h=h,successful_directions=sum(s['ok'] for s in samples),
                    same_partition_and_final_branch_directions=len(valid),
                    max_relative_error=max((s['relative_error'] for s in valid),default=None),
                    errors_by_column={s['column']:s.get('relative_error') for s in samples}))
            tangents.append(dict(identity,base_partition=event['base_partition'],
                repeat_max_stress_error=event.get('repeat_max_stress_error'),
                repeat_same_partition=event.get('repeat_same_partition'),perturbation_sweep=groups))
    denominator=totals['material']['inclusive_microseconds']
    ratio=lambda category,key='inclusive_microseconds':totals[category][key]/denominator if denominator else None
    status=json.loads((root/'execution-status.json').read_text())
    result={'job_id':manifest['job_id'],'solver_sha':manifest['solver_versions']['dam_safety_app'],
        'api_timings':status['timings'],'counts_and_inclusive_times_by_rank':dict(ranks),
        'aggregate_counters':dict(totals),'event_counts':dict(Counter(e['event'] for e in events)),
        'fractions_of_material_inclusive_time':{
            'failed_partitions':ratio('partition','failed_inclusive_microseconds'),
            'local_integrate':ratio('local'),'local_linearized':ratio('local_linearized'),
            'ad_jacobian':ratio('ad_jacobian'),'fd_jacobian':ratio('fd_jacobian'),
            'ordinary_spectrum':ratio('spectrum'),'ad_spectrum':ratio('ad_spectrum'),
            'factor':ratio('factor'),'backsolve':ratio('backsolve'),'state_chain':ratio('state_chain')},
        'failed_local_jacobian_samples':failures,'material_tangent_samples':tangents,
        'limitations':['Inclusive categories overlap; do not add them.',
            'Trace samples are real trial material evaluations, not necessarily accepted global states.',
            'Matching final branches does not prove all internal branches are smooth.',
            'Offline Newton directions identify linearized feasibility only; they do not prove a projected or active-set nonlinear update converges.',
            'Derivative perturbations are diagnostic recomputations and are excluded from production counters/times.',
            'One diagnostic timing is not a production benchmark; MPI waiting and other framework costs remain outside material timers.']}
    (root/'cost-tangent-summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args()
    result=analyze(a.root)
    print(json.dumps({k:result[k] for k in ['job_id','event_counts','fractions_of_material_inclusive_time','material_tangent_samples']},ensure_ascii=False))
