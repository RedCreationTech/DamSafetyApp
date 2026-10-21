#!/usr/bin/env python3
"""Audit passive traces; never run a material integrator or submit a simulation."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def rows(path):
    with path.open(newline='') as f:
        yield from csv.DictReader(f)


def audit(directory, prefix, job_id=None, development=False, *, labels=(1881, 767),
          mpi_ranks=4, require_reject=True, window=None, nominal_dt=None):
    labels = list(labels)
    if not labels or len(set(labels)) != len(labels) or any(x <= 0 for x in labels):
        raise ValueError('Expected unique positive source labels')
    if mpi_ranks <= 0 or (nominal_dt is not None and (not math.isfinite(nominal_dt) or nominal_dt <= 0)):
        raise ValueError('Invalid rank count or nominal dt')
    if window is not None and (len(window) != 2 or not all(map(math.isfinite, window)) or window[0] > window[1]):
        raise ValueError('Invalid observation window')
    statuses=[json.loads(p.read_text()) for p in sorted(directory.glob(prefix+'-status-rank*.json'))]
    errors=[]
    if sorted(s['rank'] for s in statuses)!=list(range(mpi_ranks)):errors.append('missing_or_duplicate_MPI_status')
    for s in statuses:
        for k in ['capture_truncated','io_failed','history_incomplete']:
            if s[k]:errors.append('rank%d:%s'%(s['rank'],k))
    points=set();owners={};sizes={};events={};weights={};matrix_calls=set();raw_calls=set()
    references={};accepted_times=set();seen={s['rank']:bytearray(s['records']+1) for s in statuses};point_times={};local_matrices=0
    files=sorted(directory.glob(prefix+'-*.csv'))
    for p in files:
        rank=int(p.stem.rsplit('rank',1)[1]); count=0
        for r in rows(p):
            count+=1; event=r['event'];events[event]=events.get(event,0)+1
            label=int(r['source_element']);qp=int(r['qp'])
            record=int(r['record'])
            if rank not in seen or record <= 0 or record >= len(seen[rank]) or seen[rank][record]:errors.append('duplicate_or_invalid_record')
            else:seen[rank][record]=1
            if event=='accepted_end':accepted_times.add(float(r['time']))
            if event=='qp' and r['field']=='x_y_z_JxW_coord':
                points.add((label,qp)); owners.setdefault(label,set()).add(rank)
                v=[float(x) for x in r['values'].split()]
                if len(v)!=5 or not all(map(math.isfinite,v)):errors.append('invalid_QP_geometry');continue
                weight=v[3]*v[4]
                if (label,qp) in weights and not math.isclose(weights[label,qp],weight,rel_tol=1e-12,abs_tol=0):errors.append('changing_QP_weight')
                weights[label,qp]=weight
                point_times.setdefault((label,qp),set()).add(float(r['time']))
            if event=='qp' and r['field']=='reference_x_y_z':
                v=tuple(float(x) for x in r['values'].split())
                if len(v)!=3 or not all(map(math.isfinite,v)):errors.append('invalid_reference_QP')
                references.setdefault((label,qp),set()).add(v)
            key=(rank,int(r['evaluation']),int(r['call']),label,qp)
            if r['field']=='column_major_6x6':
                raw_calls.add(key)
                if len(r['values'].split())!=36:errors.append('bad_6x6_shape')
            if r['field']=='RankFour_ijkl':
                matrix_calls.add(key)
                if int(r['call'])==0:errors.append('unassociated_assembly_tangent')
                if len(r['values'].split())!=81:errors.append('bad_RankFour_shape')
            if r['field']=='row_major_9x9':
                local_matrices+=1
                if len(r['values'].split())!=81:errors.append('bad_9x9_shape')
            if r['field'] in ('column_major_6x6','RankFour_ijkl','row_major_9x9') and not all(math.isfinite(float(v)) for v in r['values'].split()):errors.append('nonfinite_matrix')
        sizes[p.name]={'rows':count,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    expected={(e,q) for e in labels for q in range(4)}
    if points!=expected:errors.append('expected_declared_elements_times_four_native_QPs')
    if set(references)!=expected or any(len(v)!=1 for v in references.values()):errors.append('missing_or_changing_reference_QP')
    for e in labels:
        positions={next(iter(references[e,q])) for q in range(4) if (e,q) in references and references[e,q]}
        if len(positions)!=4:errors.append('nonunique_reference_QPs')
    if any(len(x)!=1 for x in owners.values()):errors.append('duplicate_owner')
    if any(w<=0 for w in weights.values()):errors.append('nonpositive_quadrature_weight')
    if matrix_calls-raw_calls:errors.append('assembly_tangent_has_no_matching_material_call')
    for status in statuses:
        if sum(seen[status['rank']]) != status['records']:errors.append('missing_record_sequence')
        actual={p.name[len(prefix)+1:].rsplit('-rank',1)[0]:sizes[p.name]['rows'] for p in files if p.stem.endswith('rank'+str(status['rank']))}
        if actual != status['streams']:errors.append('status_stream_count_mismatch')
    if not development:
        if not local_matrices:errors.append('missing_plastic_local_Jacobian')
        if not matrix_calls:errors.append('missing_actual_assembly_tangent')
        if require_reject and not events.get('reject_after'):errors.append('no_rejected_step_restore_observed')
        if not events.get('accepted_end'):errors.append('no_accepted_step_observed')
        if window is not None:
            required=[window[0],window[1]]
            if nominal_dt is not None:
                count=round((window[1]-window[0])/nominal_dt)
                if not math.isclose(window[0]+count*nominal_dt,window[1],abs_tol=1e-9,rel_tol=0):raise ValueError('Window is not on nominal target grid')
                required=[window[0]+i*nominal_dt for i in range(count+1)]
            if any(not any(abs(t-x)<1e-9 for x in accepted_times) for t in required):errors.append('missing_accepted_target_in_window')
            if any(not any(abs(t-x)<1e-9 for x in point_times.get(k,set())) for k in expected for t in required):errors.append('missing_native_QP_target_in_window')
    result={'schema':'ds04g-observation-audit-v1','job_id':job_id,'development_only':development,
            'complete_for_declared_scope':not errors,'errors':sorted(set(errors)),
            'statuses':statuses,'native_QPs':sorted(points),'owners':{k:sorted(v) for k,v in owners.items()},
            'normalized_weights':{str(e):[weights.get((e,q),0)/sum(weights.get((e,j),0) for j in range(4)) for q in range(4)] for e in owners if sum(weights.get((e,j),0) for j in range(4))>0},
            'events':events,'files':sizes,'paired_assembly_tangent_calls':len(matrix_calls & raw_calls),
            'declared_labels':labels,'expected_native_QPs':len(expected),
            'reference_QPs':{str(k):sorted(v) for k,v in references.items()},
            'required_window_total_s':window,'nominal_dt_s':nominal_dt,
            'rollback_observed':bool(events.get('reject_after')),
            'scope_note':'Structural trace integrity only; does not prove algorithmic tangent accuracy, global rollback, unchanged trajectory or engineering acceptance'}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--prefix',default='ds04g-observation')
    p.add_argument('--labels',nargs='+',type=int,default=[1881,767]);p.add_argument('--mpi-ranks',type=int,default=4)
    p.add_argument('--allow-no-reject',action='store_true');p.add_argument('--window',nargs=2,type=float);p.add_argument('--nominal-dt',type=float)
    group=p.add_mutually_exclusive_group(required=True);group.add_argument('--job-id');group.add_argument('--development',action='store_true')
    a=p.parse_args();r=audit(a.directory,a.prefix,a.job_id,a.development,labels=a.labels,
        mpi_ranks=a.mpi_ranks,require_reject=not a.allow_no_reject,window=a.window,nominal_dt=a.nominal_dt)
    (a.directory/'diagnostic-manifest.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['complete_for_declared_scope','errors','native_QPs','paired_assembly_tangent_calls']}))
    raise SystemExit(0 if r['complete_for_declared_scope'] else 1)
