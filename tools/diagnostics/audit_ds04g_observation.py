#!/usr/bin/env python3
"""Audit passive traces; never run a material integrator or submit a simulation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def audit(directory, prefix, job_id=None, development=False):
    statuses=[json.loads(p.read_text()) for p in sorted(directory.glob(prefix+'-status-rank*.json'))]
    errors=[]
    if sorted(s['rank'] for s in statuses)!=[0,1,2,3]:errors.append('missing_or_duplicate_MPI_status')
    for s in statuses:
        for k in ['capture_truncated','io_failed','history_incomplete']:
            if s[k]:errors.append('rank%d:%s'%(s['rank'],k))
    points=set();owners={};sizes={};events={};weights={};matrix_calls=set();raw_calls=set()
    files=sorted(directory.glob(prefix+'-*.csv'))
    for p in files:
        rank=int(p.stem.rsplit('rank',1)[1]); count=0
        for r in csv.DictReader(p.open(newline='')):
            count+=1; event=r['event'];events[event]=events.get(event,0)+1
            label=int(r['source_element']);qp=int(r['qp'])
            if event=='qp' and r['field']=='x_y_z_JxW_coord':
                points.add((label,qp)); owners.setdefault(label,set()).add(rank)
                v=[float(x) for x in r['values'].split()]
                weights[label,qp]=v[3]*v[4]
            if r['field']=='column_major_6x6':
                raw_calls.add((rank,int(r['evaluation']),int(r['call'])))
                if len(r['values'].split())!=36:errors.append('bad_6x6_shape')
            if r['field']=='RankFour_ijkl':
                matrix_calls.add((rank,int(r['evaluation']),int(r['call'])))
                if len(r['values'].split())!=81:errors.append('bad_RankFour_shape')
            if r['field']=='row_major_9x9' and len(r['values'].split())!=81:errors.append('bad_9x9_shape')
        sizes[p.name]={'rows':count,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    expected={(e,q) for e in [1881,767] for q in range(4)}
    if points!=expected:errors.append('expected_two_times_four_native_QPs')
    if any(len(x)!=1 for x in owners.values()):errors.append('duplicate_owner')
    if any(w<=0 for w in weights.values()):errors.append('nonpositive_quadrature_weight')
    if matrix_calls-raw_calls:errors.append('assembly_tangent_has_no_matching_material_call')
    if not development:
        if not any('local-jacobian' in p.name for p in files):errors.append('missing_plastic_local_Jacobian')
        if not matrix_calls:errors.append('missing_actual_assembly_tangent')
        if not events.get('reject_after'):errors.append('no_rejected_step_restore_observed')
        if not events.get('accepted_end'):errors.append('no_accepted_step_observed')
    result={'schema':'ds04g-observation-audit-v1','job_id':job_id,'development_only':development,
            'complete_for_declared_scope':not errors,'errors':sorted(set(errors)),
            'statuses':statuses,'native_QPs':sorted(points),'owners':{k:sorted(v) for k,v in owners.items()},
            'normalized_weights':{str(e):[weights.get((e,q),0)/sum(weights.get((e,j),0) for j in range(4)) for q in range(4)] for e in owners if sum(weights.get((e,j),0) for j in range(4))>0},
            'events':events,'files':sizes,'paired_assembly_tangent_calls':len(matrix_calls & raw_calls),
            'scope_note':'Structural trace integrity only; does not prove algorithmic tangent accuracy, global rollback, unchanged trajectory or engineering acceptance'}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--prefix',default='ds04g-observation')
    group=p.add_mutually_exclusive_group(required=True);group.add_argument('--job-id');group.add_argument('--development',action='store_true')
    a=p.parse_args();r=audit(a.directory,a.prefix,a.job_id,a.development)
    (a.directory/'diagnostic-manifest.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['complete_for_declared_scope','errors','native_QPs','paired_assembly_tangent_calls']}))
    raise SystemExit(0 if r['complete_for_declared_scope'] else 1)
