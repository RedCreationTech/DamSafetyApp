#!/usr/bin/env python3
"""Compare a short original-mesh Job with the sealed 18-field reference pairing.

Requires numpy/netCDF4. Preserves the full-run reference normalization and the
original near-zero eligibility rule. Does not infer unsaved IPs or later times.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import netCDF4
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(root, baseline):
    meta=json.loads((baseline/'comparison/comparison_metadata.json').read_text())
    fields=list(meta['paired_field_map'])
    scales={r['field']:float(r['global_reference_peak_scale']) for r in csv.DictReader((baseline/'comparison/field_error_curves.csv').open())}
    reference_path=baseline/'comparison/paired_all_18_fields.csv.gz'
    with netCDF4.Dataset(root/'results/uniaxial_compression.e') as ds, netCDF4.Dataset(baseline/'results/uniaxial_compression.e') as old:
        for key in ['coordx','coordy','coordz','connect1','elem_num_map','node_num_map']:
            assert np.array_equal(ds[key][:],old[key][:]),key
        ids=np.asarray(ds['elem_num_map'][:],int)
        assert np.array_equal(ids,np.arange(1,1001))
        names=list(netCDF4.chartostring(ds['name_elem_var'][:]))
        times=np.asarray(ds['time_whole'][:]); old_times=np.asarray(old['time_whole'][:])
        old_names=list(netCDF4.chartostring(old['name_elem_var'][:]))
        refs={round(float(t),8):{k:np.full(1000,np.nan) for k in fields} for t in times}
        with gzip.open(reference_path,'rt') as f:
            for r in csv.DictReader(f):
                t=round(float(r['time_s']),8)
                if t>times[-1]+1e-8: break
                if t not in refs: continue
                element=int(r['element_id'])-1
                for key in fields: refs[t][key][element]=float(r[key+'_abaqus'])
        rows=[]; old_diffs={}
        for i,t in enumerate(times):
            assert abs(old_times[np.argmin(abs(old_times-t))]-t)<1e-8
            j=int(np.argmin(abs(old_times-t)))
            for key in fields:
                scale=meta['paired_field_map'][key][1]
                c=np.asarray(ds[f'vals_elem_var{names.index(key)+1}eb1'][i,:])*scale
                previous=np.asarray(old[f'vals_elem_var{old_names.index(key)+1}eb1'][j,:])*scale
                r=refs[round(float(t),8)][key]
                assert np.isfinite(c).all() and np.isfinite(r).all()
                d=c-r; mask=abs(r)>max(scales[key]*1e-3,1e-15)
                old_diffs[key]=max(old_diffs.get(key,0.),float(np.max(abs(c-previous))))
                rows.append(dict(time_s=float(t),field=key,reference_mean=float(r.mean()),
                    computed_mean=float(c.mean()),reference_max=float(r.max()),computed_max=float(c.max()),
                    mae=float(abs(d).mean()),rmse=float(np.sqrt(np.mean(d*d))),
                    max_abs_error=float(abs(d).max()),global_reference_peak_scale=scales[key],
                    nrmse_global_peak_percent=float(np.sqrt(np.mean(d*d))/scales[key]*100) if scales[key] else None,
                    relative_eligible_elements=int(mask.sum()),
                    within_5_percent_fraction=float(np.mean(abs(d[mask]/r[mask])<=.05)) if mask.any() else None,
                    worst_abaqus_element_id=int(ids[np.argmax(abs(d))])))
    with (root/'short-field-errors.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={'job_id':json.loads((root/'manifest.json').read_text())['job_id'],
        'scope':'Original 1000 elements, 18 element-average fields, saved common times only; no later-time or raw-IP equivalence claim',
        'matched_frames':len(times),'end_time':float(times[-1]),'topology_and_element_ids_match':True,
        'reference_pairing_sha256':sha(reference_path),
        'reference_metadata_sha256':sha(baseline/'comparison/comparison_metadata.json'),
        'computed_exodus_sha256':sha(root/'results/uniaxial_compression.e'),
        'max_difference_from_original_moose_by_field':old_diffs,
        'snapshots':[r for r in rows if r['field'] in ('DamageT','DamageC','stress_zz','vonmises_stress') and any(abs(r['time_s']-t)<1e-8 for t in (.02,.05,.1))]}
    (root/'short-comparison-summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--baseline',type=Path,required=True);a=p.parse_args()
    result=analyze(a.root,a.baseline)
    print(json.dumps({'job_id':result['job_id'],'frames':result['matched_frames'],
        'final':[r for r in result['snapshots'] if abs(r['time_s']-.1)<1e-8]},ensure_ascii=False))
