#!/usr/bin/env python3
"""Read-only source comparison: expert C3D8R CSV vs computed HEX8 Exodus.

Writes analysis intermediates only in --output. Does not generate solver input or
fit material parameters. Uses original mesh centroids/connectivity to audit IDs,
not the rounded coordinate columns of the expert CSV. Missing frames stay missing.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import netCDF4
import numpy as np
import pandas as pd

FIELDS = {
    'DamageC': ('DAMAGEC', 1., '-'), 'DamageT': ('DAMAGET', 1., '-'),
    'stress_xx': ('S-S11', 1e-6, 'MPa'), 'stress_yy': ('S-S22', 1e-6, 'MPa'),
    'stress_zz': ('S-S33', 1e-6, 'MPa'), 'stress_xy': ('S-S12', 1e-6, 'MPa'),
    'stress_xz': ('S-S13', 1e-6, 'MPa'), 'stress_yz': ('S-S23', 1e-6, 'MPa'),
    'vonmises_stress': ('S-Mises', 1e-6, 'MPa'),
    'strain_xx': ('E-E11', 1., '-'), 'strain_yy': ('E-E22', 1., '-'),
    'strain_zz': ('E-E33', 1., '-'),
}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):h.update(b)
    return h.hexdigest()

def array(var):
    return np.asarray(np.ma.filled(var[:], np.nan), dtype=float)

def names(var):
    return [s.decode().strip('\x00 ') if isinstance(s,bytes) else str(s).strip('\x00 ')
            for s in netCDF4.chartostring(var[:])]

def topology(nc):
    xyz = np.column_stack([array(nc[v]) for v in ('coordx','coordy','coordz')])
    connections = [np.asarray(nc['connect'+str(b+1)][:], dtype=int)
                   for b in range(len(nc.dimensions['num_el_blk']))]
    centers = np.concatenate([xyz[c-1].mean(axis=1) for c in connections])
    return xyz, connections, centers

def element_order(mesh, result):
    _, input_conn, centers = topology(mesh)
    xyz, output_conn, output_centers = topology(result)
    ids = np.asarray(mesh['elem_num_map'][:], dtype=int)
    lookup = {tuple(np.round(c,12)):int(i) for c,i in zip(centers,ids)}
    if len(lookup)!=1000 or sum(len(c) for c in output_conn)!=1000:
        raise ValueError('Expected unique centroids for all 1000 source elements')
    mapped = np.asarray([lookup[tuple(np.round(c,12))] for c in output_centers])
    if len(set(mapped))!=1000:raise ValueError('Output topology does not map one-to-one')
    order = np.argsort(mapped)
    # Also check vertex sets; equal centroids alone are insufficient for a generic mesh.
    input_xyz = topology(mesh)[0]
    original = {int(i):set(tuple(np.round(p,12)) for p in input_xyz[c-1])
                for i,c in zip(ids,np.concatenate(input_conn))}
    for i,c in zip(mapped,np.concatenate(output_conn)):
        if original[int(i)] != set(tuple(np.round(p,12)) for p in xyz[c-1]):
            raise ValueError('Output vertices differ from original mesh')
    return order, mapped, xyz

def metrics(reference, computed, scale):
    diff = computed-reference
    # Near-zero reference values are reported by absolute error, not arbitrary division.
    eligible = np.abs(reference) > max(scale*1e-3, 1e-15)
    relative = np.abs(diff[eligible]/reference[eligible])*100
    return {
        'reference_min':float(reference.min()), 'reference_max':float(reference.max()),
        'reference_mean':float(reference.mean()), 'computed_min':float(computed.min()),
        'computed_max':float(computed.max()), 'computed_mean':float(computed.mean()),
        'mae':float(np.abs(diff).mean()), 'rmse':float(np.sqrt(np.mean(diff**2))),
        'max_abs_error':float(np.abs(diff).max()),
        'nrmse_global_peak_percent':float(np.sqrt(np.mean(diff**2))/scale*100) if scale else None,
        'max_error_global_peak_percent':float(np.abs(diff).max()/scale*100) if scale else None,
        'relative_eligible_elements':int(eligible.sum()),
        'within_5_percent_of_local_reference_fraction':float(np.mean(relative<=5)) if len(relative) else None,
        'max_local_relative_percent':float(relative.max()) if len(relative) else None,
        'worst_element':int(np.argmax(np.abs(diff))+1), 'global_reference_peak_scale':scale,
    }

def compare(reference, mesh, exodus, history, output):
    output.mkdir(parents=True, exist_ok=False)
    df = pd.read_csv(reference,skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df['time'] = df['Frame'].str.extract(r'Step Time\s*=\s*([\d.Ee+-]+)',expand=False).astype(float)
    if df.duplicated(['time','Element Label','IntPt']).any():raise ValueError('Duplicate reference rows')
    if set(df['IntPt'])!={1}:raise ValueError('Reference integration points changed')
    if not np.isfinite(df[[v[0] for v in FIELDS.values()]].to_numpy()).all():raise ValueError('Nonfinite reference')
    times = np.sort(df.time.unique())
    if len(df)!=101000 or len(times)!=101:raise ValueError('Unexpected reference frame/element count')
    scales = {key:float(df[col].abs().max()) for key,(col,_,_) in FIELDS.items()}
    summaries=[]; details=[]; matched=[]; missing=[]; qc={}
    with netCDF4.Dataset(mesh) as mesh_nc, netCDF4.Dataset(exodus) as nc:
        order, mapped, xyz = element_order(mesh_nc,nc)
        actual_times = array(nc['time_whole'])
        if np.any(np.diff(actual_times)<=0):raise ValueError('Nonmonotonic result frames')
        var_names=names(nc['name_elem_var'])
        data={}
        for key in FIELDS:
            idx=var_names.index(key)+1
            values=np.concatenate([array(nc[f'vals_elem_var{idx}eb{b+1}'])
                                  for b in range(len(nc.dimensions['num_el_blk']))],axis=1)
            if not np.isfinite(values).all():raise ValueError('Nonfinite output '+key)
            data[key]=values[:,order]*FIELDS[key][1]
        nodal_names=names(nc['name_nod_var'])
        uz=array(nc[f'vals_nod_var{nodal_names.index("disp_z")+1}'])
        if not np.isfinite(uz).all():raise ValueError('Nonfinite displacement')
        top=np.isclose(xyz[:,2],.2,atol=1e-10,rtol=0)
        bottom=np.isclose(xyz[:,2],.05,atol=1e-10,rtol=0)
        qc={'num_nodes':len(xyz), 'num_elements':len(mapped), 'element_vertex_mapping_verified':True,
            'output_elem_num_map_matches_source':bool(np.array_equal(nc['elem_num_map'][:], mapped)),
            'frames_in_exodus':len(actual_times),'last_output_time':float(actual_times[-1]),
            'last_top_displacement_m':float(uz[-1,top].mean()),
            'top_loading_max_abs_error_m':float(np.abs(uz[:,top]-actual_times[:,None]*2.5e-5).max()),
            'bottom_z_max_abs_m':float(np.abs(uz[:,bottom]).max()),
            'damage_range':{k:[float(data[k].min()),float(data[k].max())] for k in ('DamageC','DamageT')},
            'mean_E33_ramp_max_abs_error':float(np.abs(data['strain_zz'].mean(1)-actual_times*2.5e-5/.15).max()),
            'all_compared_values_finite':True}
        for time in times:
            locations=np.flatnonzero(np.isclose(actual_times,time,atol=1e-8,rtol=0))
            if len(locations)!=1:
                missing.append(float(time));continue
            index=locations[0];matched.append(float(time))
            ref=df[df.time==time].sort_values('Element Label')
            if not np.array_equal(ref['Element Label'],np.arange(1,1001)):raise ValueError('Missing reference elements')
            frame=pd.DataFrame({'time_s':time,'element_id':np.arange(1,1001)})
            for key,(col,_,unit) in FIELDS.items():
                r=ref[col].to_numpy();c=data[key][index]
                summaries.append({'time_s':float(time),'field':key,'unit':unit,**metrics(r,c,scales[key])})
                frame[key+'_abaqus']=r;frame[key+'_moose']=c;frame[key+'_error']=c-r
            details.append(frame)
    hist=pd.read_csv(history)
    if not np.isfinite(hist.to_numpy(dtype=float)).all():raise ValueError('Nonfinite solver history')
    peak=float(hist.RP1_Force.abs().max())
    qc.update({'max_axial_reaction_imbalance_N':float((hist.RP1_Force+hist.Bottom_Force).abs().max()),
               'axial_balance_percent_of_peak_force':float((hist.RP1_Force+hist.Bottom_Force).abs().max()/peak*100) if peak else None,
               'max_transverse_gauge_reaction_N':float(hist[['Top_Force_X','Top_Force_Y']].abs().to_numpy().max()),
               'max_axial_force_N':peak, 'history_last_time':float(hist.time.iloc[-1])})
    pd.DataFrame(summaries).to_csv(output/'frame_metrics.csv',index=False)
    pd.concat(details,ignore_index=True).to_csv(output/'element_comparison.csv',index=False,float_format='%.12g')
    # A small independent long-form summary is easier to inspect than the full cell table.
    report={'status':'diagnostic-only', 'matched_frames':len(matched),'matched_times':matched,
        'missing_reference_times':missing, 'quality_checks':qc,
        'method':{'identity':'original INP element ID verified against mesh vertices/centroids',
            'fields':'C3D8R IntPt=1 versus HEX8 quadrature-projected element constants',
            'stress_unit':'MPa (MOOSE Pa / 1e6)', 'time_matching':'exact within 1e-8; no result interpolation',
            'normalization':'Each field uses max(abs(reference)) over all 101 reference frames',
            'pointwise_relative_mask':'abs(reference)>0.001*field_global_peak; otherwise absolute errors only',
            'limitations':['Diagnostic: differing element integration/hourglass formulations',
                'RF3/U3 not present in reference CSV', 'Shear strains excluded pending tensor/engineering convention audit',
                'Passing a global-normalized metric is not the same as every element meeting ±5%']},
        'sources':{str(p):sha(p) for p in (reference,mesh,exodus,history)},
        'key_frames':[x for x in summaries if x['field'] in ('DamageT','DamageC','stress_zz','vonmises_stress')
                      and any(abs(x['time_s']-t)<1e-8 for t in (.01,.1,.3,.35,.4,.44,.5,.75,1))]}
    (output/'comparison-summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'matched_frames':len(matched),'missing_frames':len(missing),'qc':qc},ensure_ascii=False,indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('reference','mesh','exodus','history','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();compare(a.reference,a.mesh,a.exodus,a.history,a.output)
