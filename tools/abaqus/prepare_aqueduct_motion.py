#!/usr/bin/env python3
"""Source-hashed AQ1023 baseline-corrected, SI three-axis motion (not a solver run).

Implements the documented Newmark baseline matrix on exact piecewise-linear raw
acceleration integrals, then samples the quadratic correction at source knots.
The corrected acceleration is piecewise linear, as in the documented procedure.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as P

from prepare_aqueduct_case import SOURCE_SHA


def integrate(times,accel):
    velocity=np.zeros(len(times));displacement=np.zeros(len(times))
    for j,h in enumerate(np.diff(times)):
        slope=(accel[j+1]-accel[j])/h
        velocity[j+1]=velocity[j]+accel[j]*h+slope*h*h/2
        displacement[j+1]=displacement[j]+velocity[j]*h+accel[j]*h*h/2+slope*h*h*h/6
    return velocity,displacement


def baseline(times,accel,start=0.):
    times=np.asarray(times,dtype=float);accel=np.asarray(accel,dtype=float)
    if len(times)<2 or len(times)!=len(accel) or not np.all(np.isfinite(times)) or not np.all(np.isfinite(accel)) or np.any(np.diff(times)<=0):
        raise ValueError('finite strictly increasing complete acceleration table required')
    if start!=times[0] or times[-1]<=start:raise ValueError('explicit correction interval must span the supplied knots')
    duration=times[-1]-start;q=(times-start)/duration;v,_=integrate(times,accel);moments=np.zeros(3)
    # Integrate in segment-local coordinates to avoid subtracting almost equal
    # high-degree polynomial values at the dense earthquake knots.
    h=np.diff(q);left=q[:-1];slope=np.diff(accel)/h
    for k in range(1,4):
        for j in range(k+1):
            primitive=v[:-1]*h**(j+1)/(j+1)+duration*accel[:-1]*h**(j+2)/(j+2)+duration*slope*h**(j+3)/(2*(j+3))
            moments[k-1]+=float(np.sum(math.comb(k,j)*left**(k-j)*primitive))/duration
    matrix=np.array([[-300.,900.,-630.],[1800.,-5760.,4200.],[-1890.,6300.,-4725.]])
    coeff=matrix@moments;corrected=accel+P(coeff)(q)
    cv,cu=integrate(times,corrected)
    return dict(coefficients=coeff.tolist(),moments= moments.tolist(),corrected=corrected,velocity=cv,displacement=cu)


def prepare(contract,output):
    path=Path(contract);d=json.loads(path.read_text());out=Path(output)
    if d['source']['sha256']!=SOURCE_SHA or d['case_id']!='AQ1023-FULL-ELC':raise ValueError('wrong approved source')
    if out.exists():raise ValueError('output exists; preserve previous preparation')
    dynamic=next(s for s in d['steps_raw'] if s['name']=='ELC')
    amplitude=d['amplitudes_raw'];scale={b['amplitude']:b['value']*1e-3 for b in dynamic['boundaries'] if b['type']=='ACCELERATION'}
    if scale!={'X':3.92,'Y':3.332,'Z':2.548}:raise ValueError('source acceleration scale changed')
    corrections=[c for c in d['interfaces']['source_cards'] if c['keyword']=='baseline correction']
    if len(corrections)!=3 or any(c['data'] for c in corrections):raise ValueError('only source whole-record correction supported')
    out.mkdir(parents=True);summary=dict(source_inp_sha256=SOURCE_SHA,contract_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),static_duration_s=1.,dynamic_duration_s=dynamic['dynamic'][1],method='Documented Newmark baseline matrix; exact raw PWL velocity moments, correction sampled at source knots; corrected PWL acceleration',conventions={'record_start_s':0.,'before_first_source_knot':'hold first raw value back to t=0 before baseline correction','correction_end':'last original amplitude knot, not padded analysis end','after_last_source_knot':'hold last corrected acceleration to analysis end','step_time_to_total_time':'total=1+ELC step time','initial_dynamic_acceleration':'zero/bypass candidate, resolved HHT mode still needs 2022 check'},reference_status='custom-documented-method; exact Abaqus2022 corrected samples not supplied',axes={})
    for axis in 'XYZ':
        raw=np.asarray(amplitude[axis],float);t=raw[:,0];a=raw[:,1]*scale[axis]
        if t[0]<=0 or t[-1]>=dynamic['dynamic'][1]:raise ValueError('source domain changed')
        times=np.r_[0.,t];accel=np.r_[a[0],a];r=baseline(times,accel)
        final_times=np.r_[times,dynamic['dynamic'][1]];final_accel=np.r_[r['corrected'],r['corrected'][-1]]
        np.savetxt(out/('motion-'+axis.lower()+'.csv'),np.column_stack([final_times,final_accel]),delimiter=',',fmt='%.17g')
        np.savetxt(out/('raw-'+axis.lower()+'.csv'),np.column_stack([times,accel]),delimiter=',',fmt='%.17g')
        cv,cu=integrate(final_times,final_accel);rv,ru=integrate(times,accel)
        # Record a first-knot-start sensitivity, rather than silently asserting Abaqus interval resolution.
        alternative=baseline(t,a,start=float(t[0]));alt=alternative['corrected'];delta=float(np.max(np.abs(r['corrected'][1:]-alt)))
        summary['axes'][axis]=dict(source_knots=len(t),first_source_time_s=float(t[0]),last_source_time_s=float(t[-1]),scale_m_s2=scale[axis],correction_coefficients_m_s2=r['coefficients'],raw_final_velocity_m_s=float(rv[-1]),raw_final_displacement_m=float(ru[-1]),corrected_final_velocity_m_s=float(cv[-1]),corrected_final_displacement_m=float(cu[-1]),corrected_initial_acceleration_m_s2=float(final_accel[0]),corrected_final_acceleration_m_s2=float(final_accel[-1]),alternative_interval_start_s=float(t[0]),alternative_interval_max_acceleration_difference_m_s2=delta)
    summary['artifacts']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.glob('*.csv'))}
    (out/'motion-identity.json').write_text(json.dumps(summary,indent=2)+'\n');return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contract',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps(prepare(a.contract,a.output),indent=2))
