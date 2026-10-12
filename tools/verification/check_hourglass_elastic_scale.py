#!/usr/bin/env python3
"""Check candidate C++ source expressions offline; no build or FE execution.

This checks algebra/unchanged source, not compiled MOOSE AD or input parsing.
"""
import argparse
import ast
import csv
import hashlib
import itertools
import json
import math
import random
import re
import struct
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = 'bc5425dde7ccdb1b4c4688f2b92ed0fd8625a138'
CPP = 'src/kernels/ADQuad4HourglassControl.C'
HEADER = 'include/kernels/ADQuad4HourglassControl.h'
PARAM = '''  params.addRangeCheckedParam<Real>(
      "elastic_stiffness_scale", 1.0, "elastic_stiffness_scale > 0",
      "Scales only the elastic hourglass force; additional damping retains its original stiffness");
'''
CORRECTION = '''  // Keep the default evaluation order and the entire historical damping term unchanged.
  if (_elastic_stiffness_scale != 1.0)
    amplitude += HourglassElasticScale::correction(
        _elastic_stiffness_scale, _alpha, _dt, unweighted_amplitude, amplitude_old);

'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expression(text):
    text=' '.join(text.split())
    replacements = {'_zeta[_qp]':'zeta', '_elastic_stiffness_scale':'scale',
                    'unweighted_amplitude':'current', 'amplitude_older':'older',
                    'amplitude_old':'old', 'elastic_amplitude':'elastic',
                    'amplitude':'current', '_alpha':'alpha', '_dt':'dt'}
    for old,new in replacements.items():text=text.replace(old,new)
    tree=ast.parse(text.strip(),mode='eval')
    allowed=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Add,ast.Sub,ast.Mult,ast.Div,
             ast.USub,ast.UAdd,ast.Name,ast.Load,ast.Constant)
    assert all(isinstance(n,allowed) for n in ast.walk(tree))
    assert {n.id for n in ast.walk(tree) if isinstance(n,ast.Name)} <= {
        'current','old','older','alpha','zeta','dt','scale','elastic'}
    code=compile(tree,'<locked C++ arithmetic>','eval')
    return lambda values:eval(code,{'__builtins__':{}},values)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    p.add_argument('--native-checkpoints',required=True,type=Path)
    p.add_argument('--native-checksums',required=True,type=Path)
    args=p.parse_args();assert not args.output.exists(), 'Preserve completed check'
    original=subprocess.check_output(['git','-C',str(ROOT),'show',BASE+':'+CPP],text=True)
    current=(ROOT/CPP).read_text()
    stripped=current
    additions=(PARAM,'#include "HourglassElasticScale.h"\n','    _elastic_stiffness_scale(getParam<Real>("elastic_stiffness_scale")),\n',
               '  const ADReal unweighted_amplitude = amplitude;\n',CORRECTION)
    for part in additions:
        assert stripped.count(part)==1
        stripped=stripped.replace(part,'',1)
    assert stripped==original, 'Existing kernel, damping and correction evaluation order must remain exact'
    old_header=subprocess.check_output(['git','-C',str(ROOT),'show',BASE+':'+HEADER],text=True)
    assert (ROOT/HEADER).read_text().replace('  const Real _elastic_stiffness_scale;\n','',1)==old_header
    unchanged=['include/utils/Quad4Hourglass.h','include/utils/StressHistoryInterval.h',
               'src/materials/CDPVelocityRateMaterial.C','src/materials/AbaqusCDPStressUpdate.C']
    for name in unchanged:
        assert (ROOT/name).read_bytes()==subprocess.check_output(['git','-C',str(ROOT),'show',BASE+':'+name])
    legacy=expression(re.search(r'amplitude = (amplitude \* .*?);',original,re.S)[1])
    helper=(ROOT/'include/utils/HourglassElasticScale.h').read_text()
    delta=expression(re.search(r'if \(dt > 0\)\s*return (.*?);',helper,re.S)[1])
    static_delta=expression(re.findall(r'return (.*?);',helper,re.S)[-1])
    def legacy_residual(q):
        if q['dt']>0:
            value=legacy(q)
            if q.get('previous_interval',False):
                value+=q['alpha']*q['zeta']*(1/q['previous_dt']-1/q['dt'])*(q['older']-q['old'])
        else:value=q['current']
        return value
    def residual(q,scale):
        value=legacy_residual(q)
        if scale!=1.:
            value+=(delta if q['dt']>0 else static_delta)(dict(q,scale=scale))
        return value
    rng=random.Random(20261012);count=0;max_error=0.;max_damping_difference=0.
    for dt,alpha,zeta,paired,scale in itertools.product(
            (0.,-.01,5e-7,.000234375,.001,.005,.01),(-.05,0.),(0.,.00113),
            (False,True),(.5,.8,1.,1.2)):
        for state in range(12):
            q=dict(current=rng.uniform(-.001,.001),old=rng.uniform(-.001,.001),
                   older=rng.uniform(-.001,.001),dt=dt,alpha=alpha,zeta=zeta,
                   previous_dt=rng.choice((.00005,.003,.01)),previous_interval=paired)
            baseline=legacy_residual(q)
            # Default branch performs the unchanged arithmetic, including bit pattern.
            assert struct.pack('d',baseline)==struct.pack('d',residual(q,1.))
            k=(1+alpha)*q['current']-alpha*q['old'] if dt>0 else q['current']
            expected=baseline+(scale-1)*k
            error=abs(residual(q,scale)-expected);max_error=max(max_error,error)
            assert error<=1e-13*max(1.,abs(expected))
            c0=baseline-k;c1=residual(q,scale)-scale*k
            diff=abs(c1-c0);max_damping_difference=max(max_damping_difference,diff)
            assert diff<=1e-13*max(1.,abs(c0))
            # Linear operator coefficients from actual source expressions; no FD cancellation.
            zero=dict(q,current=0.,old=0.,older=0.)
            for name,want in [('current',(scale-1)*(1+alpha) if dt>0 else scale-1),
                              ('old',-(scale-1)*alpha if dt>0 else 0.),('older',0.)]:
                basis=dict(zero,**{name:1.})
                change=residual(basis,scale)-residual(basis,1.)
                assert abs(change-want)<=1e-10*max(1.,abs(want))
            count+=1
    # Actual stored native q values: static elastic force scales, damping power stays fixed.
    checks=json.loads(args.native_checksums.read_text())
    assert sha(args.native_checkpoints)==checks[args.native_checkpoints.name]
    samples=0;max_static_error=0.;max_power_error=0.
    with args.native_checkpoints.open(newline='') as f:
        for row in csv.DictReader(f):
            if row['round']!='H40':continue
            stiffness=float(row['modal_stiffness_N_per_m2'])
            q=[float(row['qx_m']),float(row['qy_m'])]
            rate=[float(row['qdotx_Newmark_m_per_s']),float(row['qdoty_Newmark_m_per_s'])]
            f0=stiffness*math.hypot(*q)
            fs=stiffness*math.hypot(*(residual(dict(current=x,old=0.,older=0.,dt=0.),.8) for x in q))
            error=abs(fs-.8*f0);max_static_error=max(max_static_error,error)
            assert error<1e-9*max(1.,f0)
            power0=.00113*stiffness*sum(x*x for x in rate)
            probe=dict(current=1.,old=0.,older=0.,dt=.01,alpha=-.05,zeta=.00113)
            damping_current_coefficient=residual(probe,.8)-.8*(1+probe['alpha'])
            effective_zeta=damping_current_coefficient*probe['dt']/(1+probe['alpha'])
            power_candidate=effective_zeta*stiffness*sum(x*x for x in rate)
            max_power_error=max(max_power_error,abs(power0-power_candidate))
            assert abs(power0-power_candidate)<1e-12*max(1.,power0)
            samples+=1
    result=dict(status='passed-offline',source_base_commit=BASE,
        source_branch=subprocess.check_output(['git','-C',str(ROOT),'branch','--show-current'],text=True).strip(),
        source_worktree=str(ROOT),committed=False,compiled_MOOSE_or_AD_runtime_verified=False,
        no_solver_or_network_calls=True,default_original_source_reconstruction_exact=True,
        default_scale=1.,planned_trial_scale=.8,modal_test_states=count,
        maximum_residual_delta_error=max_error,maximum_extracted_damping_roundoff=max_damping_difference,
        elastic_current_and_history_coefficient_changes_verified=True,
        damping_current_old_older_coefficients_unchanged=True,native_H40_checkpoint_states=samples,
        native_static_force_scale_max_absolute_error_N_per_m=max_static_error,
        frozen_path_candidate_power_change_W_per_m=max_power_error,
        actual_old_native_HHT_forces_reconstructed=False,
        unchanged_source_files={name:sha(ROOT/name) for name in unchanged},
        candidate_source_sha256={name:sha(ROOT/name) for name in (CPP,HEADER,'include/utils/HourglassElasticScale.h',
            'unit/src/HourglassElasticScaleTest.C','tools/verification/check_hourglass_elastic_scale.py')},
        native_checkpoints_sha256=sha(args.native_checkpoints),
        pending=['Locked-node C++ build and actual AD assembly/unit checks','Real mesh --check-input',
                 'Same candidate scale1/scale0.8 25s FE jobs if user authorizes','Seven frozen cases before promotion'])
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
