#!/usr/bin/env python3
"""Create an immutable, single-change CDP diagnostic snapshot; never run a solver.

Start with `observe` against D02, then derive one numerical/boundary probe at a
time. Native MOOSE scalar extraction samples accepted IPs, not projected fields
or failed local iterates. Every material table remains byte-identical.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def observations():
    text = ''
    properties = ['DamageT', 'DamageC', 'cdp_kappa_t', 'cdp_kappa_c',
                  'cdp_combined_damage', 'cdp_stiffness_factor']
    properties += ['cdp_'+name for name in (
        'local_iterations', 'accepted_substeps', 'jacobian_fallbacks',
        'failed_material_calls', 'attempted_partitions', 'maximum_partition_depth',
        'automatic_jacobian_evaluations', 'finite_difference_jacobian_evaluations',
        'local_factorizations', 'local_backsolves', 'integration_microseconds')]
    for tensor, prefix in [('stress', 'sigma'), ('total_strain', 'eps'),
                           ('cdp_backbone_plastic_strain', 'pb'),
                           ('cdp_viscous_plastic_strain', 'pv')]:
        for i, j, label in [(0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz'),
                            (0, 1, 'xy'), (0, 2, 'xz'), (1, 2, 'yz')]:
            name = 'probe_'+prefix+'_'+label
            properties.append(name)
            text += f'''[Materials/{name}]
  type = RankTwoCartesianComponent
  rank_two_tensor = {tensor}
  property_name = {name}
  index_i = {i}
  index_j = {j}
[]
'''
    return text+f'''[VectorPostprocessors/ip_states]
  type = ElementMaterialSampler
  property = '{' '.join(properties)}'
  execute_on = 'initial timestep_end'
[]
[Outputs/ip_csv]
  type = CSV
  file_base = ip_history
  show = ip_states
  execute_on = 'initial timestep_end'
  time_data = true
  precision = 17
[]
'''


def build(parent, output, variant, case_name):
    if (parent/'submission.json').exists():
        old = json.loads((parent/'submission.json').read_text())
    else:
        baseline = json.loads((parent/'manifest.json').read_text())
        old = dict(baseline['input_snapshot'],project_id='gmp-ise',case_name=baseline['case_name'])
    src = parent/'input'/old['input_file']
    assert sha(src) == old['input_sha256']
    text = src.read_text()
    extra_evidence = None
    extra_evidence_name = 'replay-targets.json'
    if variant == 'observe':
        assert '[VectorPostprocessors/ip_states]' not in text
        text = text.replace('[Executioner]\n',
            "[Executioner]\n  petsc_options = '-snes_converged_reason'\n", 1)
        text += '\n'+observations()
        change = 'Only add accepted-IP scalar/tensor CSV and PETSc convergence reasons'
    elif variant == 'full-short-observe':
        assert old['input_sha256']=='f8f558b93ad674192ca51f4e7fae0ce85eeb774d9127f00fdeb095e478034e9a'
        text=text.replace('[Executioner]\n',"[Executioner]\n  petsc_options = '-snes_converged_reason'\n",1)
        text,n=re.subn(r'(?m)^(\s+)end_time = 1$',r'\1end_time = 0.1',text);assert n==2
        observe=observations().replace('type = ElementMaterialSampler\n',
            "type = ElementMaterialSampler\n  elem_ids = '0 291 401 410 460 496 555 791 896'\n")
        text='# Diagnostic observation window 0-0.1s; original loading RATE and boundary conditions retained.\n'+text+'\n'+observe
        change='Original 1000 HEX8 / original BC / original numerics: only shorten to 0.1s and observe 9 fixed MOOSE element indices at all IPs'
    elif variant == 'dual-ip-path-trace':
        # Keep the formal 1000-HEX8 behavior and original adaptive stepping,
        # but bound the diagnostic to the first-divergence window and the two
        # physical points whose full-FE DamageT errors have opposite signs.
        text, n = re.subn(r'(?m)^(\s+)end_time = 0\.1$', r'\1end_time = 0.03', text)
        assert n == 2  # Executioner and output Times; loading rate is unchanged.
        text, n = re.subn(
            r"(?m)^  elem_ids = '[^']+'$",
            "  elem_ids = '54 94'",
            text,
        )
        assert n == 1
        marker = '    enable_performance_diagnostics = true\n'
        assert text.count(marker) == 1
        text = text.replace(
            marker,
            marker
            + '    enable_path_diagnostics = true\n'
            + "    diagnostic_trace_elements = '54 94'\n"
            + '    diagnostic_time_begin = 0.01\n'
            + '    diagnostic_time_end = 0.03\n'
            + '    diagnostic_max_trace_calls = 2000\n',
        )
        text = (
            '# Dual-IP accepted-path trace: Abaqus 55/IP7 -> MOOSE 54/qp6; '
            'Abaqus 95/IP1 -> MOOSE 94/qp0.\n'
            + text
        )
        extra_evidence = {
            'parent_job': json.loads((parent/'manifest.json').read_text())['job_id'],
            'window_s': [0.01, 0.03],
            'targets': [
                {'abaqus_element': 55, 'abaqus_ip': 7, 'moose_element': 54, 'moose_qp': 6,
                 'full_fe_error_sign': 'positive'},
                {'abaqus_element': 95, 'abaqus_ip': 1, 'moose_element': 94, 'moose_qp': 0,
                 'full_fe_error_sign': 'negative'},
            ],
            'purpose': 'Locate the first accepted-state strain/history divergence for two opposite-sign DamageT points; observation only',
        }
        extra_evidence_name = 'dual-ip-targets.json'
        change = ('Only shorten the observation horizon to 0.03 s and enable bounded path diagnostics '
                  'for MOOSE elements 54/94; material, mesh, boundary, loading rate, tolerances and '
                  'IterationAdaptiveDT behavior remain unchanged')
    elif variant == 'tight-abs':
        text, n = re.subn(r'(?m)^  nl_abs_tol = 1e-8$', '  nl_abs_tol = 1e-13', text)
        assert n == 1
        change = 'Only tighten global scaled nl_abs_tol from 1e-8 to 1e-13'
    elif variant == 'no-scaling':
        assert text.count('automatic_scaling = true') == 1
        text = text.replace('automatic_scaling = true', 'automatic_scaling = false')
        change = 'Only disable automatic scaling; keep all convergence tolerances'
    elif variant == 'top-clamp':
        text += '''\n[BCs/top_x_clamp]
  type = DirichletBC
  variable = disp_x
  boundary = front
  value = 0
[]
[BCs/top_y_clamp]
  type = DirichletBC
  variable = disp_y
  boundary = front
  value = 0
[]
'''
        change = 'Only add x/y clamps on loaded z-max face on the same mesh'
    elif variant == 'early-window':
        text, n = re.subn(r'(?m)^(\s+)end_time = 1$', r'\1end_time = 0.05', text)
        assert n == 2  # Executioner and output Times; loading remains -0.0025*t
        change = 'Only shorten observation to 0.05s; retain loading rate and 0.01s output times'
    elif variant in ('dt001', 'dt0001'):
        dt = '0.001' if variant == 'dt001' else '0.0001'
        text, n = re.subn(r'(?m)^(\s+)dt = [0-9.]+$', lambda m:m[1]+'dt = '+dt, text)
        assert n == 1
        text, n = re.subn(r'(?m)^(\s+)dtmax = [0-9.]+$', lambda m:m[1]+'dtmax = '+dt, text)
        assert n == 1
        change = 'Only cap initial and maximum global dt at '+dt+'s; retain adaptive cutback'
    elif variant == 'early-refined-grid':
        start,end=text.index('  [TimeStepper]'),text.index('\n[Times/field_output_times]')
        end_time=float(re.search(r'(?m)^  end_time = ([0-9.]+)$',text).group(1))
        sequence=[i/1000 for i in range(1,31)]+[i/100 for i in range(4,round(end_time*100)+1)]
        text=text[:start]+"  [TimeStepper]\n    type = TimeSequenceStepper\n    time_sequence = '"+' '.join(f'{t:.3f}' for t in sequence)+"'\n  []\n[]\n"+text[end:]
        change='Only use a prescribed time grid: dt=0.001 through 0.03s, then dt=0.01, with cutbacks retained; diagnostic, not a production performance claim'
    elif variant == 'replay-path':
        # The early-window parent contains the two accepted frames spanning the
        # first kappa_t increase. Prescribe that material-point strain segment
        # as affine displacements; expose its eight local targets as FE steps.
        targets = []
        for step in (1, 2):
            path = parent/'results'/f'ip_history_ip_states_{step:04d}.csv'
            rows = list(csv.DictReader(path.open()))
            row = next(r for r in rows if r['elem_id']=='0' and r['qp_id']=='0')
            targets.append({k:float(v) for k,v in row.items()})
        time_rows = list(csv.DictReader((parent/'results/ip_history_ip_states_time.csv').open()))
        assert float(time_rows[1]['time']) == .01 and float(time_rows[2]['time']) == .02
        assert targets[0]['cdp_kappa_t'] == targets[0]['cdp_kappa_c'] == 0
        assert targets[1]['cdp_accepted_substeps'] == 8
        functions = ''
        for axis in ('x','y','z'):
            pieces = []
            for coord in ('x','y','z'):
                label = ''.join(sorted(axis+coord))
                old_value,new_value=[r['probe_eps_'+label] for r in targets]
                pieces.append(f'if(t<=0.01,{old_value:.17g}*t/0.01,({old_value:.17g})+(t-0.01)*({new_value-old_value:.17g})/0.01)*{coord}')
            functions += f"[Functions/replay_{axis}]\n  type = ParsedFunction\n  expression = '{'+'.join(pieces)}'\n[]\n"
            functions += f"[BCs/replay_{axis}]\n  type = FunctionDirichletBC\n  variable = disp_{axis}\n  boundary = 'left right bottom top back front'\n  function = replay_{axis}\n[]\n"
        start,end=text.index('[Functions/load_path]'),text.index('[Preconditioning/smp]')
        text=text[:start]+functions+text[end:]
        text,n=re.subn(r'(?m)^(\s+)end_time = 0.05$',r'\1end_time = 0.02',text);assert n==2
        start,end=text.index('  [TimeStepper]'),text.index('\n[Times/field_output_times]')
        text=text[:start]+'''  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '0.01 0.01125 0.0125 0.01375 0.015 0.01625 0.0175 0.01875 0.02'
  []
[]
'''+text[end:]
        extra_evidence={'parent_job':json.loads((parent/'manifest.json').read_text())['job_id'],
            'element':0,'qp':0,'times':[.01,.02],'raw_targets':targets,
            'purpose':'Replay the observed affine strain segment and expose 8 substeps; not a free-lateral boundary case or reference FE solution'}
        change='Controlled material-path replay: prescribe the observed 0.01-0.02s strain segment and expose 8 equal local targets as global steps'
    else:
        raise ValueError(variant)
    output.mkdir(parents=True, exist_ok=False)
    inp = output/'input'; inp.mkdir()
    filename = case_name+'.i'; (inp/filename).write_text(text)
    for entry in old['mesh_files']+old['extra_files']:
        if entry['name'] == 'diagnostic-manifest.json':
            continue
        source = parent/'input'/entry['name']
        assert sha(source) == entry['sha256']
        (inp/entry['name']).write_bytes(source.read_bytes())
    record = {'test_variant': variant, 'parent_directory': str(parent),
        'parent_input_sha256': old['input_sha256'], 'single_change': change,
        'builder_sha256': sha(Path(__file__)),
        'classification': 'prototype diagnostic; not engineering acceptance',
        'preserved': 'four tables, E/nu, wt/wc, viscosity, loading, local integrator and solver release',
        'limitations': 'IP samples are accepted global states only; no failed local states or complete wall-time counters'}
    (inp/'diagnostic-manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    if extra_evidence:
        (inp/extra_evidence_name).write_text(json.dumps(extra_evidence,ensure_ascii=False,indent=2)+'\n')
    manifest = dict(old, case_name=case_name, input_file=filename,
        input_sha256=sha(inp/filename),
        command='DamSafetyApp-opt -i '+filename,
        extra_files=[{'name':p.name,'sha256':sha(p)} for p in sorted(inp.iterdir())
                     if p.name != filename and p.name not in {f['name'] for f in old['mesh_files']}])
    (output/'submission.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--variant',choices=['observe','tight-abs','no-scaling','top-clamp','early-window','dt001','dt0001','replay-path','early-refined-grid','full-short-observe','dual-ip-path-trace'],required=True)
    p.add_argument('--case-name',required=True)
    args = p.parse_args()
    result = build(args.parent,args.output,args.variant,args.case_name)
    print(json.dumps({k:result[k] for k in ('case_name','input_file','input_sha256')}))
