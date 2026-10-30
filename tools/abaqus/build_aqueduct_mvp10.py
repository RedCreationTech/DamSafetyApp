#!/usr/bin/env python3
"""Compose the audited AQ1030 full transient prototype; never submit or solve it.

Source stabilization and the separate prestressed eigenanalysis remain explicit
expert gates. The generated INITIAL Terminator prevents accidental FE execution.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path
import abaqus2exodus

SOURCE_SHA='9f71f0f898755431f3ba450bf1aaaa354a08ec8a6b09e5e62ab2d41b1d65e4a1'
BASELINE_SHA='ce41cfbdae51bcf70cfe2a6c037f4d8d4f2935c15ec02a66864911698a1583a6'
ALPHA=-.41421
STATIC_END=2.
END=12.

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def quote(items):
    return "'"+' '.join(map(str,items))+"'"

def parse_blocks(text):
    """Strict single-line HIT subset used by the immutable AQ1023 input."""
    blocks=collections.OrderedDict(); stack=[]
    for n,line in enumerate(text.splitlines(),1):
        s=line.strip()
        if not s or s.startswith('#'):continue
        if s.startswith('[') and s.endswith(']'):
            name=s[1:-1]
            if name in ('','../'):
                if not stack:raise ValueError('unmatched block end at '+str(n))
                stack.pop()
            else:
                if name.startswith('./'):name=name[2:]
                path='/'.join([*stack,name])
                if path in blocks:raise ValueError('repeated block '+path)
                blocks[path]=collections.OrderedDict(); stack.append(name)
        else:
            if not stack or '=' not in s:raise ValueError('unsupported HIT line '+str(n))
            key,value=map(str.strip,s.split('=',1)); path='/'.join(stack)
            if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*',key) or key in blocks[path]:
                raise ValueError('invalid or duplicate parameter '+key)
            blocks[path][key]=value
    if stack:raise ValueError('unclosed input block')
    return blocks

def render(blocks):
    return ''.join('['+name+']\n'+''.join('  '+k+' = '+v+'\n' for k,v in params.items())+'[]\n'
                   for name,params in blocks.items())

def add(blocks,path,**params):
    if path in blocks:raise ValueError('would overwrite block '+path)
    blocks[path]=collections.OrderedDict((k,str(v)) for k,v in params.items())

def drop_tree(blocks,prefix):
    for k in list(blocks):
        if k==prefix or k.startswith(prefix+'/'):del blocks[k]

def sampling_mesh(source,nodes,elements,critical,path):
    """A compact output mesh only; physical mesh/IDs/interfaces are never changed."""
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/'source.inp';p.write_text(source.read_bytes().decode('gbk'),encoding='utf-8')
        model=abaqus2exodus.parse_inp(str(p))
    parts={i.name:i.part for i in model.instances}
    node_map={(r['instance'],int(r['abaqus_node_label'])):r for r in nodes}
    used=set(); conn=[]
    for r in critical:
        part=model.parts[parts[r['instance']]]; eid=int(r['abaqus_element'])
        if part.elem_types[eid]!='C3D8R':raise ValueError('critical output expects HEX8')
        keys=[(r['instance'],n) for n in part.elems[eid]]
        used.update(keys); conn.append((r,keys))
    keys=sorted(used);local={k:i+1 for i,k in enumerate(keys)}
    lines=['$MeshFormat','2.2 0 8','$EndMeshFormat','$Nodes',str(len(keys))]
    for k in keys:
        r=node_map[k];lines.append(str(local[k])+' '+' '.join(r[c] for c in ['x_m','y_m','z_m']))
    lines+=['$EndNodes','$Elements',str(len(conn))]
    for j,(r,kk) in enumerate(conn,1):
        lines.append(f'{j} 5 2 1 1 '+' '.join(str(local[k]) for k in kk))
    lines+=['$EndElements']
    path.write_text('\n'.join(lines)+'\n')
    return {'nodes':len(keys),'elements':len(conn),'output_only':True,
            'meaning':'Sampled element-mean cloud, not direct original element IDs; native critical CSV is authoritative'}

def compose(b,controls,critical,gravity):
    """All unchanged blocks keep the frozen physical parameter values."""
    for name,p in b.items():
        if 'hht_start_time' in p:p['hht_start_time']='2'
        if p.get('type')=='IntegratedAcceleration':p['start_time']='2'
        if p.get('type') in ('NewmarkAccelAux','NewmarkVelAux'):
            p['type']='AqueductNewmarkAcceleration' if p['type']=='NewmarkAccelAux' else 'AqueductNewmarkVelocity'
            p['static_end_time']='2'
    b['Executioner'].update(end_time='12',dtmin='1e-10',dtmax='.02')
    b['Executioner/TimeIntegrator'].update(quasistatic_initialization='true',inactive_tsteps='0')
    b['Executioner/TimeStepper'].update(type='AqueductPhaseDT',static_end_time='2',
        static_dtmin='1e-10',static_dtmax='.02',dynamic_dtmin='1e-8',dynamic_dtmax='.01',growth_factor='2')
    b['Functions/dynamic_on']['expression']="'t>2'"
    b['Functions/dt_clock']=collections.OrderedDict(type='PiecewiseConstant',direction='LEFT_INCLUSIVE',
        x="'0 .001 1 1.001 2 2.001 12'",y="'.001 .02 .001 .02 .001 .0025 .0025'")
    for name,p in parse_blocks(gravity).items():
        if name in b:raise ValueError('gravity block collision')
        b[name]=p
    concrete=b['Physics/SolidMechanics/QuasiStatic/concrete']['block']
    # Read-only material diagnostics, not newly introduced constitutive updates.
    plastic=[]
    for i,j,axis in [(0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')]:
        name='cdp_pe_'+axis;plastic.append(name)
        add(b,'AuxVariables/'+name,block=concrete,family='MONOMIAL',order='CONSTANT')
        add(b,'AuxKernels/'+name,type='RankTwoAux',variable=name,rank_two_tensor='cdp_viscous_plastic_strain',
            index_i=i,index_j=j,block=concrete,execute_on='TIMESTEP_END')
    for name,prop in [('cdp_kappa_t','cdp_kappa_t'),('cdp_kappa_c','cdp_kappa_c'),
                      ('cdp_combined_damage','cdp_combined_damage')]:
        add(b,'AuxVariables/'+name,block=concrete,family='MONOMIAL',order='CONSTANT')
        add(b,'AuxKernels/'+name,type='MaterialRealAux',variable=name,property=prop,block=concrete,execute_on='TIMESTEP_END')
        plastic.append(name)
    add(b,'Postprocessors/aq_total_time',type='TimePostprocessor',execute_on="'INITIAL TIMESTEP_END'")
    pp=[]
    for row in controls:
        for field in ['ux','uy','uz','vx','vy','vz','ax','ay','az']:
            name='aq_'+row['name']+'_'+field;pp.append(name)
            add(b,'Postprocessors/'+name,type='NodalVariableValue',variable=field,nodeid=row['libmesh_node'],
                execute_on="'INITIAL TIMESTEP_END'")
    # Existing complete residual tags include gravity and constraints. This is a
    # derived HHT quantity, deliberately not named direct Abaqus RF.
    reaction=[]
    for a in 'xyz':
        name='aq_derived_support_'+a;reaction.append(name)
        add(b,'Postprocessors/'+name,type='ParsedPostprocessor',
            pp_names=quote(['aq_total_time','base_time_f'+a,'base_nontime_f'+a,'base_old_nontime_f'+a]),
            expression=quote([f'base_time_f{a}+base_nontime_f{a}-if(aq_total_time>2,{ALPHA/(1+ALPHA):.17g}*base_old_nontime_f{a},0)']),
            execute_on="'INITIAL TIMESTEP_END'")
    stress=[kind+'_'+s for kind in ['stress','strain'] for s in ['xx','yy','zz','xy','xz','yz']]
    damage=['damaget','damagec','cdp_combined_damage']
    region=stress+plastic+['damaget','damagec']
    ids=[int(r['exodus_element'])-1 for r in critical]
    if len(ids)!=len(set(ids)) or not ids:raise ValueError('critical IDs not unique')
    add(b,'VectorPostprocessors/aq_critical_native',type='AqueductElementSampler',element_ids=quote(ids),
        variable=quote(region),block=concrete,sort_by='id',execute_on='TIMESTEP_END')
    connectors=[k.split('/')[-1] for k in b if k.startswith('Postprocessors/connector_')]
    contact=[k.split('/')[-1] for k in b if k.startswith('Postprocessors/') and any(
        k.split('/')[-1].startswith(v) for v in ['normal_lm','tangent_1_','tangent_2_'])]
    overview=['ux','uy','uz']+stress+[k.split('/')[-1] for k in b if
        k.startswith('AuxVariables/') and (k.split('/')[-1].startswith('bar_') or k.split('/')[-1].startswith('shell_q'))]
    drop_tree(b,'Outputs')
    add(b,'Outputs',csv='false',exodus='false',console='false',print_linear_residuals='false')
    add(b,'Outputs/terminal',type='Console',max_rows='2')
    add(b,'Outputs/profile',type='PerfGraphOutput',level='5',heaviest_sections='25',heaviest_branch='true')
    add(b,'Outputs/gravity_fields',type='Exodus',file_base='aq1030_gravity',end_time='2',
        show=quote(overview+damage+plastic),execute_on="'INITIAL TIMESTEP_END'")
    add(b,'Outputs/gravity_history',type='CSV',file_base='aq1030_gravity_history',end_time='2',
        show=quote(['aq_total_time']+pp+reaction+connectors+contact),execute_on="'INITIAL TIMESTEP_END'")
    schedule=dict(phase_origin='2',phase_end='12',execute_on="'TIMESTEP_END FINAL'")
    add(b,'Outputs/overview',type='AqueductScheduledExodus',file_base='aq1030_overview',mark_interval='.1',
        show=quote(overview),**schedule)
    add(b,'Outputs/concrete_damage',type='AqueductScheduledExodus',file_base='aq1030_concrete_damage',
        mark_interval='.5',show=quote(damage),sampling_blocks=concrete,**schedule)
    add(b,'Outputs/critical_cloud',type='AqueductScheduledExodus',file_base='aq1030_critical',
        file='critical-output.msh',mark_interval='.05',show=quote(region),**schedule)
    add(b,'Outputs/critical_native',type='AqueductScheduledCSV',file_base='aq1030_critical_native',
        mark_interval='.05',show="'aq_critical_native'",time_data='true',**schedule)
    add(b,'Outputs/controls',type='AqueductScheduledCSV',file_base='aq1030_controls',
        mark_interval='.01',show=quote(['aq_total_time']+pp+connectors+contact),**schedule)
    add(b,'Outputs/support_reaction',type='AqueductScheduledCSV',file_base='aq1030_support',
        mark_interval='.1',show=quote(['aq_total_time']+reaction),**schedule)
    audit=[k.split('/')[-1] for k in b if k.startswith('Postprocessors/') and any(
        k.split('/')[-1].startswith(s) for s in ['global_','base_time_','base_nontime_','base_old_nontime_','spin_','cross_','max_','mean_'])]
    add(b,'Outputs/accepted_audit',type='CSV',file_base='aq1030_accepted_audit',
        show=quote(['aq_total_time']+audit),execute_on="'INITIAL TIMESTEP_END'")
    add(b,'UserObjects/aq_expert_pending',type='Terminator',expression="'1'",execute_on='INITIAL',
        fail_mode='HARD',error_level='ERROR',message="'AQ1030 prototype: expert stabilization/modal decisions and candidate binary precheck are pending. No FE run authorized.'")
    return dict(control_channels=len(pp),critical_elements=len(ids),critical_variables=region,
        overview_variables=overview,native_critical_sampler='aq_critical_native',
        derived_support_reaction='TIME+NONTIME-alpha/(1+alpha)*oldNONTIME in dynamics; not direct Abaqus RF')

def build(audited_input,baseline_input,source_inp,output):
    audited,baseline,source,out=map(Path,[audited_input,baseline_input,source_inp,output])
    if out.exists():raise ValueError('output exists; preserve previous prototype')
    if digest(source)!=SOURCE_SHA or digest(baseline)!=BASELINE_SHA:raise ValueError('unexpected source or baseline identity')
    audit=json.loads((audited/'source-audit.json').read_text())
    if audit['source']['sha256']!=SOURCE_SHA or not all(audit['geometry'].values()):raise ValueError('source audit is not applicable')
    for name,h in audit['artifacts'].items():
        if digest(audited/name)!=h:raise ValueError('audited asset hash mismatch: '+name)
    nodes=list(csv.DictReader((audited/'node-map.csv').open()))
    elements=list(csv.DictReader((audited/'element-map.csv').open()))
    controls=list(csv.DictReader((audited/'control-node-map.csv').open()))
    critical=list(csv.DictReader((audited/'critical-element-map.csv').open()))
    if len(nodes)!=36158 or len(elements)!=38534 or len(controls)!=40 or len(critical)!=798:
        raise ValueError('incomplete full assembly or output region')
    b=parse_blocks(baseline.read_text());stats=compose(b,controls,critical,(audited/'aq1030_mvp10_gravity_include.i').read_text())
    out.mkdir(parents=True)
    for p in audited.iterdir():
        if p.is_file() and p.name not in ['source-audit.json','aq1030_mvp10_output_include.i','aq1030_mvp10_gravity_include.i']:
            shutil.copyfile(p,out/p.name)
    sample=sampling_mesh(source,nodes,elements,critical,out/'critical-output.msh')
    main=out/'aq1030_mvp10_full.i'
    main.write_text('# AQ1030 full assembly: GRAV_EQ0..1 + GRAV_CHECK1..2 + ELC2..12.\n'
                    '# LOCAL PROTOTYPE ONLY; INITIAL gate enabled; no Abaqus equivalence or node precheck.\n'
                    '# FREQ_PRE20 modes are a separate unresolved branch, not a fictitious transient step.\n'+render(b))
    code=Path(__file__).resolve().parents[2]
    sources={str(p.relative_to(code)):digest(p) for directory in ['include','src']
             for p in (code/directory).rglob('Aqueduct*') if p.is_file()}
    sources[str(Path(__file__).resolve().relative_to(code))]=digest(__file__)
    identity=dict(schema='aq1030-main-input/1',status='local-prototype-not-prechecked',
        full_transient_main_generated=True,ready_for_submission=False,job_id=None,
        source_sha256=SOURCE_SHA,frozen_physical_input_sha256=BASELINE_SHA,
        preferred_future_node='192.168.1.35; re-read LIMS config before any node operation',
        timeline={'GRAV_EQ':[0,1],'GRAV_CHECK':[1,2],'ELC':[2,12],'modal':'separate20mode branch pending'},
        physics='all source assembly/material/interface objects retained; gravity added; no constitutive cache enabled',
        motion='unchanged whole-record corrected tables; shift origin1->2, use first10s, no refit',
        time_policy='FunctionDT engineering schedule with phase knots; not Abaqus automatic accuracy-controlled stepping',
        source_time_bounds={'static_initial':.001,'static_min':1e-10,'static_max':.02,'dynamic_initial':.001,'dynamic_min':1e-8,'dynamic_max':.01},
        candidate_dynamic_nominal_dt=.0025,nonlinear_tolerances='frozen 1e-5 absolute /1e-8 relative /35iterations/basic line search',
        state_handoff='continuous original solution/material/contact and accepted-force history; both Newmark vector/AD derivatives and six rate auxiliaries zero through total2s',
        output=stats,critical_cloud=sample,
        expert_pending=['A01 nonzero-gravity stabilization .0002 not implemented/approval pending',
                        'A02 separate prestressed20mode branch sequencing/implementation pending',
                        'A03 matching new reference and acceptance criteria'],
        remaining_output_gaps=['Abaqus direct RF field/history; derived HHT support force is separately named',
            'Formulation-specific PE/PEEQ and SDEG exact semantics; CDP native properties keep distinct names',
            'Abaqus10term energies; connector CU/CTF and contact CDISP/CSTRESS full native mapping',
            'Critical sampled cloud must be paired against direct ID-selected CSV in first isolated validation'],
        initial_gate='UserObjects/aq_expert_pending enabled; stops before FE until review/authorization',
        source_file_hashes=sources,source_cpp_compiled=False,moose_check_input=False,seven_cdp='deferred/not run/not passed',
        artifacts={p.name:digest(p) for p in out.iterdir() if p.is_file()})
    (out/'main-input-identity.json').write_text(json.dumps(identity,ensure_ascii=False,indent=2)+'\n')
    return identity

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['audited-input','baseline-input','source-inp','output']:p.add_argument('--'+n,required=True)
    a=p.parse_args();result=build(a.audited_input,a.baseline_input,a.source_inp,a.output)
    print(json.dumps({k:result[k] for k in ['status','full_transient_main_generated','ready_for_submission','job_id','output']},ensure_ascii=False,indent=2))
