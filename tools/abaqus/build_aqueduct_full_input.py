#!/usr/bin/env python3
"""Compose all AQ1023 source objects into a short-window full-assembly candidate.

Custom HEX8+Bbar, MITC4, mortar Tie, contact and HHT deviations are explicit in
full-input-identity.json. This is not Abaqus equivalence or final acceptance.
"""
import argparse,csv,json,shutil
from pathlib import Path
from collections import defaultdict
import numpy as np
import netCDF4
import abaqus2exodus as mesh
from map_aqueduct_interfaces import read_geometry,sha
from build_aqueduct_interface_cases import block

FIELDS=('ux','uy','uz','rx','ry','rz')
ALPHA=-.41421;BETA=(1-ALPHA)**2/4;GAMMA=.5-ALPHA

def quoted(values):return "'"+' '.join(map(str,values))+"'"

def build(preparation,motion,output,end_time=1.01,mumps_memory_retry=False,diagnostics=False,ignore_zero_entries=False,hash_assembly=False):
    root=Path(preparation);out=Path(output);motion=Path(motion)
    if out.exists():raise ValueError('output exists')
    c=json.loads((root/'assembly-contract.json').read_text());implementation=json.loads((root/'interface-implementation.json').read_text())
    for n,h in c['artifacts'].items():
        if sha(root/n)!=h:raise ValueError('preparation identity mismatch '+n)
    xyz,conn,bids,types,origins=read_geometry(root);model=mesh.parse_inp(root/'source-utf8.inp');parts={i.name:i.part for i in model.instances}
    sections={(s['part'],e):s for s in c['sections'] for e in s['element_labels']};group={}
    for origin,eid in origins.items():
        section=sections[(parts[origin[0]],origin[1])];bid=bids[eid]
        signature=(types[eid],section['material'],section.get('area_m2'),section.get('thickness_m'))
        if bid in group and group[bid]['signature']!=signature:raise ValueError('mixed physical sections in one block')
        group.setdefault(bid,dict(signature=signature,section=section,ids=[]))['ids'].append(eid)
    solids=[b for b,g in group.items() if g['signature'][0]=='HEX8'];shell=[b for b,g in group.items() if g['signature'][0]=='SHELL4'];bars=[b for b,g in group.items() if g['signature'][0]=='TRUSS'];concrete=[b for b in solids if group[b]['signature'][1]=='con'];steel=[b for b in solids if b not in concrete]
    if sum(len(group[b]['ids']) for b in solids)!=15476 or sum(len(group[b]['ids']) for b in bars)!=22770 or sum(len(group[b]['ids']) for b in shell)!=288:raise ValueError('source element coverage differs')
    out.mkdir(parents=True)
    for name in ('assembly-geometry.e','embedded-runtime-map.csv','node-map.csv','element-map.csv','interface-implementation.json'):shutil.copyfile(root/name,out/name)
    for a in 'xyz':shutil.copyfile(motion/f'motion-{a}.csv',out/f'motion-{a}.csv')
    shutil.copyfile(motion/'motion-identity.json',out/'motion-identity.json')
    rows=implementation['rows'];support=defaultdict(set)
    for r in rows:support[r['multiplier']].add(r['anchor'])
    meshblocks=block('source','type = FileMeshGenerator\nfile = assembly-geometry.e');previous='source';lmnames=[]
    for j,nodes in sorted(support.items()):
        name=f'lm_support_{j}';meshblocks+=block(name,f'type = AqueductConstraintNodesGenerator\ninput = {previous}\nnodes = {quoted(sorted(nodes))}\nblock_id = {200+j}\nblock_name = {name}');previous=name;lmnames.append(f'lm{j}')
    embedded_nodes=[]
    with (root/'embedded-runtime-map.csv').open() as f:embedded_nodes=[int(r['secondary_node']) for r in csv.DictReader(f)]
    meshblocks+=block('embedded_support',f'type = AqueductConstraintNodesGenerator\ninput = {previous}\nnodes = {quoted(embedded_nodes)}\nblock_id = 220\nblock_name = embedded_support');previous='embedded_support'
    lmnames+=['embedded_lm_x','embedded_lm_y','embedded_lm_z']
    for j in range(1,5):
        for k,role in enumerate(('secondary','primary')):
            name=f'contact_{j}_{role}';meshblocks+=block(name,f'type = LowerDBlockFromSidesetGenerator\ninput = {previous}\nsidesets = aq_contact_{j:02d}_{role}\nnew_block_id = {300+2*j+k}\nnew_block_name = {name}');previous=name
    text='# Full AQ1023 assembly; first short candidate, engineering surrogates recorded.\n'+block('Mesh',meshblocks+'allow_renumbering = false\nparallel_type = replicated\npatch_update_strategy = iteration')
    text+=block('Problem','material_coverage_check = false'+('\nignore_zeros_in_jacobian = true' if ignore_zero_entries else '')+('\nuse_hash_table_matrix_assembly = true' if hash_assembly else ''))
    text+=block('GlobalParams',"displacements = 'ux uy uz'")
    variables=''.join(block(v,'block = '+quoted(sorted(group)) if j<3 else 'block = '+quoted(shell)) for j,v in enumerate(FIELDS))
    variables+=''.join(block(f'lm{j}',f'block = {200+j}') for j in sorted(support))
    for j in range(1,5):variables+=''.join(block(v+str(j),f'block = contact_{j}_secondary\nuse_dual = true\nscaling = .01') for v in ('normal_lm','tangent_1_','tangent_2_'))
    variables+=''.join(block('embedded_lm_'+a,'block = 220') for a in 'xyz')
    text+=block('Variables',variables)
    physics=''
    for name,bb in [('concrete',concrete),('steel',steel)]:physics+=block(name,'block = '+quoted(bb)+"\nstrain = SMALL\nincremental = true\nvolumetric_locking_correction = true\nsave_in = 'fx fy fz'\ngenerate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'")
    text+=block('Physics/SolidMechanics/QuasiStatic',physics)
    aux=''.join(block(v,'block = '+quoted(sorted(group))) for v in ('fx','fy','fz','vx','vy','vz','ax','ay','az'))
    for b in bars:aux+=block(f'area_{b}',f'block = {b}\nfamily = MONOMIAL\norder = CONSTANT\ninitial_condition = {group[b]["section"]["area_m2"]}')
    aux+=''.join(block(v,'block = '+quoted(concrete)+'\nfamily = MONOMIAL\norder = CONSTANT') for v in ('damaget','damagec'))
    if diagnostics:
        for v in ('bar_stress','bar_plastic','bar_peq','bar_total_strain','bar_elastic_strain'):
            aux+=block(v,'block = '+quoted(bars)+'\nfamily = MONOMIAL\norder = CONSTANT')
        for qp in range(4):
            for layer in range(3):
                for v in ('sxx','syy','szz','sxy','peq'):
                    aux+=block(f'shell_q{qp}_l{layer}_{v}','block = '+quoted(shell)+'\nfamily = MONOMIAL\norder = CONSTANT')
        for tag in ('time','nontime'):
            for j,v in enumerate(FIELDS):aux+=block(f'{tag}_{v}','block = '+quoted(sorted(group) if j<3 else shell))
    text+=block('AuxVariables',aux)
    materials='';functions=block('steel_hardening',"type = PiecewiseLinear\nx = '0 .1'\ny = '400e6 696e6'\nextrap = false")
    for name,bb,mat in [('concrete',concrete,'con'),('steel',steel,'steel')]:
        E,nu=c['materials_si'][mat]['elastic'][0];rho=c['materials_si'][mat]['density'][0][0];bp='block = '+quoted(bb)
        materials+=block(name+'_elastic',bp+f'\ntype = ComputeIsotropicElasticityTensor\nyoungs_modulus = {E}\npoissons_ratio = {nu}')
        materials+=block(name+'_density',bp+f'\ntype = GenericConstantMaterial\nprop_names = density\nprop_values = {rho}')
        if mat=='con':
            params=c['materials_si'][mat]['cdp'][0]
            materials+=block('concrete_stress',bp+'\ntype = ComputeMultipleInelasticStress\ninelastic_models = cdp_update\nperform_finite_strain_rotations = false')
            materials+=block('cdp_update',bp+f'\ntype = AbaqusCDPStressUpdate\nyoungs_modulus = {E}\npoissons_ratio = {nu}\ndilation_angle = {params[0]}\neccentricity = {params[1]}\nbiaxial_to_uniaxial_compression_ratio = {params[2]}\ntensile_meridian_ratio = {params[3]}\nviscosity = {params[4]}\ntension_recovery = 0\ncompression_recovery = 1\nmaximum_substeps = 256\nmaximum_strain_increment = 2.5e-5\nenable_performance_diagnostics = true\n'+ '\n'.join(f'{n}_file = {n}.csv' for n in ('compression_hardening','compression_damage','tension_stiffening','tension_damage')))
            for n in ('compression_hardening','compression_damage','tension_stiffening','tension_damage'):
                with (out/f'{n}.csv').open('w',newline='') as f:
                    w=csv.writer(f);w.writerow([('damage_c' if n.startswith('compression') else 'damage_t') if n.endswith('damage') else 'stress_pa','inelastic_strain' if n.startswith('compression') else 'cracking_strain']);w.writerows(c['materials_si']['con']['cdp_'+n])
        else:
            materials+=block('steel_stress',bp+'\ntype = ComputeMultipleInelasticStress\ninelastic_models = steel_update\ntangent_operator = nonlinear')+block('steel_update',bp+'\ntype = IsotropicPlasticityStressUpdate\nyield_stress = 400e6\nhardening_function = steel_hardening')
    materials+=block('truss',f'type = AqueductSteelTruss\nblock = {quoted(bars)}\ndisplacements = \'ux uy uz\'\nyoungs_modulus = 200e9\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1')
    materials+=block('shell',f'type = AqueductFlatShellStress\nblock = {quoted(shell)}\ndisplacements = \'ux uy uz\'\nrotations = \'rx ry rz\'\nyoungs_modulus = 200e9\npoissons_ratio = .3\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1\nthickness = .1\noffset = .5\nshear_correction = .8333333333333333\ndrilling_factor = .001')
    for b in bars:materials+=block(f'line_mass_{b}',f'type = GenericConstantMaterial\nblock = {b}\nprop_names = density\nprop_values = {7850*group[b]["section"]["area_m2"]}')
    text+=block('Materials',materials)
    kernels=''.join(block(f'lm{j}',f'type = NullKernel\nvariable = lm{j}\nblock = {200+j}\ndisplacements = \'\'\njacobian_fill = 0') for j in sorted(support))
    kernels+=''.join(block('embedded_lm_'+a,f'type = NullKernel\nvariable = embedded_lm_{a}\nblock = 220\ndisplacements = \'\'\njacobian_fill = 0') for a in 'xyz')
    for b in bars:
        for j,v in enumerate(FIELDS[:3]):kernels+=block(f'bar_{b}_{v}',f'type = StressDivergenceTensorsTruss\nvariable = {v}\nblock = {b}\ncomponent = {j}\narea = area_{b}\nuse_displaced_mesh = false')
    for j,v in enumerate(FIELDS):
        kernels+=block('shell_'+v,f'type = AqueductFlatShell\nvariable = {v}\ncomponent = {j}\nblock = {quoted(shell)}\nrotations = \'rx ry rz\'')
        kernels+=block('shell_mass_'+v,f'type = AqueductFlatShellInertia\nvariable = {v}\ncomponent = {j}\nblock = {quoted(shell)}\nrotations = \'rx ry rz\'\ndensity = 7850\nthickness = .1\noffset = .5\nhht_alpha = {ALPHA}\nhht_start_time = 1')
    for v,axis in zip(FIELDS,'xyz'):kernels+=block('mass_'+v,f'type = AqueductHHTInertia\nvariable = {v}\nblock = {quoted(solids+bars)}\nbeta = {BETA}\ngamma = {GAMMA}\nvelocity = v{axis}\nacceleration = a{axis}\nhht_alpha = {ALPHA}\nhht_start_time = 1\nuse_displaced_mesh = false\nsave_in = f{axis}')
    text+=block('Kernels',kernels)
    auxkernels=''.join(block(v,f'type = MaterialRealAux\nvariable = {v}\nproperty = Damage{v[-1].upper()}\nblock = {quoted(concrete)}\nexecute_on = TIMESTEP_END') for v in ('damaget','damagec'))
    for v,axis in zip(FIELDS,'xyz'):
        auxkernels+=block('accel_'+axis,f'type = NewmarkAccelAux\nvariable = a{axis}\ndisplacement = {v}\nvelocity = v{axis}\nbeta = {BETA}\nexecute_on = timestep_end')+block('velocity_'+axis,f'type = NewmarkVelAux\nvariable = v{axis}\nacceleration = a{axis}\ngamma = {GAMMA}\nexecute_on = timestep_end')
    if diagnostics:
        for v,prop in [('bar_stress','axial_stress'),('bar_plastic','signed_plastic_stretch'),('bar_peq','equivalent_plastic_stretch'),('bar_total_strain','total_stretch'),('bar_elastic_strain','elastic_stretch')]:
            auxkernels+=block(v,f'type = MaterialRealAux\nvariable = {v}\nproperty = {prop}\nblock = {quoted(bars)}\nexecute_on = TIMESTEP_END')
        for qp in range(4):
            for layer in range(3):
                for j,v in enumerate(('sxx','syy','szz','sxy','peq')):
                    name=f'shell_q{qp}_l{layer}_{v}'
                    auxkernels+=block(name,f'type = MaterialStdVectorAux\nvariable = {name}\nproperty = aq_shell_state\nindex = {5*layer+j}\nselected_qp = {qp}\nblock = {quoted(shell)}\nexecute_on = TIMESTEP_END')
        for tag in ('time','nontime'):
            for j,v in enumerate(FIELDS):auxkernels+=block(f'{tag}_{v}',f'type = TagVectorAux\nvariable = {tag}_{v}\nv = {v}\nvector_tag = {tag}\nremove_variable_scaling = true\nblock = {quoted(sorted(group) if j<3 else shell)}')
    text+=block('AuxKernels',auxkernels)
    constraints=''
    for i,r in enumerate(rows):constraints+=block(f'mpc_{i}',f'type = AqueductLinearMPC\nvariable = ux\nanchor_node = {r["anchor"]}\nterm_nodes = {quoted(t["node"] for t in r["terms"])}\nterms = {quoted(t["variable"] for t in r["terms"])}\ncoefficients = {quoted(t["coefficient"] for t in r["terms"])}\nmultiplier = lm{r["multiplier"]}')
    allowed=sorted({b for g in implementation['embedded_groups'] for b in g['primary_blocks']})
    constraints+=block('embedded_mpc',f'type = AqueductEmbeddedMPC\nvariable = ux\ndisplacements = \'ux uy uz\'\nmultipliers = \'embedded_lm_x embedded_lm_y embedded_lm_z\'\nallowed_primary_block_ids = {quoted(allowed)}\nsource_map = embedded-runtime-map.csv')
    for j,r in enumerate(c['interfaces']['connectors']):
        a,b=[x['node']-1 for x in r['endpoints']]
        for v,k in zip(FIELDS,[605000,605000,21965000]):constraints+=block(f'connector_{j}_{v}',f'type = LinearNodalConstraint\nvariable = {v}\nprimary = {a}\nsecondary_node_ids = {b}\nweights = 1\npenalty = {k}\nformulation = penalty')
    userobjects=''
    for j in range(1,5):
        common=f'primary_boundary = aq_contact_{j:02d}_primary\nsecondary_boundary = aq_contact_{j:02d}_secondary\nprimary_subdomain = contact_{j}_primary\nsecondary_subdomain = contact_{j}_secondary\ncorrect_edge_dropping = true\nuse_displaced_mesh = true'
        userobjects+=block(f'slip_{j}',common+f'\ntype = AqueductElasticSlipUserObject\nlm_variable_normal = normal_lm{j}\nlm_variable_tangential_one = tangent_1_{j}\nlm_variable_tangential_two = tangent_2_{j}\nsecondary_variable = ux\ndisp_x = ux\ndisp_y = uy\ndisp_z = uz')
        constraints+=block(f'contact_law_{j}',common+f'\ntype = AqueductElasticSlipConstraint\nvariable = normal_lm{j}\ndisp_x = ux\ndisp_y = uy\ndisp_z = uz\nfriction_lm = tangent_1_{j}\nfriction_lm_dir = tangent_2_{j}\nweighted_gap_uo = slip_{j}\nweighted_velocities_uo = slip_{j}\ncompute_primal_residuals = false\nnormalize_c = true\nc = 1e6\nmu = .55\ncritical_slip = .00025331')
        for axis,v in zip('xyz',FIELDS):
            for nm,typ,lm,uo,extra in [('n','NormalMortarMechanicalContact',f'normal_lm{j}','weighted_gap_uo',''),('t1','TangentialMortarMechanicalContact',f'tangent_1_{j}','weighted_velocities_uo',''),('t2','TangentialMortarMechanicalContact',f'tangent_2_{j}','weighted_velocities_uo','\ndirection = direction_2')]:constraints+=block(f'contact_{j}_{nm}_{v}',common+f'\ntype = {typ}\nvariable = {lm}\nsecondary_variable = {v}\ncomponent = {axis}\ncompute_lm_residuals = false\n{uo} = slip_{j}'+extra)
        lmnames += [f'normal_lm{j}',f'tangent_1_{j}',f'tangent_2_{j}']
    text+=block('UserObjects',userobjects)+block('Constraints',constraints)
    functions+=block('dynamic_on',"type = ParsedFunction\nexpression = 't>1'")+block('dt_clock',"type = ParsedFunction\nexpression = 'if(t<1,1,.0025)'")
    bcs=''
    # Explicit source-instance nodesets (avoid truncated name guesses).
    with netCDF4.Dataset(out/'assembly-geometry.e') as nc:
        ns_names=[s for s in netCDF4.chartostring(nc['ns_names'][:]) if s.startswith('Set_bottom__')]
    if len(ns_names)!=2:raise ValueError('base boundary coverage differs')
    for v,axis in zip(FIELDS,'xyz'):
        functions+=block('motion_a'+axis,f'type = IntegratedAcceleration\ndata_file = motion-{axis}.csv\ncomponent = acceleration\nstart_time = 1')
        bcs+=block('fixed_'+axis,f'type = DirichletBC\nvariable = {v}\nboundary = {quoted(ns_names)}\nvalue = 0')+block('motion_'+axis,f'type = PresetAcceleration\nvariable = {v}\nboundary = {quoted(ns_names)}\nfunction = motion_a{axis}\nvelocity = v{axis}\nacceleration = a{axis}\nbeta = {BETA}')
    text+=block('Functions',functions)+block('BCs',bcs)
    text+=block('Controls',block('phase',f'type = ConditionalFunctionEnableControl\nconditional_function = dynamic_on\nenable_objects = {quoted("BCs::motion_"+a for a in "xyz")}\ndisable_objects = {quoted("BCs::fixed_"+a for a in "xyz")}\nexecute_on = \'initial timestep_begin\''))
    post=''
    for v in ('damaget','damagec'):
        post+=block('max_'+v,f'type = ElementExtremeValue\nvariable = {v}\nblock = {quoted(concrete)}\nvalue_type = max')+block('mean_'+v,f'type = ElementAverageValue\nvariable = {v}\nblock = {quoted(concrete)}')
    probes=[]
    # Real source nodes, one high node per instance group plus connector endpoints.
    with (root/'node-map.csv').open() as f:nmap=list(csv.DictReader(f))
    for inst in ('caoshen（dangkuai-1','paijia-1','paijia-1-lin-2-1','paijia-DI-1','paijia-DI-1-lin-2-1','gangliang-1'):
        candidates=[r for r in nmap if r['instance'].lower()==inst.lower()]
        if not candidates:raise ValueError('source probe instance missing: '+inst)
        r=max(candidates,key=lambda r:float(r['z_m']));probes.append(dict(name='probe_'+str(len(probes)),node=int(r['exodus_node'])-1,instance=r['instance'],abaqus_local_node=int(r['abaqus_node_label']),position_m=[float(r[k]) for k in ('x_m','y_m','z_m')]))
    for r in probes:
        for v in (*FIELDS[:3],'vx','vy','vz','ax','ay','az'):post+=block(r['name']+'_'+v,f'type = NodalVariableValue\nvariable = {v}\nnodeid = {r["node"]}')
    for j,r in enumerate(c['interfaces']['connectors']):
        for endp in r['endpoints']:
            n=endp['node']-1
            for v in FIELDS:post+=block(f'connector_{j}_n{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for j in range(1,5):
        for v in (f'normal_lm{j}',f'tangent_1_{j}',f'tangent_2_{j}'):
            for extreme in ('min','max'):
                post+=block(f'{v}_{extreme}',f'type = NodalExtremeValue\nvariable = {v}\nboundary = aq_contact_{j:02d}_secondary\nvalue_type = {extreme}')
    for axis in 'xyz':post+=block('base_current_force_'+axis,f'type = NodalSum\nvariable = f{axis}\nboundary = {quoted(ns_names)}')
    text+=block('Postprocessors',post)+block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner',f'type = Transient\nsolve_type = NEWTON\nend_time = {end_time}\ndtmin = 1e-8\nnl_abs_tol = 1e-5\nnl_rel_tol = 1e-8\nnl_max_its = 35\nline_search = basic\npetsc_options_iname = \'-pc_type -pc_factor_mat_solver_type\'\npetsc_options_value = \'lu mumps\'\n'+block('TimeIntegrator',f'type = AqueductHHT\nbeta = {BETA}\ngamma = {GAMMA}\nhht_alpha = {ALPHA}\nhht_start_time = 1\ninactive_tsteps = 1\nalgebraic_variables = {quoted(lmnames)}')+block('TimeStepper','type = FunctionDT\nfunction = dt_clock'))
    text+=block('Outputs','csv = true\nexodus = true\nprint_linear_residuals = false')
    if mumps_memory_retry:
        # R02 failed FACTOR_OUTMEMORY before a usable Newton correction.
        # Bound per-rank factor memory on the selected32GB node; keep all physics.
        text=text.replace("petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'", "petsc_options = '-ksp_error_if_not_converged -ksp_converged_reason'\npetsc_options_iname = '-pc_type -pc_factor_mat_solver_type -mat_mumps_icntl_14 -mat_mumps_icntl_4 -mat_mumps_icntl_23'")
        text=text.replace("petsc_options_value = 'lu mumps'", "petsc_options_value = 'lu mumps 100 2 2048'")
    if diagnostics:
        text+=block('Outputs/profile','type = PerfGraphOutput\nlevel = 5\nheaviest_sections = 25\nheaviest_branch = true')
        text=text.replace('-mat_mumps_icntl_23\'', '-mat_mumps_icntl_23 -ksp_view_mat\'').replace('lu mumps 100 2 2048\'', 'lu mumps 100 2 2048 ::ascii_info_detail\'')
    (out/'aqueduct_full_candidate.i').write_text(text)
    identity=dict(case_id=c['case_id'],full_assembly=True,end_time_total_s=end_time,dynamic_duration_s=end_time-1,source_inp=c['source'],source_contract_sha256=sha(root/'assembly-contract.json'),elements={k:sum(len(g['ids']) for g in group.values() if g['signature'][0]==k) for k in ('HEX8','TRUSS','SHELL4')},nodes=len(xyz),mpc_rows=len(rows),embedded_nodes=len(embedded_nodes),embedded_equations=3*len(embedded_nodes),embedded_formulation='explicit conjugate multipliers; all physical mass/force retained; current algebraic rows',connector_count=4,contact_pairs=4,probes=probes,hht=dict(alpha=ALPHA,beta=BETA,gamma=GAMMA),first_candidate_user_regression_waiver='2026-10-24; seven not run, not passed',
                  deviations=['HEX8 full integration Bbar and six-DOF MITC4/SPOS shell accepted engineering replacements.', 'Custom consistent primal mortar Tie; actual slave membership preserved, coefficients differ from Abaqus; only fully bound secondary facets integrated.', 'Candidate normal-plane Tie adjustment, max .198 micrometre; exact ODB adjustment pending.', 'Explicit HHT alpha=-.41421 inferred from 2025 contact default; exact Abaqus2022 parameters unverified.', 'Base DOFs4-6 are absent on solid base nodes, so source zero rotation BCs have no active DOF; no physical rotary restraint is removed.', 'Unforced static stage one step instead of Abaqus automatic increments; source stabilization .0002 not implemented; zero displacement gives zero stabilization work.', 'Short first full-assembly candidate, not complete27.5s; no matching Abaqus output or engineering threshold.', 'Consistent truss line mass rho*A; original Abaqus reduced-integration/hourglass and source truss mass equivalence pending.'])
    identity['solver_memory_policy']=dict(workspace_estimate_increase_percent=100,max_working_memory_MB_per_rank=2048,MPI_ranks=4,diagnostic_print_level=2,error_on_linear_failure=True,physics_and_tolerances_unchanged=True) if mumps_memory_retry else None
    identity['ignore_zero_jacobian_entries']=ignore_zero_entries
    identity['hash_table_matrix_assembly']=hash_assembly
    identity['native_diagnostics']=dict(enabled=diagnostics,bar_fields=5,shell_inplane_qps=4,shell_layers=3,shell_state_components=['local_sxx','local_syy','local_szz','local_sxy','peq'],residual_tags=['TIME','NONTIME'],meaning='Unmodified complete tagged residuals; reconstruct normalized HHT before nodal BCs; not direct Abaqus RF') if diagnostics else None
    identity['artifacts']={p.name:sha(p) for p in out.iterdir()};(out/'full-input-identity.json').write_text(json.dumps(identity,indent=2)+'\n');return identity

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',required=True);p.add_argument('--motion',required=True);p.add_argument('--output',required=True);p.add_argument('--end-time',type=float,default=1.01);p.add_argument('--mumps-memory-retry',action='store_true');p.add_argument('--diagnostics',action='store_true');p.add_argument('--ignore-zero-entries',action='store_true');p.add_argument('--hash-assembly',action='store_true');a=p.parse_args();print(json.dumps(build(a.preparation,a.motion,a.output,a.end_time,a.mumps_memory_retry,a.diagnostics,a.ignore_zero_entries,a.hash_assembly),indent=2))
