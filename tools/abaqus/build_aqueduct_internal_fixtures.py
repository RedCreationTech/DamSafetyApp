#!/usr/bin/env python3
"""Build source-physics interface fixtures; never a complete aqueduct input."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from abaqus2exodus import GlobalMesh,write_exodus
from build_aqueduct_interface_cases import block
from build_aqueduct_priority_fixtures import identity

FIELDS=('ux','uy','uz')
XYZ=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
SETS={'x0':[1,4,5,8],'x1':[2,3,6,7],'y0':[1,2,5,6],'y1':[3,4,7,8],'z0':[1,2,3,4],'z1':[5,6,7,8]}


def cube(out):
    gm=GlobalMesh(1e-10,merge_nodes=False)
    for n,p in enumerate(XYZ,1):gm.add(('fixture',n),p)
    write_exodus(out/'fixture.e',gm,{'solid':[(1,list(range(1,9)))]},{'solid':'HEX8'},{},SETS,{},'AQ1023 source physics single HEX8 reference')
    return block('Mesh','type = FileMesh\nfile = fixture.e\nallow_renumbering = false')+block('GlobalParams',"displacements = 'ux uy uz'")


def execution(end,dt):
    return block('Preconditioning',block('smp','type = SMP\nfull = true'))+block('Executioner',f"type = Transient\nsolve_type = NEWTON\nend_time = {end}\ndt = {dt}\ndtmin = 1e-8\nnl_abs_tol = 1e-7\nnl_rel_tol = 1e-11\nnl_max_its = 40\nline_search = basic\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'")+block('Outputs','csv = true\nexodus = true')


def steel_oracle(strains):
    E=200e9;fy=400e6;fu=696e6;H=(fu-fy)/.1;p=eq=0.;result=[]
    for e in strains:
        trial=E*(e-p);sign=1 if trial>=0 else -1;delta=max(0.,(abs(trial)-min(fu,fy+H*eq))/(E+(H if eq<.1 else 0)))
        if eq+delta>.1:delta=max(0.,(abs(trial)-fu)/E)
        p+=sign*delta;eq+=delta;s=E*(e-p)
        result.append(dict(strain=e,sxx=s,peq=eq,top_y=-.3*s/E-.5*p,right_force=s))
    return result


def steel(out):
    out.mkdir();text=cube(out)
    text+=block('Physics/SolidMechanics/QuasiStatic',block('all',"add_variables = true\nblock = solid\nstrain = SMALL\nincremental = true\nvolumetric_locking_correction = true\nsave_in = 'fx fy fz'\ngenerate_output = 'stress_xx stress_yy stress_zz'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('fx','fy','fz')))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nblock = solid\nyoungs_modulus = 200e9\npoissons_ratio = .3')+block('stress','type = ComputeMultipleInelasticStress\nblock = solid\ninelastic_models = steel_plastic\ntangent_operator = nonlinear')+block('steel_plastic','type = IsotropicPlasticityStressUpdate\nblock = solid\nyield_stress = 400e6\nhardening_function = source_hardening'))
    times=[0.,.1,.2,.5,.6,.7,1.,1.5,2.];strains=[0.,.002,.005,.12,.118,.115,0.,-.12,0.]
    text+=block('Functions',block('source_hardening',"type = PiecewiseLinear\nx = '0 .1'\ny = '400e6 696e6'\nextrap = false")+block('axial',"type = PiecewiseLinear\nx = '"+' '.join(map(str,times))+"'\ny = '"+' '.join(map(str,strains))+"'"))
    text+=block('BCs',block('right','type = FunctionDirichletBC\nvariable = ux\nboundary = x1\nfunction = axial')+''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = {axis}0\nvalue = 0') for axis,v in zip('xyz',FIELDS)))
    post=''.join(block(v,f'type = ElementAverageValue\nvariable = stress_{v}\nblock = solid') for v in ('xx','yy','zz'))
    post+=block('peq','type = ElementAverageMaterialProperty\nmat_prop = effective_plastic_strain\nblock = solid')+block('top_y','type = NodalVariableValue\nvariable = uy\nnodeid = 2')+block('right_force','type = NodalSum\nvariable = fx\nboundary = x1')
    text+=block('Postprocessors',post)+execution(2,.05);(out/'steel_solid_cycle.i').write_text(text)
    oracle=steel_oracle(strains);expected=[dict(time=t,**{'xx':o['sxx'],'yy':0.,'zz':0.,'peq':o['peq'],'top_y':o['top_y'],'right_force':o['right_force']}) for t,o in zip(times,oracle) if t>0]
    identity(out,{},expected_history=expected,scope='3D HEX8 axial strain-controlled with free lateral closure; source J2 table including constant end extrapolation',limitations=['Not a full aqueduct or CDP test.','No Abaqus steel solid result; small-strain large diagnostic strain deliberately crosses the last source plastic row.'])


def contact(out):
    out.mkdir();mesh=''
    for name,lo,hi,z0,z1,bid in [('bottom',0,.24,0,.1,1),('top',.02,.22,.16,.26,2)]:
        mesh+=block(name,f'type = GeneratedMeshGenerator\ndim = 3\nnx = 2\nny = 2\nnz = 1\nxmin = {lo}\nxmax = {hi}\nymin = {lo}\nymax = {hi}\nzmin = {z0}\nzmax = {z1}\nelem_type = HEX8')
        previous=name
        if name=='top':
            mesh+=block('top_ids',"type = RenameBoundaryGenerator\ninput = top\nold_boundary = '0 1 2 3 4 5'\nnew_boundary = '10 11 12 13 14 15'");previous='top_ids'
        ids='0 1 2 3 4 5' if name=='bottom' else '10 11 12 13 14 15'
        mesh+=block(name+'_names',f"type = RenameBoundaryGenerator\ninput = {previous}\nold_boundary = '{ids}'\nnew_boundary = '{name}_bottom {name}_back {name}_right {name}_front {name}_left {name}_top'")
        mesh+=block(name+'_id',f'type = SubdomainIDGenerator\ninput = {name}_names\nsubdomain_id = {bid}')
    mesh+=block('combine',"type = MeshCollectionGenerator\ninputs = 'bottom_id top_id'")+block('secondary',"type = LowerDBlockFromSidesetGenerator\ninput = combine\nsidesets = top_bottom\nnew_block_id = 101\nnew_block_name = secondary_lower")+block('primary',"type = LowerDBlockFromSidesetGenerator\ninput = secondary\nsidesets = bottom_top\nnew_block_id = 100\nnew_block_name = primary_lower")
    text=block('Mesh',mesh+'allow_renumbering = false\npatch_update_strategy = iteration')+block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Variables',''.join(block(v,'block = secondary_lower\nuse_dual = true\nscaling = .01') for v in ('normal_lm','tangent_1','tangent_2')))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('all',"add_variables = true\nblock = '1 2'\nstrain = SMALL\nsave_in = 'fx fy fz'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('fx','fy','fz','tx','ty','tz')))
    text+=block('Materials',block('elastic',"type = ComputeIsotropicElasticityTensor\nblock = '1 2'\nyoungs_modulus = 1e6\npoissons_ratio = 0")+block('stress',"type = ComputeLinearElasticStress\nblock = '1 2'"))
    common="primary_boundary = bottom_top\nsecondary_boundary = top_bottom\nprimary_subdomain = primary_lower\nsecondary_subdomain = secondary_lower\ncorrect_edge_dropping = true\nuse_displaced_mesh = true"
    text+=block('UserObjects',block('elastic_slip',common+'\ntype = AqueductElasticSlipUserObject\nlm_variable_normal = normal_lm\nlm_variable_tangential_one = tangent_1\nlm_variable_tangential_two = tangent_2\nsecondary_variable = ux\ndisp_x = ux\ndisp_y = uy\ndisp_z = uz'))
    constraint=block('law',common+'\ntype = AqueductElasticSlipConstraint\nvariable = normal_lm\ndisp_x = ux\ndisp_y = uy\ndisp_z = uz\nfriction_lm = tangent_1\nfriction_lm_dir = tangent_2\nweighted_gap_uo = elastic_slip\nweighted_velocities_uo = elastic_slip\ncompute_primal_residuals = false\nnormalize_c = true\nc = 1e6\nmu = .55\ncritical_slip = .0005')
    for j,v in enumerate(FIELDS):
        for name,typ,lm,uo,extra in [('normal','NormalMortarMechanicalContact','normal_lm','weighted_gap_uo',''),('t1','TangentialMortarMechanicalContact','tangent_1','weighted_velocities_uo',''),('t2','TangentialMortarMechanicalContact','tangent_2','weighted_velocities_uo','\ndirection = direction_2')]:
            constraint+=block(name+'_'+v,common+f'\ntype = {typ}\nvariable = {lm}\nsecondary_variable = {v}\ncomponent = {"xyz"[j]}\ncompute_lm_residuals = false\n{uo} = elastic_slip'+extra)
    text+=block('Constraints',constraint)
    aux_common = 'primary_boundary = bottom_top\nsecondary_boundary = top_bottom\nuse_displaced_mesh = true'
    text+=block('AuxKernels',''.join(block(v,aux_common+f'\ntype = MortarFrictionalPressureVectorAux\nvariable = {v}\ntangent_one = tangent_1\ntangent_two = tangent_2\ncomponent = {j}\nboundary = top_bottom') for j,v in enumerate(('tx','ty','tz'))))
    ts='0 1 2 3 4 5 6 7 8 9 10';hist={'ux':[0,0,.0002,.0002,.0012,.0008,0,0,.0006,.0006,.0007],'uy':[0,0,.00015,.00015,.0009,.0006,0,0,.0008,.0008,.0009],'uz':[0,-.06001,-.06001,-.06002,-.06002,-.06002,-.06002,0,0,-.06001,-.06001]}
    text+=block('Functions',''.join(block('load_'+v,"type = PiecewiseLinear\nx = '"+ts+"'\ny = '"+' '.join(map(str,values))+"'") for v,values in hist.items()))
    bcs=''
    for v in FIELDS:
        bcs+=block('bottom_'+v,f"type = DirichletBC\nvariable = {v}\nboundary = 'bottom_bottom bottom_top bottom_left bottom_right bottom_back bottom_front'\nvalue = 0")+block('top_'+v,f'type = FunctionDirichletBC\nvariable = {v}\nboundary = top_top\nfunction = load_{v}')
    text+=block('BCs',bcs)
    post=''
    for n in range(18,27):
        for v in (*FIELDS,'normal_lm','tx','ty','tz'):
            post+=block(f'n{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for v in ('fx','fy','fz'):post+=block('top_'+v,f'type = NodalSum\nvariable = {v}\nboundary = top_top')
    text+=block('Postprocessors',post)+execution(10,.1);(out/'contact_elastic_slip.i').write_text(text)
    identity(out,{},scope='native hard-normal mortar + source elastic-slip law; small-deformation planar fixture',initial_gap_m=.06,friction_coefficient=.55,critical_slip_m=.0005,source_slip_tolerance=.005,fixture_secondary_edge_length_m=.1,contact_libmesh_node_ids=list(range(18,27)),limitations=['No Abaqus reference or full source contact discretization equivalence.','Source characteristic contact length is reported in Abaqus DAT; this patch uses its own known .1m square faces.','Rotating local tangents, general large sliding and partition changes/restart remain separate validation.'])


def hex_matrices(E=30e9,nu=.2,rho=2640.):
    xyz=np.array(XYZ,float);sign=2*xyz-1;M=np.zeros((24,24));K=M.copy();lam=E*nu/((1+nu)*(1-2*nu));G=E/(2*(1+nu));C=np.zeros((6,6));C[:3,:3]=lam;C[np.diag_indices(3)]+=2*G;C[3:,3:]=np.eye(3)*G
    for q in np.array(np.meshgrid(*[[-1/np.sqrt(3),1/np.sqrt(3)]]*3)).T.reshape(-1,3):
        N=np.prod(1+sign*q,axis=1)/8;grad=np.zeros((8,3))
        for j in range(3):grad[:,j]=sign[:,j]*np.prod(1+sign[:,[k for k in range(3) if k!=j]]*q[[k for k in range(3) if k!=j]],axis=1)/4
        B=np.zeros((6,24))
        for i,g in enumerate(grad):
            B[0,3*i]=g[0];B[1,3*i+1]=g[1];B[2,3*i+2]=g[2]
            B[3,3*i:3*i+3]=[g[1],g[0],0];B[4,3*i:3*i+3]=[0,g[2],g[1]];B[5,3*i:3*i+3]=[g[2],0,g[0]]
            for j in range(3):B[:3,3*i+j]+=(sign[i,j]/4-g[j])/3
        K+=B.T@C@B/8
        for j in range(3):M[j::3,j::3]+=rho*np.outer(N,N)/8
    return M,K


def motion(out,motion_directory):
    out.mkdir();text=cube(out)
    for axis in 'xyz':shutil.copyfile(Path(motion_directory)/('motion-'+axis+'.csv'),out/('motion-'+axis+'.csv'))
    shutil.copyfile(Path(motion_directory)/'motion-identity.json',out/'motion-identity.json')
    text+=block('Physics/SolidMechanics/QuasiStatic',block('all',"add_variables = true\nblock = solid\nstrain = SMALL\nvolumetric_locking_correction = true\nsave_in = 'fx fy fz'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('fx','fy','fz','vx','vy','vz','ax','ay','az')))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nblock = solid\nyoungs_modulus = 30e9\npoissons_ratio = .2')+block('stress','type = ComputeLinearElasticStress\nblock = solid')+block('rho',"type = GenericConstantMaterial\nblock = solid\nprop_names = density\nprop_values = 2640"))
    functions=block('dynamic_on',"type = ParsedFunction\nexpression = 't>1'")+block('dt_clock',"type = ParsedFunction\nexpression = 'if(t<1,1,.0025)'")
    for axis in 'xyz':
        for component in ('acceleration','velocity','displacement'):
            prefix = {'acceleration': 'a', 'velocity': 'v', 'displacement': 'u'}[component]
            functions+=block('motion_'+prefix+axis,f'type = IntegratedAcceleration\ndata_file = motion-{axis}.csv\ncomponent = {component}\nstart_time = 1')
    text+=block('Functions',functions)
    text+=block('Kernels',''.join(block('inertia_'+v,f'type = InertialForce\nvariable = {v}\nblock = solid\nbeta = .25\ngamma = .5\nvelocity = v{axis}\nacceleration = a{axis}\nuse_displaced_mesh = false\nsave_in = f{axis}') for v,axis in zip(FIELDS,'xyz')))
    aux=''
    for v,axis in zip(FIELDS,'xyz'):
        aux+=block('accel_'+axis,f'type = NewmarkAccelAux\nvariable = a{axis}\ndisplacement = {v}\nvelocity = v{axis}\nbeta = .25\nexecute_on = timestep_end')+block('velocity_'+axis,f'type = NewmarkVelAux\nvariable = v{axis}\nacceleration = a{axis}\ngamma = .5\nexecute_on = timestep_end')
    text+=block('AuxKernels',aux)
    bcs=''
    for v,axis in zip(FIELDS,'xyz'):
        bcs+=block('fixed_'+axis,f'type = DirichletBC\nvariable = {v}\nboundary = z0\nvalue = 0')+block('motion_'+axis,f'type = PresetAcceleration\nvariable = {v}\nboundary = z0\nfunction = motion_a{axis}\nvelocity = v{axis}\nacceleration = a{axis}\nbeta = .25')
    text+=block('BCs',bcs)
    enable=' '.join(['Kernels::inertia_'+v for v in FIELDS]+['BCs::motion_'+a for a in 'xyz']);disable=' '.join('BCs::fixed_'+a for a in 'xyz')
    text+=block('Controls',block('phase',f"type = ConditionalFunctionEnableControl\nconditional_function = dynamic_on\nenable_objects = '{enable}'\ndisable_objects = '{disable}'\nexecute_on = 'initial timestep_begin'"))
    post=''
    for n in range(8):
        for v in (*FIELDS,'vx','vy','vz','ax','ay','az'):post+=block(f'n{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for axis in 'xyz':
        post+=block('base_force_'+axis,f'type = NodalSum\nvariable = f{axis}\nboundary = z0')
        for c in 'auv':post+=block('ref_'+c+axis,f'type = FunctionValuePostprocessor\nfunction = motion_{c}{axis}')
    text+=block('Postprocessors',post)
    text+=block('Preconditioning',block('smp','type = SMP\nfull = true'))+block('Executioner',"type = Transient\nsolve_type = NEWTON\nend_time = 1.05\ndtmin = 1e-8\nnl_abs_tol = 1e-6\nnl_rel_tol = 1e-11\nline_search = basic\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'\n"+block('TimeIntegrator','type = NewmarkBeta\nbeta = .25\ngamma = .5\ninactive_tsteps = 1')+block('TimeStepper','type = FunctionDT\nfunction = dt_clock'))+block('Outputs','csv = true\nexodus = true')
    (out/'three_axis_history.i').write_text(text);M,K=hex_matrices();np.savez(out/'independent-hex8-matrices.npz',M=M,K=K,coordinates=np.array(XYZ))
    identity(out,{},scope='1s unforced static to three-axis source-baseline-corrected base acceleration; 50ms free elastic HEX8 motion',reference='independent full 24-DOF HEX8 Bbar stiffness and consistent mass with discrete Newmark .25/.5 recurrence',limitations=['Not the complete aqueduct or CDP validation.','Diagnostic uses Newmark alpha0; source 2022 contact-dependent HHT-MD is not asserted reproduced.','Baseline interval start convention remains explicit with sensitivity; no Abaqus corrected record supplied.'])


def build(output,motion_directory):
    out=Path(output)
    if out.exists():raise ValueError('output exists')
    out.mkdir(parents=True);steel(out/'steel-solid');contact(out/'elastic-slip-contact');motion(out/'three-axis-history',motion_directory)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);p.add_argument('--motion',required=True);a=p.parse_args();build(a.output,a.motion)
