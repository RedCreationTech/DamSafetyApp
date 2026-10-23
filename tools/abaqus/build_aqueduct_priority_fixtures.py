#!/usr/bin/env python3
"""Tie-independent three-axis embedding, shell-solid and shell plastic fixtures."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from abaqus2exodus import GlobalMesh, write_exodus
from build_aqueduct_interface_cases import block
from build_aqueduct_history_cases import execution, FIELDS
from map_aqueduct_interfaces import project_quad, rotation_operator, cross_matrix


def identity(out, expected, **extra):
    data=dict(expected=expected,**extra,artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()})
    (out/'fixture-identity.json').write_text(json.dumps(data,indent=2)+'\n')


def embedded(out):
    out.mkdir();gm=GlobalMesh(1e-10,merge_nodes=False)
    xyz=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1),
         (.25,.5,.5),(.75,.5,.5),(.5,.25,.5),(.5,.75,.5),(.5,.5,.25),(.5,.5,.75)]
    for n,p in enumerate(xyz,1):gm.add(('fixture',n),p)
    blocks={'host':[(1,list(range(1,9)))]};types={'host':'HEX8'}
    for j,axis in enumerate('xyz'):
        blocks['steel_'+axis]=[(j+2,[9+2*j,10+2*j])];types['steel_'+axis]='TRUSS'
    ns={'x0':[1,4,5,8],'x1':[2,3,6,7],'y0':[1,2,5,6],'y1':[3,4,7,8],'z0':[1,2,3,4],'z1':[5,6,7,8]}
    write_exodus(out/'fixture.e',gm,blocks,types,{},ns,{},'Three-axis embedded force-transfer fixture')
    with (out/'embedded-map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['secondary_node','primary_element','secondary_block','x_m','y_m','z_m']+[f'host_node_{j}' for j in range(8)])
        for j in range(3):
            for n in (9+2*j,10+2*j):w.writerow([n-1,0,j+2,*xyz[n-1],*range(8)])
    text=block('Mesh','type = FileMesh\nfile = fixture.e\nallow_renumbering = false')
    text+=block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Variables',''.join(block(v,'') for v in FIELDS[:3]))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('host',"block = host\nstrain = SMALL\nsave_in = 'host_rx host_ry host_rz'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('host_rx','host_ry','host_rz'))+
                block('area',"block = 'steel_x steel_y steel_z'\nfamily = MONOMIAL\norder = CONSTANT\ninitial_condition = 1e-6")+
                block('axial_stress',"block = 'steel_x steel_y steel_z'\nfamily = MONOMIAL\norder = CONSTANT"))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nblock = host\nyoungs_modulus = 1e6\npoissons_ratio = 0')+
                block('stress','type = ComputeLinearElasticStress\nblock = host')+
                block('steel',"type = AqueductSteelTruss\nblock = 'steel_x steel_y steel_z'\ndisplacements = 'ux uy uz'\nyoungs_modulus = 200e9\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1"))
    text+=block('Kernels',''.join(block(v,f"type = StressDivergenceTensorsTruss\nvariable = {v}\ncomponent = {j}\nblock = 'steel_x steel_y steel_z'\ndisplacements = 'ux uy uz'\narea = area\nuse_displaced_mesh = false") for j,v in enumerate(FIELDS[:3])))
    text+=block('AuxKernels',block('stress',"type = MaterialRealAux\nvariable = axial_stress\nproperty = axial_stress\nblock = 'steel_x steel_y steel_z'"))
    text+=block('Constraints',''.join(block(f'{axis}_{v}',f'type = AqueductEmbeddedConstraint\nvariable = {v}\nprimary_variable = {v}\nprimary = host\nsecondary = steel_{axis}\nallowed_primary_block_ids = 1\nsource_map = embedded-map.csv\nformulation = kinematic\npenalty = 1e10') for axis in 'xyz' for v in FIELDS[:3]))
    text+=block('BCs',''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = {axis}0\nvalue = 0') for axis,v in zip('xyz',FIELDS)))
    forces=(1100.,825.,550.);strains=np.array(forces)/1.1e6;expected={}
    text+=block('NodalKernels',''.join(block(v,f'type = UserForcingFunctorNodalKernel\nvariable = {v}\nboundary = {axis}1\nfunctor = {force/4}') for v,axis,force in zip(FIELDS,'xyz',forces)))
    post=''
    for n in [6,*range(8,14)]:
        for j,v in enumerate(FIELDS[:3]):
            name=f'n{n}_{v}';expected[name]=float(xyz[n][j]*strains[j])
            post+=block(name,f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for j,axis in enumerate('xyz'):
        expected['steel_'+axis]=float(200e9*strains[j]);expected['host_'+axis]=float(1e6*strains[j])
        post+=block('steel_'+axis,f'type = ElementAverageValue\nvariable = axial_stress\nblock = steel_{axis}')
        post+=block('host_'+axis,f'type = NodalSum\nvariable = host_r{axis}\nboundary = {axis}1')
    text+=block('Postprocessors',post)+block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner',"type = Steady\nsolve_type = NEWTON\nnl_abs_tol = 1e-7\nnl_rel_tol = 1e-12\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'")+block('Outputs','csv = true\nexodus = true')
    (out/'embedded_three_axis.i').write_text(text)
    identity(out,expected,reference='Independent affine energy K_j=1e6+200e9*1e-6*.5=1.1e6 N/m for every axis',limitations=['Three orthogonal bars; no oblique bar, shear or multiple host interface.','No source Tie adjustment or CDP material.'])


def shell_plastic(out):
    out.mkdir();text=block('Mesh','type = GeneratedMesh\ndim = 2\nnx = 1\nny = 1\nxmax = .2\nymax = .2\nelem_type = QUAD4')
    text+=block('GlobalParams',"displacements = 'ux uy uz'\nrotations = 'rx ry rz'")
    text+=block('Variables',''.join(block(v,'') for v in FIELDS))
    props={f'sxx_{j}':5*j for j in range(3)}|{f'peq_{j}':5*j+4 for j in range(3)}|{'syy':6,'szz':7}
    text+=block('AuxVariables',block('moment_ry','')+''.join(block(v,'family = MONOMIAL\norder = CONSTANT') for v in props))
    text+=block('Kernels',''.join(block(v,f'type = AqueductFlatShell\nvariable = {v}\ncomponent = {j}'+ ('\nsave_in = moment_ry' if v=='ry' else '')) for j,v in enumerate(FIELDS)))
    text+=block('Materials',block('shell','type = AqueductFlatShellStress\nyoungs_modulus = 200e9\npoissons_ratio = .3\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1\nthickness = .1\noffset = .5\nshear_correction = .8333333333333333\ndrilling_factor = .001'))
    text+=block('AuxKernels',''.join(block(v,f'type = MaterialStdVectorAux\nvariable = {v}\nproperty = aq_shell_state\nindex = {index}\nexecute_on = TIMESTEP_END') for v,index in props.items()))
    text+=block('BCs',block('ux','type = DirichletBC\nvariable = ux\nboundary = left\nvalue = 0')+block('uy','type = DirichletBC\nvariable = uy\nboundary = bottom\nvalue = 0')+
                ''.join(block(v,f"type = DirichletBC\nvariable = {v}\nboundary = 'left right'\nvalue = 0") for v in FIELDS[2:]))
    text+=block('Functions',block('load',"type = PiecewiseLinear\nx = '0 .5 1 1.5 2'\ny = '0 4500000 0 -5000000 0'"))
    text+=block('NodalKernels',block('load','type = UserForcingFunctorNodalKernel\nvariable = ux\nboundary = right\nfunctor = load'))
    post=''.join(block(v,f'type = ElementAverageValue\nvariable = {v}') for v in props)
    post+=block('right_u','type = PointValue\nvariable = ux\npoint = ".2 .1 0"')+block('top_v','type = PointValue\nvariable = uy\npoint = ".1 .2 0"')
    post+=block('right_moment','type = NodalSum\nvariable = moment_ry\nboundary = right')
    text+=block('Postprocessors',post)+execution(2,.05)
    (out/'shell_plastic_cycle.i').write_text(text)
    H=2.96e9;history=[]
    for t,s,p,q in [(.5,450e6,50e6/H,50e6/H),(1,0,50e6/H,50e6/H),(1.5,-500e6,0,100e6/H),(2,0,0,100e6/H)]:
        row=dict(time=t,right_u=.2*(s/200e9+p),top_v=.2*(-.3*s/200e9-.5*p),right_moment=-.05*.02*s,syy=0.,szz=0.)
        row.update({f'sxx_{j}':s for j in range(3)});row.update({f'peq_{j}':q for j in range(3)});history.append(row)
    tolerances={**{k:10. for k in props if k.startswith('s')},**{f'peq_{j}':1e-9 for j in range(3)},'right_u':1e-9,'top_v':1e-9,'right_moment':.001}
    identity(out,{},expected_history=history,absolute_tolerances=tolerances,reference='Uniaxial isotropic hardening with plane stress and zero curvature; plastic lateral strain=-p/2; offset moment=-hF/2',limitations=['Curvature and transverse motion prescribed zero; axial and lateral in-plane response are free.','Not plastic bending, torsion, dynamic free response or Abaqus S4R equivalence.'])


def shell_solid(out):
    out.mkdir();gm=GlobalMesh(1e-10,merge_nodes=False)
    xyz=np.array([[0,0,-.2],[.2,0,-.2],[.2,.2,-.2],[0,.2,-.2],[0,0,0],[.2,0,0],[.2,.2,0],[0,.2,0],
                  [.2,0,0],[.4,0,0],[.4,.2,0],[.2,.2,0]])
    for n,p in enumerate(xyz,1):gm.add(('fixture',n),tuple(p))
    write_exodus(out/'fixture.e',gm,{'host':[(1,list(range(1,9)))],'shell':[(2,[9,10,11,12])]}, {'host':'HEX8','shell':'QUAD4'},{},
                 {'host_all':list(range(1,9)),'shell_right':[10,11]}, {},'Shell-edge to solid-face force and moment fixture')
    # Source-restricted reference face. Exactly the same P/A operator used by the converter.
    face_nodes=[1,2,6,5];face=xyz[face_nodes];rows=[]
    for node in (8,11):
        p=project_quad(xyz[node],face);A=rotation_operator(p['gradients'],p['normal']);r=xyz[node]-p['point']
        P=np.array([N*np.eye(3)-cross_matrix(r)@a for N,a in zip(p['weights'],A)])
        for j,v in enumerate(FIELDS):
            terms=[dict(node=node,variable=v,coefficient=1.)]
            op=P[:,j,:] if j<3 else A[:,j-3,:]
            for n,coeff in zip(face_nodes,op):
                for k,c in enumerate(coeff):
                    if abs(c)>1e-14:terms.append(dict(node=n,variable=FIELDS[k],coefficient=-float(c)))
            rows.append(dict(anchor=node,component=j,terms=terms))
    (out/'coupling-rows.json').write_text(json.dumps(dict(rows=rows,coordinates=xyz.tolist(),loads_at_2=dict(force=[10000.,2000.,500.],moment=[10.,1000.,5.])),indent=2)+'\n')
    text=block('Mesh',block('file','type = FileMeshGenerator\nfile = fixture.e')+block('anchors',"type = AqueductConstraintNodesGenerator\ninput = file\nnodes = '8 11'\nblock_id = 10\nblock_name = anchors")+'allow_renumbering = false')
    text+=block('Problem','material_coverage_check = SKIP_LIST\nmaterial_coverage_block_list = 10')
    text+=block('GlobalParams',"displacements = 'ux uy uz'\nrotations = 'rx ry rz'")
    text+=block('Variables',''.join(block(v,"block = 'host shell'" if j<3 else 'block = shell') for j,v in enumerate(FIELDS))+''.join(block('lm'+str(j),'block = anchors') for j in range(6)))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('host','block = host\nstrain = SMALL'))
    text+=block('Kernels',''.join(block(v,f'type = AqueductFlatShell\nvariable = {v}\ncomponent = {j}\nblock = shell') for j,v in enumerate(FIELDS))+''.join(block('lm'+str(j),f"type = NullKernel\nvariable = lm{j}\nblock = anchors\ndisplacements = ''") for j in range(6)))
    text+=block('Materials',block('host_elastic','type = ComputeIsotropicElasticityTensor\nblock = host\nyoungs_modulus = 200e9\npoissons_ratio = 0')+block('host_stress','type = ComputeLinearElasticStress\nblock = host')+
                block('shell','type = AqueductFlatShellStress\nblock = shell\nyoungs_modulus = 200e9\npoissons_ratio = 0\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1\nthickness = .1\noffset = .5\nshear_correction = .8333333333333333\ndrilling_factor = .001'))
    text+=block('Constraints',''.join(block(f'c_{r["anchor"]}_{r["component"]}',f'type = AqueductLinearMPC\nvariable = ux\nanchor_node = {r["anchor"]}\nterm_nodes = \'{" ".join(str(t["node"]) for t in r["terms"])}\'\nterms = \'{" ".join(t["variable"] for t in r["terms"])}\'\ncoefficients = \'{" ".join(str(t["coefficient"]) for t in r["terms"])}\'\nmultiplier = lm{r["component"]}') for r in rows))
    text+=block('BCs',''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = host_all\nvalue = 0') for v in FIELDS[:3]))
    text+=block('Functions',''.join(block(v,f'type = PiecewiseLinear\nx = \'0 1 2\'\ny = \'0 {500 if v=="ry" else 0} {load}\'') for v,load in zip(FIELDS,[5000,1000,250,5,500,2.5])))
    text+=block('NodalKernels',''.join(block(v,f'type = UserForcingFunctorNodalKernel\nvariable = {v}\nboundary = shell_right\nfunctor = {v}') for v in FIELDS))
    post=''.join(block(f'lm_{r["anchor"]}_{r["component"]}',f'type = NodalVariableValue\nvariable = lm{r["component"]}\nnodeid = {r["anchor"]}') for r in rows)
    post+=''.join(block(f'anchor_{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}') for n in (8,11) for v in FIELDS)
    post+=''.join(block(name,f'type = PointValue\nvariable = {v}\npoint = ".4 .1 0"') for name,v in [('tip_ux','ux'),('tip_uz','uz'),('tip_ry','ry')])
    text+=block('Postprocessors',post)+execution(2,.1)
    (out/'shell_solid_transfer.i').write_text(text)
    identity(out,{},expected_history=[dict(time=1.,tip_ux=3e-6,tip_uz=-6e-6,tip_ry=6e-5)],absolute_tolerances=dict(tip_ux=1e-12,tip_uz=1e-12,tip_ry=1e-11),
             reference='Fixed solid reference face, actual six-DOF MPC shell edge; pure end My=1000 at t1; six-component equilibrium at t2',limitations=['Solid nodes fixed: tests force/moment transfer to a rigid support, not compliant solid response.','Custom surface-gradient rotation is not asserted Abaqus shell-solid discretization equivalent.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
    embedded(out/'embedded');shell_plastic(out/'shell-plastic');shell_solid(out/'shell-solid')
