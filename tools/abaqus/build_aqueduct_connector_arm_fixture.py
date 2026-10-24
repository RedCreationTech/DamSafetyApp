#!/usr/bin/env python3
"""Static Cartesian+Align reference-arm fixture; supports have no bulk stiffness/mass."""
import argparse,json
from pathlib import Path
import numpy as np
from abaqus2exodus import GlobalMesh,write_exodus
from build_aqueduct_interface_cases import block
from build_aqueduct_full_input import quoted,FIELDS
from map_aqueduct_interfaces import cartesian_connector_terms


def build(output):
    out=Path(output);out.mkdir();gm=GlobalMesh(1e-10,merge_nodes=False)
    gm.add(('fixture',1),(0,0,0));gm.add(('fixture',2),(0,0,.1))
    write_exodus(out/'fixture.e',gm,{'reference':[(1,[1,2])]},{'reference':'TRUSS'},{},{'fixed':[1],'loaded':[2]}, {},'Geometric support edge only: no bulk stiffness or mass')
    text=block('Mesh',block('source','type = FileMeshGenerator\nfile = fixture.e')+block('support','type = AqueductConstraintNodesGenerator\ninput = source\nnodes = 0\nblock_id = 100\nblock_name = multiplier_support')+'allow_renumbering = false\nparallel_type = replicated')
    text+=block('Problem','material_coverage_check = false')
    text+=block('Variables',''.join(block(v,'block = reference') for v in FIELDS)+''.join(block('force_'+a,'block = multiplier_support')+block('align_'+a,'block = multiplier_support') for a in 'xyz'))
    text+=block('Kernels',''.join(block(v,f'type = NullKernel\nvariable = {v}\nblock = reference\njacobian_fill = 0') for v in FIELDS)+''.join(block(v+'_'+a,f'type = NullKernel\nvariable = {v}_{a}\nblock = multiplier_support\njacobian_fill = 0') for a in 'xyz' for v in ('force','align')))
    constraints='';K=[605000.,605000.,21965000.]
    for a,k,terms in zip('xyz',K,cartesian_connector_terms(0,1,[0,0,.1])):
        constraints+=block('spring_'+a,f'type = AqueductLinearMPC\nvariable = ux\nanchor_node = 0\nterm_nodes = {quoted(t["node"] for t in terms)}\nterms = {quoted(t["variable"] for t in terms)}\ncoefficients = {quoted(t["coefficient"] for t in terms)}\nmultiplier = force_{a}\nmultiplier_compliance = {1/k}')
        constraints+=block('align_'+a,f'type = AqueductLinearMPC\nvariable = ux\nanchor_node = 0\nterm_nodes = \'0 1\'\nterms = \'r{a} r{a}\'\ncoefficients = \'-1 1\'\nmultiplier = align_{a}')
    text+=block('Constraints',constraints)
    force=np.array([12.,-18.,9.]);rotation=np.array([.001,-.002,.003]);functions=''
    for a,f,w in zip('xyz',force,rotation):functions+=block('load_'+a,f'type = ParsedFunction\nexpression = \'{f}*t\'')+block('rotation_'+a,f'type = ParsedFunction\nexpression = \'{w}*t\'')
    text+=block('Functions',functions)+block('NodalKernels',''.join(block('load_'+a,f'type = UserForcingFunctorNodalKernel\nvariable = u{a}\nboundary = loaded\nfunctor = load_{a}') for a in 'xyz'))
    text+=block('BCs',''.join(block('fixed_'+a,f'type = DirichletBC\nvariable = u{a}\nboundary = fixed\nvalue = 0')+block('rotation_'+a,f'type = FunctionDirichletBC\nvariable = r{a}\nboundary = fixed\nfunction = rotation_{a}') for a in 'xyz'))
    text+=block('AuxVariables',''.join(block('reaction_'+v,'block = reference') for v in FIELDS))
    text+=block('AuxKernels',''.join(block(v,f'type = TagVectorAux\nvariable = reaction_{v}\nv = {v}\nvector_tag = nontime\nremove_variable_scaling = true\nblock = reference') for v in FIELDS))
    post=''
    for v in FIELDS:post+=block('tip_'+v,f'type = NodalVariableValue\nvariable = {v}\nnodeid = 1')+block('reaction_'+v,f'type = NodalVariableValue\nvariable = reaction_{v}\nnodeid = 0')
    for a in 'xyz':post+=block('force_'+a,f'type = NodalVariableValue\nvariable = force_{a}\nnodeid = 0')+block('align_'+a,f'type = NodalVariableValue\nvariable = align_{a}\nnodeid = 0')
    text+=block('Postprocessors',post)+block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner',"type = Transient\nsolve_type = NEWTON\nend_time = 1\ndt = .25\nnl_abs_tol = 1e-10\nnl_rel_tol = 1e-12\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'")+block('Outputs','csv = true\nexodus = true')
    (out/'cartesian_arm.i').write_text(text);history=[]
    for t in np.arange(0,1.001,.25):
        F=t*force;theta=t*rotation;u=np.cross(theta,[0,0,.1])+F/np.array(K);moment=-np.cross([0,0,.1],F)
        row={'time':float(t)}
        for j,a in enumerate('xyz'):row.update({f'tip_u{a}':float(u[j]),f'tip_r{a}':float(theta[j]),f'force_{a}':float(F[j]),f'align_{a}':0.,f'reaction_u{a}':float(-F[j]),f'reaction_r{a}':float(moment[j])})
        history.append(row)
    (out/'fixture-identity.json').write_text(json.dumps(dict(expected_history=history,reference='Independent direct s=ub-ua-theta_a cross d0; force=K*s; node-a moment=-d0 cross force; Align=theta_b-theta_a',reference_arm_m=[0,0,.1],stiffness_N_per_m=K,reference_geometry_only_no_bulk_stiffness_or_mass=True,limitations=['First-order small rotations, default initial global axes; not finite-rotation or Abaqus2022 measured reference.']),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();build(a.output)
