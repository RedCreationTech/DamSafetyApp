#!/usr/bin/env python3
"""Build small independent interface fixtures; never a full aqueduct input."""
import argparse
import hashlib
import json
from pathlib import Path


def block(name, entries):
    return '['+name+']\n'+''.join('  '+line+'\n' for line in entries.splitlines())+'[]\n'


def connector_case(stiffness):
    fields=('ux','uy','uz','rx','ry','rz')
    text=block('Problem','kernel_coverage_check = false\nmaterial_coverage_check = false')
    text+=block('Mesh','type = GeneratedMesh\ndim = 1\nnx = 1\nxmax = .1\nallow_renumbering = false')
    text+=block('Variables',''.join(block(v,'') for v in fields))
    constraints=''
    for j,v in enumerate(fields):
        constraints+=block(v,f'type = LinearNodalConstraint\nvariable = {v}\nprimary = 0\nsecondary_node_ids = 1\nweights = 1\npenalty = {stiffness[j] if j<3 else 1e9}\nformulation = {"penalty" if j<3 else "kinematic"}')
    text+=block('Constraints',constraints)
    text+=block('BCs',''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = left\nvalue = {0 if j<3 else (.01,.02,.03)[j-3]}') for j,v in enumerate(fields)))
    text+=block('NodalKernels',''.join(block(v,f'type = UserForcingFunctorNodalKernel\nvariable = {v}\nboundary = right\nfunctor = {f}') for v,f in zip(fields,(stiffness[0]*.001,stiffness[1]*.002,stiffness[2]*.001))))
    text+=block('Postprocessors',''.join(block(v,f'type = NodalVariableValue\nvariable = {v}\nnodeid = 1') for v in fields))
    text+=block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner','type = Steady\nsolve_type = NEWTON\nnl_abs_tol = 1e-8\nnl_rel_tol = 1e-12\npetsc_options_iname = \'-pc_type\'\npetsc_options_value = \'lu\'')
    return text+block('Outputs','csv = true\nexodus = true')


def contact_case():
    mesh=''
    for name,z0,z1,id in (('bottom',0,.1,1),('top',.16,.26,2)):
        mesh+=block(name,f'type = GeneratedMeshGenerator\ndim = 3\nnx = 2\nny = 2\nnz = 1\nxmax = .2\nymax = .2\nzmin = {z0}\nzmax = {z1}\nelem_type = HEX8')
        mesh+=block(name+'_names',f'type = RenameBoundaryGenerator\ninput = {name}\nold_boundary = \'0 1 2 3 4 5\'\nnew_boundary = \'{name}_bottom {name}_back {name}_right {name}_front {name}_left {name}_top\'')
        mesh+=block(name+'_id',f'type = SubdomainIDGenerator\ninput = {name}_names\nsubdomain_id = {id}')
    mesh+=block('combine','type = MeshCollectionGenerator\ninputs = \'bottom_id top_id\'')
    text=block('Mesh',mesh+'patch_update_strategy = iteration\n')
    text+=block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Physics/SolidMechanics/QuasiStatic',block('all',"add_variables = true\nstrain = SMALL\nblock = '1 2'\nsave_in = 'force_x force_y force_z'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('force_x','force_y','force_z')))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nyoungs_modulus = 1e6\npoissons_ratio = 0')+block('stress','type = ComputeLinearElasticStress'))
    contact="primary = bottom_top\nsecondary = top_bottom\nmodel = coulomb\nformulation = mortar\nfriction_coefficient = .55\nc_normal = 100\nc_tangential = 100\nnormalize_c = true\nuse_dual = true\ncorrect_edge_dropping = true"
    text+=block('Contact',block('pair',contact))
    bcs=''
    for v in ('ux','uy','uz'):
        bcs+=block('bottom_'+v,f'type = DirichletBC\nvariable = {v}\nboundary = bottom_bottom\nvalue = 0')
    for v,fn in (('ux','if(t<1,0,.04*(t-1))'),('uy','0'),('uz','if(t<1,-.08*t,-.08)')):
        bcs+=block('top_'+v,f'type = FunctionDirichletBC\nvariable = {v}\nboundary = top_top\nfunction = \'{fn}\'')
    text+=block('BCs',bcs)
    text+=block('Postprocessors',''.join(block('top_'+v,f'type = NodalSum\nvariable = force_{v}\nboundary = top_top') for v in ('x','y','z')))
    text+=block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner','type = Transient\nsolve_type = NEWTON\nend_time = 2\ndt = .1\ndtmin = .0001\nnl_abs_tol = 1e-7\nnl_rel_tol = 1e-10\nnl_max_its = 40\nline_search = basic\npetsc_options_iname = \'-pc_type\'\npetsc_options_value = \'lu\'')
    return text+block('Outputs','csv = true\nexodus = true')


def build(contract,output):
    d=json.loads(Path(contract).read_text());output=Path(output)
    if d.get('case_id')!='AQ1023-FULL-ELC':raise ValueError('wrong source contract')
    behavior=d['interfaces']['connector_behaviors_si'];values=list(behavior.values())
    if len(values)!=4 or any(v!=values[0] for v in values):raise ValueError('fixture requires the four identical source connectors')
    k=[values[0][str(j)] for j in (1,2,3)]
    if k!=[605000.,605000.,21965000.]:raise ValueError('source connector parameters changed')
    if output.exists():raise ValueError('output exists')
    output.mkdir(parents=True)
    (output/'cartesian_align.i').write_text(connector_case(k))
    (output/'contact_open_close_slip.i').write_text(contact_case())
    shell=Path(__file__).resolve().parents[2]/'test/tests/aqueduct_interfaces/shell_bending.i'
    (output/'shell_bending.i').write_bytes(shell.read_bytes())
    summary=dict(status='interface-fixtures-only-not-full-aqueduct',no_solver_executed=True,contract_sha256=hashlib.sha256(Path(contract).read_bytes()).hexdigest(),artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.glob('*.i'))},connector_expected=dict(ux=.001,uy=.002,uz=.001,rx=.01,ry=.02,rz=.03),contact_initial_gap_m=.06,contact_friction=.55,contact_limitation='Mortar ideal Coulomb; source elastic-slip tolerance is not yet reproduced')
    (output/'interface-case-manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contract',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    print(json.dumps(build(a.contract,a.output),indent=2))
