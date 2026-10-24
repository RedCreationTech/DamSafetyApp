#!/usr/bin/env python3
"""Three-direction dynamic embedded mass/HHT fixture with independent 60-DOF solve."""
import argparse,csv,json,shutil
from pathlib import Path
import numpy as np
from build_aqueduct_priority_fixtures import embedded
from build_aqueduct_internal_fixtures import hex_matrices
from build_aqueduct_interface_cases import block
from build_aqueduct_full_input import quoted,ALPHA,BETA,GAMMA


def build(output):
    out=Path(output);embedded(out)
    xyz=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[.25,.5,.5],[.75,.5,.5],[.5,.25,.5],[.5,.75,.5],[.5,.5,.25],[.5,.5,.75]])
    text=block('Mesh',block('source','type = FileMeshGenerator\nfile = fixture.e')+block('support',"type = AqueductConstraintNodesGenerator\ninput = source\nnodes = '8 9 10 11 12 13'\nblock_id = 10\nblock_name = support")+'allow_renumbering = false\nparallel_type = replicated')
    text+=block('Problem','material_coverage_check = false')+block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Variables',''.join(block(v,"block = '1 2 3 4'") for v in ('ux','uy','uz'))+''.join(block('em_'+a,'block = support') for a in 'xyz'))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('host','block = host\nstrain = SMALL\nvolumetric_locking_correction = true'))
    text+=block('AuxVariables',block('area',"block = '2 3 4'\nfamily = MONOMIAL\norder = CONSTANT\ninitial_condition = 1e-6"))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nblock = host\nyoungs_modulus = 1e6\npoissons_ratio = 0')+block('stress','type = ComputeLinearElasticStress\nblock = host')+block('host_mass','type = GenericConstantMaterial\nblock = host\nprop_names = density\nprop_values = 1')+block('steel',"type = AqueductSteelTruss\nblock = '2 3 4'\nyoungs_modulus = 200e9\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1")+block('line_mass',"type = GenericConstantMaterial\nblock = '2 3 4'\nprop_names = density\nprop_values = .00785"))
    kernels=''
    for j,v in enumerate(('ux','uy','uz')):
        kernels+=block('bar_'+v,f'type = StressDivergenceTensorsTruss\nvariable = {v}\ncomponent = {j}\nblock = \'2 3 4\'\narea = area\nuse_displaced_mesh = false')
        kernels+=block('mass_'+v,f'type = AqueductHHTInertia\nvariable = {v}\nblock = \'1 2 3 4\'\nhht_alpha = {ALPHA}\nhht_start_time = 0\nuse_displaced_mesh = false')
        kernels+=block('em_'+v[-1],f'type = NullKernel\nvariable = em_{v[-1]}\nblock = support\ndisplacements = \'\'\njacobian_fill = 0')
    text+=block('Kernels',kernels)+block('Constraints',block('embed',"type = AqueductEmbeddedMPC\nvariable = ux\nmultipliers = 'em_x em_y em_z'\nsource_map = embedded-map.csv\nallowed_primary_block_ids = 1"))
    text+=block('BCs',''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = {a}0\nvalue = 0') for v,a in zip(('ux','uy','uz'),'xyz')))
    text+=block('Functions',''.join(block('load_'+a,f'type = ParsedFunction\nexpression = \'{f/4}*t\'') for a,f in zip('xyz',[1100,825,550])))
    text+=block('NodalKernels',''.join(block(v,f'type = UserForcingFunctorNodalKernel\nvariable = {v}\nboundary = {a}1\nfunctor = load_{a}') for v,a in zip(('ux','uy','uz'),'xyz')))
    post=''
    for n in range(14):
        for v in ('ux','uy','uz'):post+=block(f'n{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for n in range(8,14):
        for a in 'xyz':post+=block(f'em_{n}_{a}',f'type = NodalVariableValue\nvariable = em_{a}\nnodeid = {n}')
    text+=block('Postprocessors',post)+block('Preconditioning',block('smp','type = SMP\nfull = true'))+block('Executioner',"type = Transient\nsolve_type = NEWTON\nend_time = .03\ndt = .01\nnl_abs_tol = 1e-9\nnl_rel_tol = 1e-11\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'\n"+block('TimeIntegrator',f'type = AqueductHHT\nhht_alpha = {ALPHA}\nbeta = {BETA}\ngamma = {GAMMA}\nalgebraic_variables = \'em_x em_y em_z\''))+block('Outputs','csv = true\nexodus = true')
    (out/'embedded_dynamic_hht.i').write_text(text);(out/'embedded_three_axis.i').unlink()
    M=np.zeros((42,42));K=M.copy();mh,kh=hex_matrices(1e6,0.,1.);M[:24,:24]=mh;K[:24,:24]=kh
    for a in range(3):
        n1,n2=8+2*a,9+2*a;L=.5
        for j in range(3):
            ids=[3*n1+j,3*n2+j];M[np.ix_(ids,ids)]+=.00785*L/6*np.array([[2,1],[1,2]])
            if j==a:K[np.ix_(ids,ids)]+=200e9*1e-6/L*np.array([[1,-1],[-1,1]])
    C=np.zeros((18,42));fixed=[];load=np.zeros(42)
    for a in range(3):
        for n in range(8):
            if xyz[n,a]==0:fixed.append(3*n+a)
            if xyz[n,a]==1:load[3*n+a]=[1100,825,550][a]/4
    for n in range(8,14):
        N=np.prod(np.where(xyz[:8]==1,xyz[n],1-xyz[n]),axis=1)
        for j in range(3):C[(n-8)*3+j,3*n+j]=1;C[(n-8)*3+j,j:24:3]=-N
    np.savez(out/'independent-matrices.npz',M=M,K=K,C=C,coordinates=xyz,fixed=np.array(fixed),load_vector=load)
    (out/'fixture-identity.json').write_text(json.dumps(dict(alpha=ALPHA,beta=BETA,gamma=GAMMA,host_density_kg_m3=1,truss_density_kg_m3=7850,truss_area_m2=1e-6,masses_kg=dict(host=1,truss_total=.011775),reference='Independent host Bbar HEX8 and initial-axis 3 truss K/M plus 18 interpolation rows; all42 physical DOF and18 multipliers',limitations=['Host density1 is a diagnostic fixture, not original aqueduct material.', 'Small elastic axial bars; not multi-host, oblique bars, plastic evolution or full aqueduct.']),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();build(a.output)
