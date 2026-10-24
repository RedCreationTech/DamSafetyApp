#!/usr/bin/env python3
"""Free dynamic nonmatching solid/solid mortar with an independent global matrix."""
import argparse,json
from pathlib import Path
import numpy as np
from abaqus2exodus import GlobalMesh,write_exodus
from aqueduct_s2s import mortar
from build_aqueduct_interface_cases import block
from build_aqueduct_full_input import quoted,ALPHA,BETA,GAMMA
from build_aqueduct_internal_fixtures import hex_matrices

def build(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False);gm=GlobalMesh(1e-10,merge_nodes=False)
    xyz=np.array([[x,y,z] for z in (0.,1.) for y in (0.,.5,1.) for x in (0.,.5,1.)]+[[x,y,z] for z in (1.,2.) for x,y in [(0,0),(1,0),(1,1),(0,1)]])
    for n,p in enumerate(xyz,1):gm.add(('fixture',n),tuple(p))
    lower=[(i+1,[a,a+1,a+4,a+3,a+9,a+10,a+13,a+12]) for i,a in enumerate((1,2,4,5))];upper=[(5,list(range(19,27)))]
    write_exodus(out/'fixture.e',gm,{'lower':lower,'upper':upper},{'lower':'HEX8','upper':'HEX8'},{},{'base':list(range(1,10)),'load':list(range(23,27))},{},'Nonmatching 9-to-4 node primal S2S mortar free dynamics')
    primary=[dict(nodes=e[1][4:]) for e in lower];secondary=[dict(nodes=list(range(19,23)))];scalar,meta=mortar(xyz,secondary,primary,list(range(19,23)))
    rows=[]
    for r in scalar:
        for j,v in enumerate(('ux','uy','uz')):rows.append(dict(anchor=r['anchor'],component=j,terms=[dict(variable=v,**t) for t in r['scalar_terms']]))
    text=block('Mesh',block('source','type = FileMeshGenerator\nfile = fixture.e')+block('support',"type = AqueductConstraintNodesGenerator\ninput = source\nnodes = '18 19 20 21'\nblock_id = 10\nblock_name = support")+'allow_renumbering = false\nparallel_type = replicated')+block('Problem','material_coverage_check = false')+block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Variables',''.join(block(v,"block = '1 2'") for v in ('ux','uy','uz'))+''.join(block('lm'+str(j),'block = support') for j in range(3)))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('solid',"block = '1 2'\nstrain = SMALL\nvolumetric_locking_correction = true"))
    text+=block('Materials',block('elastic',"type = ComputeIsotropicElasticityTensor\nblock = '1 2'\nyoungs_modulus = 1e6\npoissons_ratio = 0")+block('stress',"type = ComputeLinearElasticStress\nblock = '1 2'")+block('mass',"type = GenericConstantMaterial\nblock = '1 2'\nprop_names = density\nprop_values = 2000"))
    text+=block('Kernels',''.join(block('lm'+str(j),f'type = NullKernel\nvariable = lm{j}\nblock = support\ndisplacements = \'\'\njacobian_fill = 0') for j in range(3))+''.join(block('mass_'+v,f'type = AqueductHHTInertia\nvariable = {v}\nblock = \'1 2\'\nhht_alpha = {ALPHA}\nhht_start_time = 0\nuse_displaced_mesh = false') for v in ('ux','uy','uz')))
    text+=block('Constraints',''.join(block('row'+str(i),f'type = AqueductLinearMPC\nvariable = ux\nanchor_node = {r["anchor"]}\nterm_nodes = {quoted(t["node"] for t in r["terms"])}\nterms = {quoted(t["variable"] for t in r["terms"])}\ncoefficients = {quoted(t["coefficient"] for t in r["terms"])}\nmultiplier = lm{r["component"]}') for i,r in enumerate(rows)))
    text+=block('BCs',''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = base\nvalue = 0') for v in ('ux','uy','uz')))
    text+=block('Functions',block('load',"type = ParsedFunction\nexpression = '1000*t'"))
    text+=block('NodalKernels',block('load','type = UserForcingFunctorNodalKernel\nvariable = uz\nboundary = load\nfunctor = load'))
    post=''
    for n in range(26):
        for v in ('ux','uy','uz'):post+=block(f'n{n}_{v}',f'type = NodalVariableValue\nvariable = {v}\nnodeid = {n}')
    for r in rows:post+=block(f'lm_{r["anchor"]}_{r["component"]}',f'type = NodalVariableValue\nvariable = lm{r["component"]}\nnodeid = {r["anchor"]}')
    text+=block('Postprocessors',post)+block('Preconditioning',block('smp','type = SMP\nfull = true'))+block('Executioner',"type = Transient\nsolve_type = NEWTON\nend_time = .03\ndt = .01\nnl_abs_tol = 1e-9\nnl_rel_tol = 1e-11\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'\n"+block('TimeIntegrator',f'type = AqueductHHT\nhht_alpha = {ALPHA}\nbeta = {BETA}\ngamma = {GAMMA}\nalgebraic_variables = \'lm0 lm1 lm2\''))+block('Outputs','csv = true\nexodus = true')
    (out/'s2s_hht.i').write_text(text)
    # Independent quadrature on physical subcells, no MOOSE assembly values used.
    nd=78;M=np.zeros((nd,nd));K=M.copy();sign=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
    Cmat=np.diag([1e6,1e6,1e6,5e5,5e5,5e5])
    for _,nodes in lower+upper:
        ids=np.array([3*(n-1)+j for n in nodes for j in range(3)]);ex=xyz[np.array(nodes)-1]
        for q in np.array(np.meshgrid(*[[-1/np.sqrt(3),1/np.sqrt(3)]]*3)).T.reshape(-1,3):
            N=np.prod(1+sign*q,axis=1)/8;D=np.empty((8,3))
            for j in range(3):D[:,j]=sign[:,j]*np.prod(1+sign[:,[k for k in range(3) if k!=j]]*q[[k for k in range(3) if k!=j]],axis=1)/8
            J=ex.T@D;grad=D@np.linalg.inv(J);mean=sign/8@np.linalg.inv(J);B=np.zeros((6,24))
            for i,g in enumerate(grad):
                B[0,3*i]=g[0];B[1,3*i+1]=g[1];B[2,3*i+2]=g[2]
                B[3,3*i:3*i+3]=[g[1],g[0],0];B[4,3*i:3*i+3]=[0,g[2],g[1]];B[5,3*i:3*i+3]=[g[2],0,g[0]]
                for j in range(3):B[:3,3*i+j]+=(mean[i,j]-g[j])/3
            K[np.ix_(ids,ids)]+=B.T@Cmat@B*np.linalg.det(J)
            for j in range(3):M[np.ix_(ids[j::3],ids[j::3])]+=2000*np.outer(N,N)*np.linalg.det(J)
    C=np.zeros((12,nd))
    for i,r in enumerate(rows):
        for t in r['terms']:C[i,3*t['node']+r['component']]+=t['coefficient']
    np.savez(out/'independent-matrices.npz',M=M,K=K,C=C,coordinates=xyz,fixed=np.arange(27),load_dofs=np.arange(22*3+2,26*3,3))
    (out/'fixture-identity.json').write_text(json.dumps(dict(mortar=meta,rows=rows,alpha=ALPHA,beta=BETA,gamma=GAMMA,limitations=['Independent global solve uses the frozen custom C; area quadrature and affine/virtual-work properties are separately tested.', 'Not Abaqus proprietary S2S coefficients or full aqueduct.']),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();build(a.output)
