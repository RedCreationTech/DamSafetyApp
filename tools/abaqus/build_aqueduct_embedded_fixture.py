#!/usr/bin/env python3
"""A free-solve HEX8/interior-truss force-transfer check; not the aqueduct."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from abaqus2exodus import GlobalMesh, write_exodus
from build_aqueduct_interface_cases import block


def build(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    gm=GlobalMesh(1e-10,merge_nodes=False)
    coordinates=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1),(.25,.5,.5),(.75,.5,.5)]
    for i,point in enumerate(coordinates,1):gm.add(('fixture',i),point)
    write_exodus(out/'embedded_fixture.e',gm,{'host':[(1,list(range(1,9)))],'steel':[(2,[9,10])]},
                 {'host':'HEX8','steel':'TRUSS'},{},
                 {'left':[1,4,5,8],'right':[2,3,6,7],'all':list(range(1,11))}, {},
                 'Independent interior embedded steel axial force transfer')
    with (out/'embedded-map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['secondary_node','primary_element','secondary_block','x_m','y_m','z_m']+[f'host_node_{j}' for j in range(8)])
        for n in (9,10):w.writerow([n-1,0,2,*coordinates[n-1],*range(8)])
    text=block('Mesh','type = FileMesh\nfile = embedded_fixture.e\nallow_renumbering = false')
    text+=block('GlobalParams',"displacements = 'ux uy uz'")
    text+=block('Variables',''.join(block(v,'') for v in ('ux','uy','uz')))
    text+=block('Physics/SolidMechanics/QuasiStatic',block('host',"block = host\nstrain = SMALL\nsave_in = 'host_rx host_ry host_rz'"))
    text+=block('AuxVariables',''.join(block(v,'') for v in ('host_rx','host_ry','host_rz'))+
                block('area','block = steel\nfamily = MONOMIAL\norder = CONSTANT\ninitial_condition = 1e-6')+
                block('axial_stress','block = steel\nfamily = MONOMIAL\norder = CONSTANT'))
    text+=block('Materials',block('elastic','type = ComputeIsotropicElasticityTensor\nblock = host\nyoungs_modulus = 1e6\npoissons_ratio = 0')+
                block('stress','type = ComputeLinearElasticStress\nblock = host')+
                block('steel',"type = AqueductSteelTruss\nblock = steel\ndisplacements = 'ux uy uz'\nyoungs_modulus = 200e9\nyield_stress = 400e6\nultimate_stress = 696e6\nplastic_strain_at_ultimate = .1"))
    text+=block('Kernels',''.join(block(v,f"type = StressDivergenceTensorsTruss\nvariable = {v}\ncomponent = {j}\nblock = steel\ndisplacements = 'ux uy uz'\narea = area\nuse_displaced_mesh = false") for j,v in enumerate(('ux','uy','uz'))))
    text+=block('AuxKernels',block('stress','type = MaterialRealAux\nvariable = axial_stress\nproperty = axial_stress\nblock = steel'))
    text+=block('Constraints',block('embedded','type = AqueductEmbeddedConstraint\nvariable = ux\nprimary_variable = ux\nprimary = host\nsecondary = steel\nallowed_primary_block_ids = 1\nsource_map = embedded-map.csv\nformulation = kinematic\npenalty = 1e10'))
    text+=block('BCs',block('left','type = DirichletBC\nvariable = ux\nboundary = left\nvalue = 0')+
                ''.join(block(v,f'type = DirichletBC\nvariable = {v}\nboundary = all\nvalue = 0') for v in ('uy','uz')))
    text+=block('NodalKernels',block('load','type = UserForcingFunctorNodalKernel\nvariable = ux\nboundary = right\nfunctor = 275'))
    post=''.join(block(name,f'type = NodalVariableValue\nvariable = ux\nnodeid = {n}') for name,n in [('right_u',1),('bar_left_u',8),('bar_right_u',9)])
    post+=block('steel_stress','type = ElementAverageValue\nvariable = axial_stress\nblock = steel')
    post+=block('host_right_force','type = NodalSum\nvariable = host_rx\nboundary = right')
    text+=block('Postprocessors',post)
    text+=block('Preconditioning',block('smp','type = SMP\nfull = true'))
    text+=block('Executioner',"type = Steady\nsolve_type = NEWTON\nnl_abs_tol = 1e-7\nnl_rel_tol = 1e-12\npetsc_options_iname = '-pc_type'\npetsc_options_value = 'lu'")
    text+=block('Outputs','csv = true\nexodus = true')
    (out/'embedded_transfer.i').write_text(text)
    # Independent energy: K=Eh Ah/Lh + Es As Ls/Lh^2 = 1.1e6 N/m.
    expected=dict(right_u=.001,bar_left_u=.00025,bar_right_u=.00075,steel_stress=200e6,host_right_force=1000.)
    identity=dict(scope='independent axial host/interior-truss fixture',expected=expected,
                  independent_equation='K=Eh*Ah/Lh+Es*As*Ls/Lh^2; F=1100 N; bar force=200 N; weighted host transfer=100 N',
                  limitations=['No Abaqus reference or CDP material.','Only axial embedded DOF is constrained; transverse DOFs prescribed zero.','Not a Tie or shell-solid test; not a full aqueduct input.'],
                  artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()})
    (out/'fixture-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    return identity


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);print(json.dumps(build(p.parse_args().output),indent=2))
