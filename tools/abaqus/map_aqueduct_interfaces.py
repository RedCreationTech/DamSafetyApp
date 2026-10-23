#!/usr/bin/env python3
"""Build explicit source-restricted interface rows; no runnable full input or solve."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import netCDF4
import numpy as np
from scipy.spatial import cKDTree

from audit_aqueduct_geometry import locate_embedded
from build_aqueduct_interface_cases import block


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def quad_shape(x):
    sx=np.array([-1,1,1,-1]);sy=np.array([-1,-1,1,1]);r,s=x
    return (1+sx*r)*(1+sy*s)/4,np.column_stack([sx*(1+sy*s),sy*(1+sx*r)])/4


def project_quad(point,quad):
    x=np.zeros(2)
    for _ in range(30):
        N,D=quad_shape(x);J=quad.T@D
        score=J.T@(point-N@quad)
        free=[j for j in range(2) if not ((x[j]>=1-1e-13 and score[j]>0) or (x[j]<=-1+1e-13 and score[j]<0))]
        if not free:break
        delta=np.zeros(2);delta[free]=np.linalg.lstsq(J[:,free],point-N@quad,rcond=None)[0]
        new=np.clip(x+delta,-1,1)
        if np.max(np.abs(new-x))<1e-13:x=new;break
        x=new
    N,D=quad_shape(x);J=quad.T@D
    gradient=D@np.linalg.inv(J.T@J)@J.T
    normal=np.cross(J[:,0],J[:,1]);normal/=np.linalg.norm(normal)
    return dict(point=N@quad,weights=N,gradients=gradient,normal=normal,natural=x,
                distance=float(np.linalg.norm(point-N@quad)))


def project_surface(points,facets,coords,max_distance=1e-5):
    if any(len(f['nodes'])!=4 for f in facets):raise ValueError('primary surface must contain only quadrilateral solid faces')
    xyz=np.array([coords[np.array(f['nodes'])-1] for f in facets]);centers=xyz.mean(1)
    radius=float(np.linalg.norm(xyz-centers[:,None,:],axis=2).max());tree=cKDTree(centers)
    results=[]
    for point in points:
        candidates=[]
        for j in tree.query_ball_point(point,radius+max_distance):
            p=project_quad(point,xyz[j])
            candidates.append((p['distance'],facets[j]['instance'],facets[j]['element_label'],facets[j]['face'],j,p))
        if not candidates:raise ValueError('no source primary facet within the bounded search')
        best=min(candidates,key=lambda x:x[:5]);p=best[-1]
        if p['distance']>max_distance:raise ValueError(f'source interface projection exceeds analysis cap: {p["distance"]}')
        p['facet']=facets[best[-2]];results.append(p)
    return results


def rotation_operator(gradients,normal):
    return np.array([np.outer(np.cross(g,normal),normal)+.5*np.outer(normal,np.cross(normal,g)) for g in gradients])


def cross_matrix(r):return np.array([[0,-r[2],r[1]],[r[2],0,-r[0]],[-r[1],r[0],0]])


def read_geometry(root):
    with netCDF4.Dataset(root/'assembly-geometry.e') as nc:
        coords=np.column_stack([nc[k][:] for k in ('coordx','coordy','coordz')])
        ids=np.array(nc['elem_num_map'][:],int);blocks={};connectivity={};types={};offset=0
        for j,block in enumerate(np.array(nc['eb_prop1'][:],int),1):
            conn=np.array(nc[f'connect{j}'][:],int);typ=nc[f'connect{j}'].elem_type
            for k,row in enumerate(conn):
                eid=int(ids[offset+k]);connectivity[eid]=row;blocks[eid]=int(block);types[eid]=typ
            offset+=len(conn)
    with (root/'element-map.csv').open() as f:
        origins={(r['instance'],int(r['abaqus_element_label'])):int(r['exodus_element']) for r in csv.DictReader(f)}
    return coords,connectivity,blocks,types,origins


def build_embedded_precheck(preparation,audit,output):
    """Validate original-geometry host identities without choosing any Tie policy."""
    root=Path(preparation);audit=Path(audit);out=Path(output)
    if out.exists():raise ValueError('output already exists')
    contract=json.loads((root/'assembly-contract.json').read_text())
    for name,digest in contract['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('source preparation hash mismatch: '+name)
    summary=json.loads((audit/'geometry-audit.json').read_text())
    if summary['contract_sha256']!=sha(root/'assembly-contract.json'):raise ValueError('audit identity differs')
    coords,connectivity,blocks,types,origins=read_geometry(root);rows=[];groups=[]
    for j,g in enumerate(contract['interfaces']['embedded'],1):
        detail=audit/f'embedded-{j:02d}.json'
        if sha(detail)!=summary['embedded'][j-1]['sha256']:raise ValueError('source host audit hash mismatch')
        owner={}
        for origin in g['embedded']:
            eid=origins[tuple(origin)]
            if types[eid]!='TRUSS':raise ValueError('non-truss embedded child')
            for n in connectivity[eid]:owner[int(n)]=min(owner.get(int(n),blocks[eid]),blocks[eid])
        d=json.loads(detail.read_text())
        if d['unresolved'] or {r['node'] for r in d['found']}!=set(owner):raise ValueError('incomplete audited node coverage')
        for r in d['found']:
            n=r['node'];eid=origins[tuple(r['host'])]
            rows.append([n-1,eid-1,owner[n],*coords[n-1],*[int(v)-1 for v in connectivity[eid]]])
        groups.append(dict(primary_blocks=sorted({blocks[origins[tuple(x)]] for x in g['host']}),owning_secondary_blocks=sorted(set(owner.values()))))
    if len({r[0] for r in rows})!=len(rows):raise ValueError('duplicate embedded node owner')
    out.mkdir(parents=True);shutil.copyfile(root/'assembly-geometry.e',out/'assembly-geometry.e')
    with (out/'embedded-runtime-map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['secondary_node','primary_element','secondary_block','x_m','y_m','z_m']+[f'host_node_{j}' for j in range(8)]);w.writerows(rows)
    text=block('Mesh','type = FileMesh\nfile = assembly-geometry.e\nallow_renumbering = false')
    text+=block('Variables',''.join(block(v,'') for v in ('ux','uy','uz')))
    text+=block('Kernels',''.join(block(v,f'type = NullKernel\nvariable = {v}') for v in ('ux','uy','uz')))
    constraints=''
    for j,g in enumerate(groups):
        for owner in g['owning_secondary_blocks']:
            for v in ('ux','uy','uz'):
                constraints+=block(f'embedded_{j}_{owner}_{v}',f'type = AqueductEmbeddedConstraint\nvariable = {v}\nprimary_variable = {v}\nprimary = {g["primary_blocks"][0]}\nsecondary = {owner}\nallowed_primary_block_ids = \'{" ".join(map(str,g["primary_blocks"]))}\'\nsource_map = embedded-runtime-map.csv\nformulation = kinematic\npenalty = 1e10')
    text+=block('Constraints',constraints)
    text+=block('UserObjects',block('no_solve',"type = Terminator\nexpression = '1'\nexecute_on = INITIAL\nfail_mode = HARD\nerror_level = ERROR\nmessage = 'Original geometry host audit only. No physical solve defined. Use --check-input.'"))
    text+=block('Executioner','type = Steady')
    (out/'embedding_check_only.i').write_text(text)
    record=dict(status='original-geometry-host-precheck-only',embedded_nodes=len(rows),Tie_policy_applied=False,physical_solve=False,source_contract_sha256=sha(root/'assembly-contract.json'),artifacts={p.name:sha(p) for p in out.iterdir()})
    (out/'precheck-identity.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def build(preparation,audit,output):
    root=Path(preparation);audit=Path(audit);out=Path(output)
    if out.exists():raise ValueError('output already exists')
    contract=json.loads((root/'assembly-contract.json').read_text())
    for name,digest in contract['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('source preparation hash mismatch: '+name)
    audit_summary=json.loads((audit/'geometry-audit.json').read_text())
    if audit_summary['contract_sha256']!=sha(root/'assembly-contract.json'):raise ValueError('audit identity differs')
    coords,connectivity,blocks,types,origins=read_geometry(root)
    rows=[];adjusted=[];tied_nodes=set();faces=contract['surfaces']
    for tie in contract['interfaces']['ties']:
        slave,primary=tie['surfaces'];nodes=sorted({n for f in faces[slave] for n in f['nodes']})
        projected=project_surface([coords[n-1] for n in nodes],faces[primary],coords)
        for node,p in zip(nodes,projected):
            if node in tied_nodes:raise ValueError('multiple Tie definitions own one slave node')
            tied_nodes.add(node)
            if tie['options'].get('adjust','').lower()!='yes':raise ValueError('source adjust semantics changed')
            adjusted.append(dict(node=node,original=coords[node-1].tolist(),adjusted=p['point'].tolist(),distance_m=p['distance']))
            coords[node-1]=p['point']
            for j in range(3):
                terms=[dict(node=node-1,variable=('ux','uy','uz')[j],coefficient=1.)]
                terms += [dict(node=n-1,variable=('ux','uy','uz')[j],coefficient=-float(w)) for n,w in zip(p['facet']['nodes'],p['weights']) if abs(w)>1e-15]
                rows.append(dict(kind='tie',source_line=tie['source_line'],anchor=node-1,multiplier=j,terms=terms))
    shell_distances=[]
    for coupling in contract['interfaces']['shell_solid']:
        slave,primary=coupling['surfaces'];nodes=sorted({n for f in faces[slave] for n in f['nodes']})
        for node,p in zip(nodes,project_surface([coords[n-1] for n in nodes],faces[primary],coords)):
            shell_distances.append(p['distance']);A=rotation_operator(p['gradients'],p['normal']);r=coords[node-1]-p['point']
            P=np.array([N*np.eye(3)-cross_matrix(r)@a for N,a in zip(p['weights'],A)])
            for component in range(6):
                terms=[dict(node=node-1,variable=('ux','uy','uz','rx','ry','rz')[component],coefficient=1.)]
                operator=P[:,component,:] if component<3 else A[:,component-3,:]
                for n,coeff in zip(p['facet']['nodes'],operator):
                    for j,c in enumerate(coeff):
                        if abs(c)>1e-15:terms.append(dict(node=n-1,variable=('ux','uy','uz')[j],coefficient=-float(c)))
                rows.append(dict(kind='shell_solid',source_line=coupling['source_line'],anchor=node-1,multiplier=component,terms=terms))
    for connector in contract['interfaces']['connectors']:
        a,b=[x['node']-1 for x in connector['endpoints']]
        for j,rotation in enumerate(('rx','ry','rz')):
            rows.append(dict(kind='align',source_line=connector['source_line'],anchor=a,multiplier=j+3,terms=[dict(node=a,variable=rotation,coefficient=1.),dict(node=b,variable=rotation,coefficient=-1.)]))
    if len({(r['anchor'],r['multiplier']) for r in rows})!=len(rows):raise ValueError('duplicate multiplier row ownership')
    embedded=[];embedded_groups=[]
    for j,g in enumerate(contract['interfaces']['embedded'],1):
        source_details=audit/f'embedded-{j:02d}.json'
        if sha(source_details)!=audit_summary['embedded'][j-1]['sha256']:raise ValueError('source host audit hash mismatch')
        host_origins=[tuple(x) for x in g['host']];host_ids=[origins[x] for x in host_origins]
        child_ids=[origins[tuple(x)] for x in g['embedded']]
        owner={}
        for eid in child_ids:
            if types[eid]!='TRUSS':raise ValueError('non-truss embedded child')
            for n in connectivity[eid]:owner[int(n)]=min(owner.get(int(n),blocks[eid]),blocks[eid])
        found,missing=locate_embedded([(n,coords[n-1]) for n in sorted(owner)],np.array([coords[connectivity[e]-1] for e in host_ids]),host_origins)
        if missing:raise ValueError(f'Tie adjustment invalidated {len(missing)} embedded nodes')
        for p in found:
            eid=origins[tuple(p['host'])];n=p['node']
            embedded.append([n-1,eid-1,owner[n],*coords[n-1],*[int(x)-1 for x in connectivity[eid]]])
        embedded_groups.append(dict(source_line=g['source_line'],nodes=len(found),primary_blocks=sorted({blocks[e] for e in host_ids}),owning_secondary_blocks=sorted(set(owner.values()))))
    if len({x[0] for x in embedded})!=len(embedded):raise ValueError('one embedded node appears in multiple source groups')
    out.mkdir(parents=True)
    for name in ('assembly-geometry.e','element-map.csv','source-utf8.inp'):shutil.copyfile(root/name,out/name)
    with netCDF4.Dataset(out/'assembly-geometry.e','r+') as nc:
        for j,k in enumerate(('coordx','coordy','coordz')):nc[k][:]=coords[:,j]
        if 'num_side_sets' in nc.dimensions:raise ValueError('preparation sidesets already exist; do not replace them implicitly')
        nss=2*len(contract['interfaces']['contacts']);nc.createDimension('num_side_sets',nss)
        nc.createVariable('ss_status','i4',('num_side_sets',))[:]=1
        ids=nc.createVariable('ss_prop1','i4',('num_side_sets',));ids.setncattr('name','ID');ids[:]=np.arange(10001,10001+nss)
        names=[];order=np.array(nc['elem_num_map'][:],int);internal={eid:i+1 for i,eid in enumerate(order)}
        for j,contact in enumerate(contract['interfaces']['contacts'],1):
            for role,surf in zip(('secondary','primary'),contact['surfaces']):
                index=len(names)+1;names.append(f'aq_contact_{j:02d}_{role}')
                facets=faces[surf];nc.createDimension(f'num_side_ss{index}',len(facets))
                es=nc.createVariable(f'elem_ss{index}','i4',(f'num_side_ss{index}',));ss=nc.createVariable(f'side_ss{index}','i4',(f'num_side_ss{index}',))
                es[:]=[internal[origins[(f['instance'],f['element_label'])]] for f in facets]
                ss[:]=[{'S1':5,'S2':6,'S3':1,'S4':2,'S5':3,'S6':4}[f['face']] for f in facets]
        names_var=nc.createVariable('ss_names','S1',('num_side_sets','len_name'))
        width=len(nc.dimensions['len_name']);names_var[:]=netCDF4.stringtochar(np.array(names,dtype=f'S{width}'))
    with (root/'node-map.csv').open() as f:node_rows=list(csv.DictReader(f))
    with (out/'node-map.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(node_rows[0]));w.writeheader()
        for r in node_rows:
            for j,k in enumerate(('x_m','y_m','z_m')):r[k]=coords[int(r['exodus_node'])-1,j]
            w.writerow(r)
    with (out/'embedded-runtime-map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['secondary_node','primary_element','secondary_block','x_m','y_m','z_m']+[f'host_node_{i}' for i in range(8)]);w.writerows(embedded)
    check='# Geometry/host identity precheck only. INITIAL guard rejects any actual solve.\n'
    check+=block('Mesh',"type = FileMesh\nfile = assembly-geometry.e\nallow_renumbering = false")
    check+=block('GlobalParams',"displacements = 'ux uy uz'")
    check+=block('Variables',''.join(block(v,'') for v in ('ux','uy','uz')))
    check+=block('Kernels',''.join(block(v,f'type = NullKernel\nvariable = {v}') for v in ('ux','uy','uz')))
    constraints=''
    for j,g in enumerate(embedded_groups):
        for owner in g['owning_secondary_blocks']:
            for v in ('ux','uy','uz'):
                constraints+=block(f'embedded_{j}_{owner}_{v}',f'type = AqueductEmbeddedConstraint\nvariable = {v}\nprimary_variable = {v}\nprimary = {g["primary_blocks"][0]}\nsecondary = {owner}\nallowed_primary_block_ids = \'{" ".join(map(str,g["primary_blocks"]))}\'\nsource_map = embedded-runtime-map.csv\nformulation = kinematic\npenalty = 1e10')
    check+=block('Constraints',constraints)
    check+=block('UserObjects',block('audit_only',"type = Terminator\nexpression = '1'\nexecute_on = INITIAL\nfail_mode = HARD\nerror_level = ERROR\nmessage = 'Geometry audit only: use --check-input; no physical solve is defined.'"))
    check+=block('Executioner','type = Steady')
    (out/'embedding_check_only.i').write_text(check)
    implementation=dict(status='interface-rows-not-full-solve-input',runtime_ids='Exodus original maps minus one; must pass runtime identity checks',rows=rows,embedded_groups=embedded_groups,tie_adjustments=adjusted,projection_analysis_cap_m=1e-5,projection_cap_is_not_abaqus_default=True,shell_solid_max_projection_distance_m=max(shell_distances),shell_solid_operator='surface gradient rotation with offset translation; rigid-mode and virtual-work preserving custom operator, not asserted Abaqus-discretization-equivalent',pending=['source elastic-slip contact','full materials/time-history input','runtime coverage and combination mechanics','same-binary frozen CDP regressions'])
    (out/'interface-implementation.json').write_text(json.dumps(implementation,indent=2)+'\n')
    contract['derivation']=dict(source_preparation_contract_sha256=sha(root/'assembly-contract.json'),Tie_adjust_yes_applied=True,original_source_unchanged=True)
    contract['artifacts']={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'assembly-contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2)+'\n')
    summary=dict(status=implementation['status'],rows=len(rows),embedded_nodes=len(embedded),tie_adjusted_nodes=len(adjusted),max_tie_adjustment_m=max(x['distance_m'] for x in adjusted),shell_solid_max_projection_distance_m=implementation['shell_solid_max_projection_distance_m'],contact_sidesets=nss,contract_sha256=sha(out/'assembly-contract.json'))
    (out/'interface-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',required=True);p.add_argument('--audit',required=True);p.add_argument('--output',required=True);p.add_argument('--embedded-precheck-only',action='store_true');a=p.parse_args();print(json.dumps((build_embedded_precheck if a.embedded_precheck_only else build)(a.preparation,a.audit,a.output),indent=2))
