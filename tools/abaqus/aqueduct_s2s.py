#!/usr/bin/env python3
"""Conservative primal mortar Tie on planar QUAD4 source surfaces.

The weak equations integrate N_secondary*(u_secondary-u_primary) on the exact
polygon overlap. No nearest-node/nearest-face interpolation is used for weights.
Only faces whose vertices are all in the actual DAT-bound set participate. This
is a custom S2S discretization, not a reconstruction of proprietary Abaqus rows.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from map_aqueduct_interfaces import read_geometry, quad_shape, sha


def cross(a,b):return float(a[0]*b[1]-a[1]*b[0])

def area(p):return sum(cross(a,b) for a,b in zip(p,np.roll(p,-1,axis=0)))/2

def ccw(p):return p if area(p)>0 else p[::-1]

def clip(subject,boundary):
    result=list(subject)
    for a,b in zip(boundary,np.roll(boundary,-1,axis=0)):
        old=result;result=[]
        if not old:break
        for p,q in zip(old,old[1:]+old[:1]):
            fp=cross(b-a,p-a);fq=cross(b-a,q-a);ip=fp>=-1e-14;iq=fq>=-1e-14
            if ip:result.append(p)
            if ip!=iq:result.append(p+(q-p)*fp/(fp-fq))
    return np.array(result)

def shape_at(point,quad):
    natural=np.zeros(2)
    for _ in range(20):
        N,D=quad_shape(natural);delta=np.linalg.solve(quad.T@D,point-N@quad);natural+=delta
        if np.linalg.norm(delta,np.inf)<1e-13:break
    else:raise ValueError('inverse QUAD4 map failed')
    if np.max(np.abs(natural))>1+1e-7:raise ValueError('overlap quadrature outside original facet')
    return quad_shape(natural)[0]

def triangle_quadrature(p,order):
    g,w=leggauss(order);g=(g+1)/2;w=w/2
    for i in range(1,len(p)-1):
        a,b,c=p[0],p[i],p[i+1];jac=abs(cross(b-a,c-a))
        for r,wr in zip(g,w):
            for s,ws in zip(g,w):yield a+r*(b-a)+(1-r)*s*(c-a),jac*(1-r)*wr*ws

def mortar(coords,secondary,primary,bound,order=5):
    coords=np.asarray(coords,float);bound=set(bound)
    if any(len(f['nodes'])!=4 for f in secondary+primary):raise ValueError('S2S requires QUAD4 faces')
    pn=sorted({n for f in primary for n in f['nodes']});origin=coords[pn[0]-1]
    _,_,V=np.linalg.svd(coords[np.array(pn)-1]-origin,full_matrices=False);basis=V[:2].T;normal=V[2]
    error=float(np.max(np.abs((coords[np.array(pn)-1]-origin)@normal)))
    if error>1e-8:raise ValueError('primary surface is not planar: '+str(error))
    xy=(coords-origin)@basis
    sf=[f for f in secondary if set(f['nodes'])<=bound]
    if {n for f in sf for n in f['nodes']}!=bound:raise ValueError('bound node has no fully bound secondary face')
    D={n:{} for n in bound};B={n:{} for n in bound};slave_area=overlap_area=0.;pairs=0
    master=[(f,xy[np.array(f['nodes'])-1]) for f in primary]
    for f in sf:
        sn=f['nodes'];sq=xy[np.array(sn)-1];slave_area+=abs(area(sq))
        for mf,mq in master:
            if np.any(sq.max(0)<mq.min(0)-1e-13) or np.any(mq.max(0)<sq.min(0)-1e-13):continue
            poly=clip(ccw(sq),ccw(mq))
            if len(poly)<3 or abs(area(poly))<1e-18:continue
            overlap_area+=abs(area(poly));pairs+=1
            for point,weight in triangle_quadrature(poly,order):
                
                try:Ns=shape_at(point,sq);Nm=shape_at(point,mq)
                except ValueError as exc:raise ValueError(str((f,mf,sq.tolist(),mq.tolist(),point.tolist()))) from exc
                for i,n in enumerate(sn):
                    for j,k in enumerate(sn):D[n][k]=D[n].get(k,0.)+weight*Ns[i]*Ns[j]
                    for j,k in enumerate(mf['nodes']):B[n][k]=B[n].get(k,0.)+weight*Ns[i]*Nm[j]
    rows=[];max_constant=max_affine=0.
    for n in sorted(bound):
        lump=sum(D[n].values())
        if lump<=1e-16:raise ValueError('bound node has no integrated overlap')
        coeff={k:v/lump for k,v in D[n].items()}
        for k,v in B[n].items():coeff[k]=coeff.get(k,0.)-v/lump
        coeff={k:v for k,v in sorted(coeff.items()) if abs(v)>1e-14}
        max_constant=max(max_constant,abs(sum(coeff.values())))
        max_affine=max(max_affine,float(np.linalg.norm(sum((v*xy[k-1] for k,v in coeff.items()),np.zeros(2)),np.inf)))
        rows.append(dict(anchor=n-1,scalar_terms=[dict(node=k-1,coefficient=v) for k,v in coeff.items()],integrated_test_area_m2=lump,primary_nodes=sorted(B[n]),secondary_nodes=sorted(D[n])))
    if max_constant>1e-10 or max_affine>1e-10:raise ValueError('mortar affine/constant patch failed')
    return rows,dict(bound_nodes=len(bound),fully_bound_faces=len(sf),excluded_faces=len(secondary)-len(sf),intersection_pairs=pairs,
                     secondary_area_m2=slave_area,overlap_area_m2=overlap_area,coverage_fraction=overlap_area/slave_area,
                     max_constant_residual=max_constant,max_affine_residual_m=max_affine,primary_planarity_error_m=error,
                     plane_origin=origin.tolist(),plane_normal=normal.tolist(),quadrature_order=order)

def build(preparation,membership,output):
    root=Path(preparation);out=Path(output)
    if out.exists():raise ValueError('output exists')
    contract=json.loads((root/'assembly-contract.json').read_text())
    for name,digest in contract['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('source preparation hash mismatch: '+name)
    members=json.loads(Path(membership).read_text());coords,*_=read_geometry(root);coords=np.asarray(coords)
    rows=[];summaries=[];adjustments=[]
    for tie in contract['interfaces']['ties']:
        name=tie['options']['name'];group=[r for r in members if r['tie'].lower()==name.lower()];bound={r['exodus_node_1based'] for r in group if r['bound_in_printed_table']}
        s,m=[contract['surfaces'][n] for n in tie['surfaces']]
        if {r['exodus_node_1based'] for r in group}!={n for f in s for n in f['nodes']}:raise ValueError('DAT candidate membership differs from source surface')
        scalar,summary=mortar(coords,s,m,bound)
        origin=np.array(summary['plane_origin']);normal=np.array(summary['plane_normal'])
        for n in sorted(bound):
            point=coords[n-1].copy();distance=float((point-origin)@normal);coords[n-1]-=distance*normal
            adjustments.append(dict(node=n,original=point.tolist(),adjusted=coords[n-1].tolist(),distance_m=abs(distance),provenance='candidate-normal-plane-adjustment-not-ODB-confirmed'))
        support_differences=0
        for row in scalar:
            actual=next(r['main_exodus_nodes_1based'] for r in group if r['exodus_node_1based']==row['anchor']+1)
            support_differences+=set(actual)!=set(row['primary_nodes'])
            for j,v in enumerate(('ux','uy','uz')):
                rows.append(dict(kind='tie_s2s',tie=name,source_line=tie['source_line'],anchor=row['anchor'],multiplier=j,
                                 terms=[dict(variable=v,**t) for t in row['scalar_terms']]))
        summaries.append(dict(name=name,unbound_nodes=len(group)-len(bound),Abaqus_primary_support_difference_rows=support_differences,**summary))
    out.mkdir(parents=True)
    record=dict(schema_version='aq-s2s-primal-mortar-v1',operator='custom consistent primal mortar; not Abaqus coefficient reconstruction',source_contract_sha256=sha(root/'assembly-contract.json'),membership_sha256=sha(membership),ties=summaries,rows=rows,adjustments=adjustments,
                limitations=['Only all-bound secondary faces participate; partial-bound faces are explicitly excluded.', 'Abaqus finite-area averaging coefficients are unavailable. Consistent mortar couples several secondary nodes per weak row; stencil and discretization differ.', 'Normal plane adjustments are geometry-derived candidates. Exact Abaqus adjusted coordinates await existing ODB read-only export.'])
    (out/'s2s-operator.json').write_text(json.dumps(record,indent=2)+'\n')
    np.save(out/'adjusted-coordinates.npy',coords)
    return record

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',required=True);p.add_argument('--membership',required=True);p.add_argument('--output',required=True);a=p.parse_args();d=build(a.preparation,a.membership,a.output);print(json.dumps(dict(rows=len(d['rows']),ties=d['ties']),indent=2))
