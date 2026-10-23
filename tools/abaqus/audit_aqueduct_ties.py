#!/usr/bin/env python3
"""Measure source Tie geometry. Distances are NOT Abaqus tied-node admission."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from map_aqueduct_interfaces import project_quad, read_geometry, sha


def audit(preparation, output):
    root=Path(preparation);out=Path(output)
    if out.exists():raise ValueError('output exists')
    contract=json.loads((root/'assembly-contract.json').read_text())
    for name,digest in contract['artifacts'].items():
        if sha(root/name)!=digest:raise ValueError('source artifact identity differs')
    coords,*_=read_geometry(root);records=[];summaries=[]
    for tie in contract['interfaces']['ties']:
        secondary,primary=tie['surfaces'];faces=contract['surfaces'][primary]
        xyz=np.array([coords[np.array(f['nodes'])-1] for f in faces])
        centers=xyz.mean(1);radius=np.linalg.norm(xyz-centers[:,None,:],axis=2).max();tree=cKDTree(centers)
        diagonals=np.maximum(np.linalg.norm(xyz[:,0]-xyz[:,2],axis=1),np.linalg.norm(xyz[:,1]-xyz[:,3],axis=1))
        nodes=sorted({n for f in contract['surfaces'][secondary] for n in f['nodes']});distances=[]
        for n in nodes:
            point=coords[n-1];_,seed=tree.query(point);bound=project_quad(point,xyz[seed])['distance']
            # Triangle inequality excludes all facets whose closest possible point is farther.
            candidates=tree.query_ball_point(point,radius+bound+1e-12)
            distance,index,p=min((p['distance'],i,p) for i in candidates for p in [project_quad(point,xyz[i])])
            assert abs(sum(p['weights'])-1)<1e-12 and min(p['weights'])>=-1e-12
            assert np.linalg.norm(p['point']-p['weights']@xyz[index])<1e-12
            distances.append(distance);records.append([tie['options']['name'],n,distance,*point,*p['point'],faces[index]['instance'],faces[index]['element_label'],faces[index]['face']])
        summaries.append(dict(name=tie['options']['name'],source_line=tie['source_line'],secondary_nodes=len(nodes),max_closest_distance_m=max(distances),main_facet_diagonal_min_m=float(min(diagonals)),main_facet_diagonal_median_m=float(np.median(diagonals)),main_facet_diagonal_max_m=float(max(diagonals)),counts_below_distance_m={str(v):sum(d<=v for d in distances) for v in (1e-5,.005,.01,.02,.03,.05,.1,.2)},abaqus_tied_nodes=None,abaqus_untied_nodes=None))
    out.mkdir(parents=True)
    with (out/'tie-closest-geometry.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['tie_name','exodus_node_1_based','closest_distance_m','x_m','y_m','z_m','projected_x_m','projected_y_m','projected_z_m','primary_instance','primary_element_label','primary_face']);w.writerows(records)
    result=dict(status='geometry-only-not-abaqus-tie-admission',source_contract_sha256=sha(root/'assembly-contract.json'),source_sha256=contract['source']['sha256'],ties=summaries,limitations=['No nodes bound or moved.','Abaqus Standard default S2S averages over a finite region and may include adjacent secondary face nodes.','Typical diagonal and significant secondary-face coverage are not reconstructed by thresholding this diagnostic.','Current nodal projection mapper is not asserted to reproduce S2S coefficients.'],artifacts={'tie-closest-geometry.csv':sha(out/'tie-closest-geometry.csv')})
    (out/'tie-geometry-diagnostic.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps(audit(a.preparation,a.output),ensure_ascii=False,indent=2))
