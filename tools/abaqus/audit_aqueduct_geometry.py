#!/usr/bin/env python3
"""Audit prepared geometry and host-restricted embedded interpolation, no FE solve."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import netCDF4
import numpy as np
from scipy.spatial import cKDTree

import abaqus2exodus as mesh
from prepare_aqueduct_case import sha

SIGNS = np.array([[-1,-1,-1], [1,-1,-1], [1,1,-1], [-1,1,-1],
                  [-1,-1,1], [1,-1,1], [1,1,1], [-1,1,1]], dtype=float)


def hex_volumes(x):
    """Eight Gauss-point Jacobians and integrated volume in original node order."""
    gradients = []
    for xi, eta, zeta in SIGNS / np.sqrt(3):
        sx, sy, sz = SIGNS.T
        gradients.append(np.stack([sx*(1+sy*eta)*(1+sz*zeta),
                                   sy*(1+sx*xi)*(1+sz*zeta),
                                   sz*(1+sx*xi)*(1+sy*eta)], axis=1) / 8)
    jac = np.einsum('eni,qnj->eqij', x, np.array(gradients))
    det = np.linalg.det(jac)
    return det, det.sum(axis=1)


def locate_embedded(points, hosts, identities, tolerance=1e-8):
    """Only declared host elements are searched; no extrapolation or fallback.

    tolerance is dimensionless natural-coordinate roundoff allowance, not the
    Abaqus exterior tolerance. Out-of-host points are reported for later review.
    """
    centers = hosts.mean(axis=1)
    lower, upper = hosts.min(axis=1), hosts.max(axis=1)
    radii = np.linalg.norm(hosts-centers[:, None, :], axis=2).max(axis=1)
    tree = cKDTree(centers)
    found, missing = [], []
    for node, pt in points:
        candidates = tree.query_ball_point(pt, r=float(radii.max()) + 1e-8)
        candidates.sort(key=lambda j: identities[j])
        chosen = None
        for j in candidates:
            if np.any(pt < lower[j]-1e-8) or np.any(pt > upper[j]+1e-8):
                continue
            natural = mesh.point_in_hex(pt, hosts[j], tol=tolerance)
            if natural is None:
                continue
            weights = np.array(mesh._hex_shape(*natural))
            residual = float(np.max(np.abs(weights @ hosts[j]-pt)))
            if residual > 1e-8 or abs(float(weights.sum())-1) > 1e-12:
                continue
            chosen = dict(node=int(node), host_index=j, host=list(identities[j]),
                          natural=list(natural), weights=weights.tolist(),
                          coordinate_residual_m=residual)
            break
        if chosen is None:
            missing.append(dict(node=int(node), xyz_m=pt.tolist()))
        else:
            found.append(chosen)
    return found, missing


def audit(preparation, output):
    root, output = Path(preparation), Path(output)
    if output.exists():
        raise ValueError('audit output already exists')
    d = json.loads((root / 'assembly-contract.json').read_text())
    for filename, digest in d['artifacts'].items():
        if sha(root / filename) != digest:
            raise ValueError(f'artifact hash mismatch: {filename}')
    nodes = {}
    with (root / 'node-map.csv').open() as f:
        for row in csv.DictReader(f):
            key = (row['instance'], int(row['abaqus_node_label']))
            if key in nodes:
                raise ValueError('duplicate node origin')
            nodes[key] = int(row['exodus_node'])
    origins = {}
    with (root / 'element-map.csv').open() as f:
        for row in csv.DictReader(f):
            key = (row['instance'], int(row['abaqus_element_label']))
            if key in origins:
                raise ValueError('duplicate element origin')
            origins[key] = int(row['exodus_element'])
    connectivity, types = {}, {}
    with netCDF4.Dataset(root / 'assembly-geometry.e') as nc:
        coords = np.column_stack([nc[n][:] for n in ('coordx', 'coordy', 'coordz')])
        labels = nc['elem_num_map'][:]
        offset = 0
        for j in range(1, len(nc.dimensions['num_el_blk']) + 1):
            v = nc[f'connect{j}']; conn = v[:]
            for k, row in enumerate(conn):
                eid = int(labels[offset+k])
                connectivity[eid] = np.array(row, dtype=int)-1
                types[eid] = v.elem_type
            offset += len(conn)
    if len(nodes) != len(coords) or len(set(nodes.values())) != len(coords):
        raise ValueError('node origin mapping is not bijective')
    if set(origins.values()) != set(connectivity):
        raise ValueError('element origin mapping is not bijective')
    if not np.isfinite(coords).all():
        raise ValueError('nonfinite geometry')
    hex_ids = sorted(e for e in types if types[e] == 'HEX8')
    det, volumes = hex_volumes(np.array([coords[connectivity[e]] for e in hex_ids]))
    truss_ids = sorted(e for e in types if types[e] == 'TRUSS')
    lengths = [float(np.linalg.norm(coords[connectivity[e]][1]-coords[connectivity[e]][0]))
               for e in truss_ids]
    shell_ids = sorted(e for e in types if types[e] == 'SHELL4')
    # Triangle area checks detect collapsed shell geometry; no shell formula claim.
    areas = []
    for eid in shell_ids:
        a,b,c,e = coords[connectivity[eid]]
        areas.append(float((np.linalg.norm(np.cross(b-a,c-a)) +
                            np.linalg.norm(np.cross(c-a,e-a))) / 2))
    groups = []
    output.mkdir(parents=True)
    for j, group in enumerate(d['interfaces']['embedded'], 1):
        host_origins = [tuple(x) for x in group['host']]
        host_ids = [origins[x] for x in host_origins]
        if any(types[e] != 'HEX8' for e in host_ids):
            raise ValueError('embedded host is not HEX8')
        embedded_ids = [origins[tuple(x)] for x in group['embedded']]
        if any(types[e] != 'TRUSS' for e in embedded_ids):
            raise ValueError('embedded child is not TRUSS')
        node_ids = sorted(set(int(n) for e in embedded_ids for n in connectivity[e]))
        found, missing = locate_embedded([(n+1, coords[n]) for n in node_ids],
                                         np.array([coords[connectivity[e]] for e in host_ids]),
                                         host_origins)
        path = output / f'embedded-{j:02d}.json'
        path.write_text(json.dumps(dict(found=found, unresolved=missing), indent=2)+'\n')
        groups.append(dict(source_line=group['source_line'], host_elements=len(host_ids),
                           embedded_elements=len(embedded_ids), nodes=len(node_ids),
                           located=len(found), unresolved=len(missing), file=path.name,
                           sha256=sha(path),
                           max_coordinate_residual_m=max((x['coordinate_residual_m']
                                                          for x in found), default=None)))
    summary = dict(status='offline-geometry-and-embedding-audit', no_solver_executed=True,
                   contract_sha256=sha(root / 'assembly-contract.json'),
                   node_map_bijective=True, element_map_bijective=True,
                   nodes=len(coords), mesh_elements=len(connectivity),
                   element_types=dict(Counter(types.values())),
                   hex_gauss_points=det.size, min_hex_detJ_m3=float(det.min()),
                   nonpositive_hex_gauss_points=int((det <= 0).sum()),
                   integrated_hex_volume_m3=float(volumes.sum()),
                   minimum_truss_length_m=min(lengths, default=None),
                   zero_length_trusses=sum(v <= 0 for v in lengths),
                   minimum_shell_area_m2=min(areas, default=None),
                   zero_area_shells=sum(v <= 0 for v in areas), embedded=groups,
                   embedding_policy='declared_host_only_no_extrapolation_no_nearest_node',
                   natural_coordinate_tolerance=1e-8, coordinate_tolerance_m=1e-8,
                   does_not_prove='Abaqus exterior-tolerance, embedded/Tie stiffness or FE response equivalence')
    (output / 'geometry-audit.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preparation', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.preparation, a.output), indent=2))
