#!/usr/bin/env python3
"""Build the approved 2026-08-30 HEX8 diagnostic, never reading Abaqus results.

This is deliberately source-hash locked, not a general kinematic-coupling converter.
Run with --inp EXPERT.inp --output NEW_DIRECTORY. Requires numpy/netCDF4.
The output is a complete LIMS multipart snapshot; this script does not submit jobs.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from abaqus2exodus import build_global_mesh, parse_inp, write_exodus
from abaqus_cdp import choose_material, parse_materials, validate_material, write_bundle

SOURCE_SHA = "5fbd5e1697874a1b7307173b1d60588bc3d7de66df1f65f5fca2dfb0e7bfbfd7"
CASE = "tpl-bj-uniaxial-tension-cdp-hex8-v1"
INPUT = CASE + ".i"
MESH = "uniaxial_tension_mesh.e"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def material_si(source: Path) -> tuple[dict, dict]:
    material = choose_material(parse_materials(source), "concrete")
    original = validate_material(material)
    converted = copy.deepcopy(material)
    for key in ("elastic", "compression_hardening", "tension_stiffening"):
        for _, row in converted.blocks[key].rows:
            row[0] *= 1e6
    return original, validate_material(converted)


def make_mesh(source: Path, destination: Path, *, expected_type: str = "C3D8R",
              expected_displacement: float = 0.025, title: str | None = None) -> dict:
    model = parse_inp(source)
    check(len(model.parts) == len(model.instances) == len(model.steps) == 1, "Expected one part/instance/step")
    part = model.parts["concrete_cube"]
    check(len(part.nodes) == 1331 and len(part.elems) == 1000, "Unexpected topology")
    check(set(part.elem_types.values()) == {expected_type}, f"Expected {expected_type}")
    check(model.steps[0]["static"] == [0.01, 1.0, 1e-15, 1.0], "Unexpected Static parameters")
    boundaries = model.steps[0]["boundaries"]
    check(all(b["op"] == "NEW" for b in boundaries), "OP=NEW must be preserved")
    translations = {(b["set"], b["dof1"], b["dof2"], b["value"])
                    for b in boundaries if b["dof1"] <= 3}
    check(translations == {("_PickedSet14", 3, 3, expected_displacement), ("_PickedSet16", 3, 3, 0.0)},
          "Unexpected translations: do not silently clamp transverse DOFs")
    gm, blocks, types, meta, node_sets, _ = build_global_mesh(model, 1e-9)
    check(len(gm.coords) == 1331, "No disconnected RP node may enter the solve")
    check(all(gm.node_map[("concrete_cube-1", n)] == n for n in part.nodes), "Node IDs changed")
    check(all(origin == ("concrete_cube-1", eid) for eid, origin in gm.elem_origin.items()), "Element IDs changed")
    block = next(iter(blocks))
    check(all(conn == part.elems[eid] for eid, conn in blocks[block]), "Connectivity changed")
    gm.coords = [tuple(c * 1e-3 for c in xyz) for xyz in gm.coords]
    xyz = np.asarray(gm.coords)
    top = node_sets["SURF__PickedSurf8"]
    bottom = node_sets["_PickedSet16__concrete_cube_"]
    check(len(top) == len(bottom) == 121, "Expected 121 nodes per face")
    check(np.allclose(xyz[np.asarray(top)-1, 2], .2), "Surface S3 is not the top")
    check(np.allclose(xyz[np.asarray(bottom)-1, 2], .05), "Bottom node set mismatch")
    # Verify the mapped HEX8 Jacobian at all eight Gauss points, not just its bounds.
    signs = np.asarray([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                        [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
    cells = xyz[np.asarray([c for _, c in blocks[block]]) - 1]
    determinants = []
    for point in signs / np.sqrt(3):
        factors = 1 + signs * point
        grad = np.column_stack([signs[:, k] * factors[:, (k+1)%3] * factors[:, (k+2)%3] / 8
                                for k in range(3)])
        determinants.append(np.linalg.det(np.einsum("eni,nj->eij", cells, grad)))
    determinants = np.asarray(determinants)
    check(np.all(determinants > 0), "Nonpositive Jacobian")
    check(np.isclose(determinants.sum(), .15**3, rtol=1e-10), "Wrong total volume")
    write_exodus(destination, gm, blocks, types, meta, {"top": top, "bottom": bottom}, {},
                  title or "uniaxial_tension: mm->m; original IDs; diagnostic HEX8, not C3D8R")
    return {"nodes": 1331, "elements": 1000, "block": block, "element_id_mapping": "identity",
            "node_id_mapping": "identity", "minimum_gauss_detJ_m3": float(determinants.min()),
            "volume_m3": float(determinants.sum()), "bounds_m": [xyz.min(0).tolist(), xyz.max(0).tolist()],
            "top_nodes": top, "bottom_nodes": bottom}


def moose_input(si: dict, fragment: str) -> str:
    block = "concrete_cube__concrete"
    properties = {"DamageC": "DamageC", "DamageT": "DamageT", "kappa_c": "cdp_kappa_c",
                  "kappa_t": "cdp_kappa_t", "local_iterations": "cdp_local_iterations",
                  "accepted_substeps": "cdp_accepted_substeps",
                  "jacobian_fallbacks": "cdp_jacobian_fallbacks",
                  "integration_microseconds": "cdp_integration_microseconds"}
    out = [f"""# Template ID : {CASE}
# Status      : prototype diagnostic; not engineering acceptance
# Source INP  : uniaxial_tension.inp SHA256={SOURCE_SHA}
# Units       : user-approved mm,N,s,MPa -> m,N,s,Pa; lengths x1e-3, stresses x1e6.
# Loading     : original *Static 0.01,1,1e-15,1; U3=+0.025mm -> +2.5e-5*t m.
# Mapping     : OP=NEW clears initial clamps. Bottom U3=0 only; bottom U1/U2 free.
# Mapping     : top kinematic coupling with zero RP rotations gives uniform top translations.
# Mapping     : top U1=U2=0 fixes two free rigid-body translation gauges, no bottom XY clamp.
# Mapping     : C3D8R IDs/connectivity retained as standard HEX8, not reduced/hourglass equivalent.
# Mapping     : four material tables are from INP only, no reference result CSV is a solver input.
# Mapping     : default recovery wt=0,wc=1; viscosity=0.0005s; no automatic stabilization.
# Mapping     : IterationAdaptiveDT differs from Abaqus algorithm; original increment bounds retained.
# Mapping     : elemental fields are quadrature-projected constants, not raw C3D8R IntPt=1 values.
# Mapping     : reference CSV has no RP RF3/U3, so the RP history cannot yet be accepted by comparison.

[Mesh/file]
  type = FileMeshGenerator
  file = {MESH}
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[Physics/SolidMechanics/QuasiStatic/concrete]
  add_variables = true
  incremental = true
  block = {block}
  strain = SMALL
  generate_output = 'stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz
                     strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz
                     max_principal_stress mid_principal_stress min_principal_stress vonmises_stress'
  save_in = 'resid_x resid_y resid_z'
[]
[AuxVariables]
  [resid_x]
  []
  [resid_y]
  []
  [resid_z]
  []
"""]
    for name in properties:
        out.append(f"  [{name}]\n    order = CONSTANT\n    family = MONOMIAL\n  []\n")
    out.append("[]\n[AuxKernels]\n")
    for name, prop in properties.items():
        out.append(f"  [{name}]\n    type = MaterialRealAux\n    variable = {name}\n    property = {prop}\n    block = {block}\n    execute_on = 'initial timestep_end'\n  []\n")
    out.append("[]\n[Functions/top_displacement]\n  type = ParsedFunction\n  expression = '2.5e-5*t'\n[]\n[BCs]\n")
    for name, var, face in [("bottom_z", "disp_z", "bottom"), ("top_x_gauge", "disp_x", "top"),
                             ("top_y_gauge", "disp_y", "top")]:
        out.append(f"  [{name}]\n    type = DirichletBC\n    variable = {var}\n    boundary = {face}\n    value = 0\n  []\n")
    out.append("  [top_z]\n    type = FunctionDirichletBC\n    variable = disp_z\n    boundary = top\n    function = top_displacement\n  []\n[]\n")
    # Fragment is a valid standalone audit artifact; merge its Materials entries explicitly.
    inner = fragment.split("[Materials]\n", 1)[1].rsplit("[]", 1)[0]
    out.append(f"""[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    block = {block}
    youngs_modulus = {si['elastic']['youngs_modulus_pa']:.17g}
    poissons_ratio = {si['elastic']['poissons_ratio']:.17g}
  []
  [stress]
    type = ComputeMultipleInelasticStress
    block = {block}
    inelastic_models = cdp_stress_update
    perform_finite_strain_rotations = false
  []
{inner}[]
[Postprocessors]
  [RP1_Force]
    type = NodalSum
    variable = resid_z
    boundary = top
  []
  [Bottom_Force]
    type = NodalSum
    variable = resid_z
    boundary = bottom
  []
  [Top_Force_X]
    type = NodalSum
    variable = resid_x
    boundary = top
  []
  [Top_Force_Y]
    type = NodalSum
    variable = resid_y
    boundary = top
  []
  [RP1_Displacement]
    type = AverageNodalVariableValue
    variable = disp_z
    boundary = top
  []
""")
    for name, var in [("damagec", "DamageC"), ("damaget", "DamageT"), ("mises", "vonmises_stress"),
                      ("stress_zz", "stress_zz"), ("local_iterations", "local_iterations"),
                      ("accepted_substeps", "accepted_substeps"), ("jacobian_fallbacks", "jacobian_fallbacks")]:
        out.append(f"  [max_{name}]\n    type = ElementExtremeValue\n    variable = {var}\n    value_type = max\n  []\n")
    out.append("""[]
[Preconditioning/smp]
  type = SMP
  full = true
[]
[Executioner]
  type = Transient
  start_time = 0
  end_time = 1
  solve_type = NEWTON
  line_search = bt
  automatic_scaling = true
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 50
  num_steps = 100000
  dtmin = 1e-15
  dtmax = 1
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 8
    iteration_window = 3
    growth_factor = 1.15
    cutback_factor = 0.5
  []
[]
[Times/field_output_times]
  type = TimeIntervalTimes
  start_time = 0
  end_time = 1
  time_interval = 0.01
[]
[Outputs]
  [field_exodus]
    type = Exodus
    execute_on = 'initial timestep_end'
    sync_times_object = field_output_times
    sync_only = true
    file_base = uniaxial_tension
  []
  [history_csv]
    type = CSV
    execute_on = 'initial timestep_end'
    sync_times_object = field_output_times
    sync_only = true
    file_base = uniaxial_tension
  []
[]
""")
    return "".join(out)


def build(source: Path, output: Path) -> dict:
    check(sha(source) == SOURCE_SHA, "Source changed: re-review units, coupling, boundaries and tables first")
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "input"
    inputs.mkdir()
    raw, si = material_si(source)
    manifest = write_bundle(source, inputs, si)
    mesh = make_mesh(source, inputs / MESH)
    fragment_path = inputs / "abaqus_cdp_material.i"
    fragment = fragment_path.read_text().replace("fb0_fc0 =", "biaxial_to_uniaxial_compression_ratio =").replace("kc =", "tensile_meridian_ratio =")
    fragment = fragment.replace("    type = AbaqusCDPStressUpdate", "    type = AbaqusCDPStressUpdate\n    block = concrete_cube__concrete\n    maximum_substeps = 256\n    maximum_strain_increment = 2.5e-5\n    enable_performance_diagnostics = true")
    fragment_path.write_text(fragment, encoding="utf-8")
    (inputs / INPUT).write_text(moose_input(si, fragment), encoding="utf-8")
    # Copy only the original model. The Abaqus output CSV is intentionally never opened here.
    (inputs / source.name).write_bytes(source.read_bytes())
    manifest.update({"status": "prototype / user-approved HEX8 diagnostic",
        "source_units": {"length": "mm", "stress": "MPa", "force": "N", "time": "s"},
        "unit_conversion": {"length_displacement_factor": .001, "stress_modulus_factor": 1e6,
                            "strain_damage_time_factor": 1, "approval": "user confirmation 2026-08-30"},
        "source_elastic_raw": raw["elastic"], "mesh_mapping": mesh,
        "builder": {"file": Path(__file__).name, "sha256": sha(Path(__file__))},
        "assumptions": ["用户批准源单位 mm,N,s,MPa；长度与位移乘1e-3，应力与E乘1e6；应变/损伤/时间不变",
            "保留1000单元/1331节点原ID和连通；装配平移也换算到m",
            "OP=NEW：底面仅U3=0；顶面刚性运动学耦合、RP转角为0",
            "顶面U1=U2=0仅固定两个自由整体平移零基准；底面横向自由",
            "默认刚度恢复wt=0,wc=1；保持黏性0.0005s；未加入额外稳定化",
            "材料四表只来自原始INP；piecewise-linear interpolation；constant endpoint extrapolation",
            "保持原Static初始/总时长/最小/最大增量0.01/1/1e-15/1；自适应算法与Abaqus不同",
            "每0.01s同步输出；所有场为真实求解，无结果CSV插值驱动"],
        "unsupported": ["Abaqus C3D8R减缩积分与沙漏控制未等价复现",
            "HEX8多积分点投影单元常数场不能视为与C3D8R IntPt=1逐积分点等价",
            "专家场CSV未包含RF3/U3历史；不可宣称力位移已验证",
            "本次是诊断，不是自定义本构与Abaqus完全等价或工程验收"],
        "artifacts": {p.name: sha(p) for p in sorted(inputs.iterdir()) if p.name != "cdp-mapping-manifest.json"}})
    (inputs / "cdp-mapping-manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    files = {p.name: sha(p) for p in sorted(inputs.iterdir())}
    submission = {"project_id": "gmp-ise", "case_name": CASE,
        "input_file": INPUT, "input_sha256": files[INPUT],
        "mesh_files": [{"name": MESH, "sha256": files[MESH]}],
        "extra_files": [{"name": n, "sha256": h} for n,h in files.items() if n not in (INPUT,MESH)],
        "command": f"DamSafetyApp-opt -i {INPUT}"}
    (output / "submission.json").write_text(json.dumps(submission, indent=2)+"\n")
    return {"output": str(output), "source_sha256": SOURCE_SHA, "input_sha256": files[INPUT],
            "files": len(files), "nodes": mesh["nodes"], "elements": mesh["elements"],
            "table_counts": manifest["table_point_counts"], "elastic_si": si["elastic"],
            "recovery": si["recovery"], "end_time": 1, "end_displacement_m": 2.5e-5}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--inp", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    print(json.dumps(build(args.inp.resolve(), args.output.resolve()), ensure_ascii=False, indent=2))
