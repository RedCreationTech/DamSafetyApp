#!/usr/bin/env python3
"""Hash-locked 2026-08-31 C3D8 compression -> custom-CDP HEX8 input. No result CSV input."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import build_uniaxial_tension_diagnostic as base
from abaqus_cdp import write_bundle
from c3d8_formulation import apply_c3d8_bbar

SOURCE_SHA = "032449357ed66108b0f11967d62cb729ebf1c970d00efbb4ba41ef5bf00a1e26"
CASE = "tpl-bj-uniaxial-compression-cdp-c3d8-hex8-v1"
INPUT = CASE + ".i"
MESH = "uniaxial_compression_mesh.e"
SOLVER_SHA = "f99f18d1e24e5168ef6b88857d06b8c01f0c0fee"


def moose_input(si: dict, fragment: str) -> str:
    text = base.moose_input(si, fragment)
    for old, new in [
        (base.CASE, CASE), (base.SOURCE_SHA, SOURCE_SHA),
        ("uniaxial_tension", "uniaxial_compression"),
        ("U3=+0.025mm -> +2.5e-5*t m", "U3=-2.5mm -> -0.0025*t m"),
        ("expression = '2.5e-5*t'", "expression = '-0.0025*t'"),
        ("C3D8R IDs/connectivity retained as standard HEX8, not reduced/hourglass equivalent.",
         "C3D8 IDs/connectivity retained as HEX8 with element-averaged volumetric strain (B-bar)."),
        ("not raw C3D8R IntPt=1 values.", "compare with C3D8 volume-averaged 8-IP fields, not raw IP maxima."),
    ]:
        base.check(old in text, f"Shared template changed; review replacement: {old}")
        text = text.replace(old, new)
    text = apply_c3d8_bbar(text)
    marker = "[Postprocessors]\n"
    base.check(text.count(marker) == 1, "Postprocessor anchor changed")
    return text.replace(marker, marker + "  [min_stress_zz]\n    type = ElementExtremeValue\n"
                        "    variable = stress_zz\n    value_type = min\n  []\n")


def build(source: Path, output: Path) -> dict:
    base.check(base.sha(source) == SOURCE_SHA, "Source changed: re-review units, coupling, boundaries and tables first")
    raw, si = base.material_si(source)
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "input"
    inputs.mkdir()
    manifest = write_bundle(source, inputs, si)
    mesh = base.make_mesh(source, inputs / MESH, expected_type="C3D8", expected_displacement=-2.5,
                          title="uniaxial_compression: mm->m; original C3D8 IDs; standard HEX8")
    fragment_path = inputs / "abaqus_cdp_material.i"
    fragment = fragment_path.read_text().replace("fb0_fc0 =", "biaxial_to_uniaxial_compression_ratio =").replace("kc =", "tensile_meridian_ratio =")
    fragment = fragment.replace("    type = AbaqusCDPStressUpdate", "    type = AbaqusCDPStressUpdate\n    block = concrete_cube__concrete\n    maximum_substeps = 256\n    maximum_strain_increment = 2.5e-5\n    enable_performance_diagnostics = true")
    fragment_path.write_text(fragment, encoding="utf-8")
    (inputs / INPUT).write_text(moose_input(si, fragment), encoding="utf-8")
    (inputs / source.name).write_bytes(source.read_bytes())
    manifest.update({
        "status": "prototype / expert C3D8 compression diagnostic",
        "source_units": {"length": "mm", "stress": "MPa", "force": "N", "time": "s"},
        "unit_conversion": {"length_displacement_factor": .001, "stress_modulus_factor": 1e6,
                            "strain_damage_time_factor": 1,
                            "approval": "2026-08-30 confirmed mm/MPa; 2026-08-31 same geometry/material compression counterpart"},
        "source_elastic_raw": raw["elastic"], "mesh_mapping": mesh,
        "loading": {"source_U3_mm": -2.5, "target_U3_m": -.0025, "end_time_s": 1},
        "builder": {"file": Path(__file__).name, "sha256": base.sha(Path(__file__)),
                    "shared_file": Path(base.__file__).name, "shared_sha256": base.sha(Path(base.__file__)),
                    "source_state": "local reviewed converter changes; solver binary remains locked"},
        "assumptions": [
            "单位沿用已确认mm,N,s,MPa；长度/位移乘1e-3，应力/E乘1e6；应变/损伤/时间不变",
            "原C3D8网格1000单元1331节点，原ID/连通/装配平移保留；HEX8启用physics层B-bar体积应变修正",
            "OP=NEW清除初始夹持：底面仅U3=0，横向自由；顶面RP零转角对应刚性平移",
            "顶面U1=U2=0固定两个自由整体平移零基准，不添加底部横向约束",
            "U3=-2.5mm在1s内线性加载；本次不是拉伸+0.025mm，也不是旧压缩-5mm",
            "材料四表仅来自本INP；默认wt=0,wc=1；黏性0.0005s；没有额外稳定化",
            "原Static初始/总时长/最小/最大增量0.01/1/1e-15/1保留，自适应算法非Abaqus同算法",
            "每0.01s输出真实求解场；压缩负应力最小值另记min_stress_zz；无结果CSV驱动",
        ],
        "element_formulation": {"source": "C3D8", "moose": "HEX8", "volumetric_locking_correction": True, "scope": "incremental small strain", "validation": "see damASR E15/E17 evidence; candidate, not final acceptance"},
        "unsupported": ["B-bar兼容修正不代表全部本构及单元内部实现完全等价",
            "MOOSE元素常量投影不是原始IP值；只能对齐8IP体积平均后对比",
            "参考场CSV没有RF3/U3时，不宣称RP力位移对标通过", "诊断性验证，最终接受由专家评审"],
        "artifacts": {p.name: base.sha(p) for p in sorted(inputs.iterdir()) if p.name != "cdp-mapping-manifest.json"},
    })
    (inputs / "cdp-mapping-manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n")
    files = {p.name: base.sha(p) for p in sorted(inputs.iterdir())}
    submission = {"project_id": "gmp-ise", "case_name": CASE, "input_file": INPUT,
        "input_sha256": files[INPUT], "mesh_files": [{"name": MESH, "sha256": files[MESH]}],
        "extra_files": [{"name": n, "sha256": h} for n, h in files.items() if n not in (INPUT, MESH)],
        "command": f"DamSafetyApp-opt -i {INPUT}"}
    (output / "submission.json").write_text(json.dumps(submission, indent=2)+"\n")
    return {"source_sha256": SOURCE_SHA, "input_sha256": files[INPUT], "file_count": len(files),
            "mesh": mesh, "table_counts": manifest["table_point_counts"], "elastic_si": si["elastic"],
            "end_time_s": 1, "end_displacement_m": -.0025}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.inp.resolve(), args.output.resolve())
    (args.output / "build-summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "mesh"}, ensure_ascii=False, indent=2))
