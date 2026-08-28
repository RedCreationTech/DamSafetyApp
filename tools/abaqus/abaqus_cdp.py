#!/usr/bin/env python3
"""Extract and validate Abaqus Concrete Damaged Plasticity material data.

The converter deliberately treats the Abaqus input deck as the source of material
data. ODB/CSV result exports are validation evidence and are never read here.

P0 supports one elastic row, the five scalar CDP parameters, and strain-based
two-column compression/tension stress and damage tables without temperature,
field-variable, or strain-rate dependencies.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


CDP_KEYWORDS = {
    "elastic": "elastic",
    "concrete damaged plasticity": "plasticity",
    "concrete compression hardening": "compression_hardening",
    "concrete compression damage": "compression_damage",
    "concrete tension stiffening": "tension_stiffening",
    "concrete tension damage": "tension_damage",
}

TABLE_SPECS = {
    "compression_hardening": ("stress_pa", "inelastic_strain"),
    "compression_damage": ("damage_c", "inelastic_strain"),
    "tension_stiffening": ("stress_pa", "cracking_strain"),
    "tension_damage": ("damage_t", "cracking_strain"),
}


class CDPInputError(ValueError):
    """A traceable, user-facing CDP input validation error."""


@dataclass(frozen=True)
class SourceLine:
    source: Path
    number: int
    text: str


@dataclass
class KeywordBlock:
    keyword: str
    options: dict[str, str | bool]
    source: Path
    line: int
    rows: list[tuple[int, list[float]]] = field(default_factory=list)


@dataclass
class CDPMaterial:
    name: str
    blocks: dict[str, KeywordBlock] = field(default_factory=dict)

    def require(self, key: str) -> KeywordBlock:
        try:
            return self.blocks[key]
        except KeyError as exc:
            raise CDPInputError(
                f"材料 {self.name!r} 缺少 *{_display_keyword(key)}") from exc


def _display_keyword(key: str) -> str:
    for keyword, canonical in CDP_KEYWORDS.items():
        if canonical == key:
            return keyword.title()
    return key


def _parse_keyword(line: str) -> tuple[str, dict[str, str | bool]]:
    tokens = [token.strip() for token in line[1:].split(",")]
    keyword = tokens[0].lower()
    options: dict[str, str | bool] = {}
    for token in tokens[1:]:
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            options[key.strip().lower()] = value.strip().strip("\"'")
        else:
            options[token.lower()] = True
    return keyword, options


def _expanded_lines(path: Path, stack: tuple[Path, ...] = ()) -> list[SourceLine]:
    source = path.expanduser().resolve()
    if source in stack:
        chain = " -> ".join(str(item) for item in (*stack, source))
        raise CDPInputError(f"*Include 循环引用: {chain}")
    if not source.is_file():
        raise CDPInputError(f"Abaqus 输入文件不存在: {source}")

    result: list[SourceLine] = []
    for number, raw in enumerate(
        source.read_text(encoding="utf-8", errors="ignore").splitlines(), 1
    ):
        stripped = raw.strip()
        if stripped.lower().startswith("*include"):
            _, options = _parse_keyword(stripped)
            child = options.get("input") or options.get("file")
            if not isinstance(child, str) or not child:
                raise CDPInputError(f"{source}:{number}: *Include 缺少 input/file")
            result.extend(_expanded_lines(source.parent / child, (*stack, source)))
        else:
            result.append(SourceLine(source, number, raw))
    return result


def _numeric_row(item: SourceLine) -> list[float]:
    values = []
    for token in item.text.rstrip(",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError as exc:
            raise CDPInputError(
                f"{item.source}:{item.number}: 非数值材料数据 {token!r}") from exc
        if not math.isfinite(value):
            raise CDPInputError(
                f"{item.source}:{item.number}: 材料数据必须是有限数值")
        values.append(value)
    return values


def parse_materials(inp: Path) -> dict[str, CDPMaterial]:
    """Parse CDP-related keyword blocks grouped by Abaqus material name."""
    materials: dict[str, CDPMaterial] = {}
    current_material: CDPMaterial | None = None
    current_block: KeywordBlock | None = None

    for item in _expanded_lines(inp):
        stripped = item.text.strip()
        if not stripped or stripped.startswith("**"):
            continue
        if stripped.startswith("*"):
            keyword, options = _parse_keyword(stripped)
            current_block = None
            if keyword == "material":
                name = options.get("name")
                if not isinstance(name, str) or not name:
                    raise CDPInputError(f"{item.source}:{item.number}: *Material 缺少 name")
                current_material = materials.setdefault(name, CDPMaterial(name))
                continue
            canonical = CDP_KEYWORDS.get(keyword)
            if canonical and current_material is not None:
                if canonical in current_material.blocks:
                    previous = current_material.blocks[canonical]
                    raise CDPInputError(
                        f"{item.source}:{item.number}: 材料 {current_material.name!r} "
                        f"重复定义 *{keyword}；首次位于 {previous.source}:{previous.line}")
                current_block = KeywordBlock(
                    keyword=keyword,
                    options=options,
                    source=item.source,
                    line=item.number,
                )
                current_material.blocks[canonical] = current_block
            continue

        if current_block is not None:
            values = _numeric_row(item)
            if values:
                current_block.rows.append((item.number, values))

    return {
        name: material
        for name, material in materials.items()
        if "plasticity" in material.blocks
    }


def choose_material(materials: dict[str, CDPMaterial], name: str | None) -> CDPMaterial:
    if name:
        try:
            return materials[name]
        except KeyError as exc:
            choices = ", ".join(sorted(materials)) or "无"
            raise CDPInputError(f"未找到 CDP 材料 {name!r}；可选材料: {choices}") from exc
    if len(materials) != 1:
        choices = ", ".join(sorted(materials)) or "无"
        raise CDPInputError(f"必须用 --material 选择 CDP 材料；可选材料: {choices}")
    return next(iter(materials.values()))


def _require_shape(block: KeywordBlock, rows: int | None, columns: int) -> None:
    if rows is not None and len(block.rows) != rows:
        raise CDPInputError(
            f"{block.source}:{block.line}: *{block.keyword} 需要 {rows} 行，"
            f"实际 {len(block.rows)} 行")
    if not block.rows:
        raise CDPInputError(f"{block.source}:{block.line}: *{block.keyword} 没有数据")
    for line, values in block.rows:
        if len(values) != columns:
            raise CDPInputError(
                f"{block.source}:{line}: P0 的 *{block.keyword} 仅支持 {columns} 列，"
                f"实际 {len(values)} 列；温度/场变量/率依赖尚不支持")


def _table(block: KeywordBlock) -> list[tuple[float, float]]:
    _require_shape(block, None, 2)
    return [(values[0], values[1]) for _, values in block.rows]


def _validate_abscissa(name: str, table: list[tuple[float, float]]) -> None:
    if table[0][1] != 0:
        raise CDPInputError(f"{name} 首行横坐标必须为 0，实际 {table[0][1]}")
    previous = None
    for index, (_, x_value) in enumerate(table, 1):
        if x_value < 0:
            raise CDPInputError(f"{name} 第 {index} 行横坐标为负: {x_value}")
        if previous is not None and x_value <= previous:
            raise CDPInputError(
                f"{name} 横坐标必须严格递增；第 {index} 行 {x_value} <= {previous}")
        previous = x_value


def _interpolate(table: list[tuple[float, float]], x_value: float) -> float:
    if x_value <= table[0][1]:
        return table[0][0]
    for (y0, x0), (y1, x1) in zip(table, table[1:]):
        if x_value <= x1:
            ratio = (x_value - x0) / (x1 - x0)
            return y0 + ratio * (y1 - y0)
    return table[-1][0]


def _plastic_strains(
    stress_table: list[tuple[float, float]],
    damage_table: list[tuple[float, float]],
    youngs_modulus: float,
    label: str,
) -> list[dict[str, float]]:
    derived = []
    previous = -math.inf
    for index, (stress, strain_measure) in enumerate(stress_table, 1):
        damage = _interpolate(damage_table, strain_measure)
        if not 0 <= damage < 1:
            raise CDPInputError(
                f"{label} 第 {index} 行插值得到非法损伤 {damage}，要求 0 <= d < 1")
        plastic_strain = strain_measure - damage / (1 - damage) * stress / youngs_modulus
        tolerance = 1e-14 * max(1.0, abs(strain_measure))
        if plastic_strain < -tolerance:
            raise CDPInputError(
                f"{label} 第 {index} 行得到负等效塑性应变 {plastic_strain:.16g}")
        plastic_strain = max(0.0, plastic_strain)
        if plastic_strain + tolerance < previous:
            raise CDPInputError(
                f"{label} 第 {index} 行等效塑性应变下降: "
                f"{plastic_strain:.16g} < {previous:.16g}")
        previous = plastic_strain
        derived.append(
            {
                "strain_measure": strain_measure,
                "stress_pa": stress,
                "damage": damage,
                "equivalent_plastic_strain": plastic_strain,
            }
        )
    return derived


def validate_material(material: CDPMaterial) -> dict[str, object]:
    elastic = material.require("elastic")
    plasticity = material.require("plasticity")
    _require_shape(elastic, 1, 2)
    _require_shape(plasticity, 1, 5)

    for key in TABLE_SPECS:
        block = material.require(key)
        dependencies = block.options.get("dependencies", "0")
        if str(dependencies) not in ("0", "0.0"):
            raise CDPInputError(
                f"{block.source}:{block.line}: P0 不支持 *{block.keyword} DEPENDENCIES"
            )
    for key in ("tension_stiffening", "tension_damage"):
        block = material.blocks[key]
        table_type = str(block.options.get("type", "strain")).lower()
        if table_type != "strain":
            raise CDPInputError(
                f"{block.source}:{block.line}: P0 仅支持 *{block.keyword}, TYPE=STRAIN"
            )

    youngs_modulus, poissons_ratio = elastic.rows[0][1]
    if youngs_modulus <= 0:
        raise CDPInputError("Young's modulus 必须大于 0")
    if not -1 < poissons_ratio < 0.5:
        raise CDPInputError("Poisson's ratio 必须满足 -1 < nu < 0.5")

    dilation_angle, eccentricity, fb0_fc0, kc, viscosity = plasticity.rows[0][1]
    if not 0 <= dilation_angle < 90:
        raise CDPInputError("dilation angle 必须位于 [0, 90) 度")
    if eccentricity <= 0 or fb0_fc0 <= 0 or not 0.5 < kc <= 1 or viscosity < 0:
        raise CDPInputError("eccentricity、fb0/fc0、Kc 或 viscosity 超出 P0 合法范围")

    tables = {key: _table(material.require(key)) for key in TABLE_SPECS}
    for key, table in tables.items():
        _validate_abscissa(key, table)
        if key.endswith("damage"):
            if table[0][0] != 0:
                raise CDPInputError(
                    f"{key} 首行损伤必须为 0，实际 {table[0][0]}")
            previous_damage = -math.inf
            for index, (damage, _) in enumerate(table, 1):
                if not 0 <= damage < 1:
                    raise CDPInputError(
                        f"{key} 第 {index} 行损伤 {damage} 不满足 0 <= d < 1")
                if damage < previous_damage:
                    raise CDPInputError(
                        f"{key} 第 {index} 行损伤下降: {damage} < {previous_damage}")
                previous_damage = damage
        else:
            for index, (stress, _) in enumerate(table, 1):
                if stress <= 0:
                    raise CDPInputError(f"{key} 第 {index} 行应力必须为正绝对值")

    compression_plastic = _plastic_strains(
        tables["compression_hardening"],
        tables["compression_damage"],
        youngs_modulus,
        "compression",
    )
    tension_plastic = _plastic_strains(
        tables["tension_stiffening"],
        tables["tension_damage"],
        youngs_modulus,
        "tension",
    )

    compression_options = material.blocks["compression_damage"].options
    tension_options = material.blocks["tension_damage"].options
    return {
        "material": material.name,
        "elastic": {"youngs_modulus_pa": youngs_modulus, "poissons_ratio": poissons_ratio},
        "plasticity": {
            "dilation_angle_degrees": dilation_angle,
            "eccentricity": eccentricity,
            "fb0_fc0": fb0_fc0,
            "kc": kc,
            "viscosity_seconds": viscosity,
        },
        "recovery": {
            "tension_recovery": _option_float(compression_options, "tension recovery", 0.0),
            "compression_recovery": _option_float(
                tension_options, "compression recovery", 1.0
            ),
        },
        "tables": tables,
        "derived": {
            "compression": compression_plastic,
            "tension": tension_plastic,
        },
        "provenance": {
            key: {
                "keyword": block.keyword,
                "source": str(block.source),
                "keyword_line": block.line,
                "data_lines": [line for line, _ in block.rows],
                "options": block.options,
            }
            for key, block in material.blocks.items()
        },
    }


def _option_float(options: dict[str, str | bool], key: str, default: float) -> float:
    value = options.get(key, default)
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise CDPInputError(f"关键字选项 {key!r} 不是数值: {value!r}") from exc
    if not 0 <= result <= 1:
        raise CDPInputError(f"关键字选项 {key!r} 必须位于 [0, 1]")
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_table(path: Path, headers: tuple[str, str], rows: Iterable[tuple[float, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(headers)
        for row in rows:
            writer.writerow((format(row[0], ".17g"), format(row[1], ".17g")))


def _moose_fragment(validated: dict[str, object], files: dict[str, str]) -> str:
    elastic = validated["elastic"]
    plasticity = validated["plasticity"]
    recovery = validated["recovery"]
    return f"""# Generated prototype input fragment. Requires DamSafetyApp AbaqusCDPStressUpdate.
[Materials]
  [cdp_stress_update]
    type = AbaqusCDPStressUpdate
    youngs_modulus = {elastic['youngs_modulus_pa']:.17g}
    poissons_ratio = {elastic['poissons_ratio']:.17g}
    dilation_angle = {plasticity['dilation_angle_degrees']:.17g}
    eccentricity = {plasticity['eccentricity']:.17g}
    fb0_fc0 = {plasticity['fb0_fc0']:.17g}
    kc = {plasticity['kc']:.17g}
    viscosity = {plasticity['viscosity_seconds']:.17g}
    tension_recovery = {recovery['tension_recovery']:.17g}
    compression_recovery = {recovery['compression_recovery']:.17g}
    compression_hardening_file = {files['compression_hardening']}
    compression_damage_file = {files['compression_damage']}
    tension_stiffening_file = {files['tension_stiffening']}
    tension_damage_file = {files['tension_damage']}
  []
[]
"""


def write_bundle(inp: Path, out_dir: Path, validated: dict[str, object]) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, str] = {}
    for key, headers in TABLE_SPECS.items():
        filename = f"{key}.csv"
        _write_table(out_dir / filename, headers, validated["tables"][key])
        files[key] = filename

    fragment_name = "abaqus_cdp_material.i"
    (out_dir / fragment_name).write_text(
        _moose_fragment(validated, files), encoding="utf-8"
    )

    source = inp.expanduser().resolve()
    source_files = []
    seen_sources = set()
    for item in _expanded_lines(source):
        if item.source not in seen_sources:
            seen_sources.add(item.source)
            source_files.append(
                {"path": str(item.source), "sha256": _sha256(item.source)}
            )
    artifact_hashes = {
        name: _sha256(out_dir / name) for name in (*files.values(), fragment_name)
    }
    manifest = {
        "schema_version": "cdp-mapping-manifest/0.1",
        "status": "prototype",
        "source": {"path": str(source), "sha256": _sha256(source)},
        "source_files": source_files,
        "unit_system": {"length": "m", "force": "N", "time": "s", "stress": "Pa"},
        "material": validated["material"],
        "elastic": validated["elastic"],
        "plasticity": validated["plasticity"],
        "recovery": validated["recovery"],
        "table_point_counts": {
            key: len(validated["tables"][key]) for key in TABLE_SPECS
        },
        "derived_plastic_strain": validated["derived"],
        "provenance": validated["provenance"],
        "artifacts": artifact_hashes,
        "assumptions": [
            "SI consistent units supplied by project decision",
            "piecewise-linear interpolation",
            "constant endpoint extrapolation",
        ],
        "unsupported": [
            "temperature/field-variable/rate-dependent material tables",
            "TYPE=DISPLACEMENT tension stiffening",
            "Abaqus automatic stabilization equivalence",
            "C3D8R hourglass equivalence",
        ],
    }
    manifest_path = out_dir / "cdp-mapping-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inp", required=True, type=Path)
    parser.add_argument("--material")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--check", action="store_true", help="仅解析和验证，不写文件")
    args = parser.parse_args()

    try:
        material = choose_material(parse_materials(args.inp), args.material)
        validated = validate_material(material)
        if not args.check:
            if args.out_dir is None:
                parser.error("非 --check 模式必须提供 --out-dir")
            write_bundle(args.inp, args.out_dir, validated)
        summary = {
            "status": "valid",
            "material": material.name,
            "table_point_counts": {
                key: len(validated["tables"][key]) for key in TABLE_SPECS
            },
        }
        print(json.dumps(summary, ensure_ascii=False))
        return 0
    except CDPInputError as exc:
        parser.exit(2, f"CDP input error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
