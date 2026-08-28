import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[2]
TOOL_PATH = PROJECT_DIR / "tools" / "abaqus" / "abaqus_cdp.py"
SPEC = importlib.util.spec_from_file_location("abaqus_cdp", TOOL_PATH)
CDP = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CDP
SPEC.loader.exec_module(CDP)


VALID_INPUT = """\
*Material, name=CONCRETE
*Elastic
3.04e10, 0.2
*Concrete Damaged Plasticity
36.31, 0.1, 1.16, 0.667, 0.0005
*Concrete Compression Hardening
2.0e7, 0.0
2.2e7, 0.001
*Concrete Compression Damage, tension recovery=1.
0.0, 0.0
0.1, 0.001
*Concrete Tension Stiffening
3.0e6, 0.0
2.0e6, 0.001
*Concrete Tension Damage, compression recovery=0.75
0.0, 0.0
0.1, 0.001
"""


class AbaqusCDPTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.inp = self.work / "case.inp"
        self.inp.write_text(VALID_INPUT, encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def _validated(self):
        materials = CDP.parse_materials(self.inp)
        material = CDP.choose_material(materials, None)
        return CDP.validate_material(material)

    def test_extracts_scalars_tables_recovery_and_plastic_strain(self):
        validated = self._validated()
        self.assertEqual(validated["material"], "CONCRETE")
        self.assertEqual(validated["elastic"]["youngs_modulus_pa"], 3.04e10)
        self.assertEqual(validated["plasticity"]["eccentricity"], 0.1)
        self.assertEqual(validated["plasticity"]["kc"], 0.667)
        self.assertEqual(validated["plasticity"]["viscosity_seconds"], 0.0005)
        self.assertEqual(validated["recovery"]["tension_recovery"], 1.0)
        self.assertEqual(validated["recovery"]["compression_recovery"], 0.75)
        self.assertEqual(len(validated["tables"]["compression_hardening"]), 2)
        self.assertGreater(
            validated["derived"]["compression"][1]["equivalent_plastic_strain"], 0
        )

    def test_piecewise_linear_interpolation_and_constant_endpoints(self):
        table = [(0.0, 0.0), (0.5, 1.0), (1.0, 2.0)]
        self.assertEqual(CDP._interpolate(table, -1.0), 0.0)
        self.assertEqual(CDP._interpolate(table, 0.5), 0.25)
        self.assertEqual(CDP._interpolate(table, 1.0), 0.5)
        self.assertEqual(CDP._interpolate(table, 3.0), 1.0)

    def test_writes_traceable_bundle_and_prototype_input_fragment(self):
        validated = self._validated()
        out_dir = self.work / "bundle"
        manifest = CDP.write_bundle(self.inp, out_dir, validated)

        self.assertEqual(manifest["schema_version"], "cdp-mapping-manifest/0.1")
        self.assertEqual(
            manifest["table_point_counts"],
            {
                "compression_hardening": 2,
                "compression_damage": 2,
                "tension_stiffening": 2,
                "tension_damage": 2,
            },
        )
        self.assertEqual(len(manifest["source"]["sha256"]), 64)
        self.assertEqual(len(manifest["source_files"]), 1)
        self.assertEqual(manifest["provenance"]["plasticity"]["keyword_line"], 4)
        self.assertTrue((out_dir / "compression_hardening.csv").is_file())
        fragment = (out_dir / "abaqus_cdp_material.i").read_text(encoding="utf-8")
        self.assertIn("type = AbaqusCDPStressUpdate", fragment)
        loaded = json.loads(
            (out_dir / "cdp-mapping-manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(loaded["source"]["sha256"], manifest["source"]["sha256"])

    def test_rejects_negative_derived_plastic_strain(self):
        invalid = VALID_INPUT.replace(
            "0.1, 0.001\n*Concrete Tension Stiffening",
            "0.9, 0.001\n*Concrete Tension Stiffening",
            1,
        )
        self.inp.write_text(invalid, encoding="utf-8")
        material = CDP.choose_material(CDP.parse_materials(self.inp), None)
        with self.assertRaisesRegex(CDP.CDPInputError, "负等效塑性应变"):
            CDP.validate_material(material)

    def test_rejects_rate_or_field_dependent_extra_columns(self):
        invalid = VALID_INPUT.replace("2.2e7, 0.001", "2.2e7, 0.001, 20")
        self.inp.write_text(invalid, encoding="utf-8")
        material = CDP.choose_material(CDP.parse_materials(self.inp), None)
        with self.assertRaisesRegex(CDP.CDPInputError, "仅支持 2 列"):
            CDP.validate_material(material)

    def test_rejects_displacement_based_tension_table_in_p0(self):
        invalid = VALID_INPUT.replace(
            "*Concrete Tension Stiffening",
            "*Concrete Tension Stiffening, type=DISPLACEMENT",
        )
        self.inp.write_text(invalid, encoding="utf-8")
        material = CDP.choose_material(CDP.parse_materials(self.inp), None)
        with self.assertRaisesRegex(CDP.CDPInputError, "TYPE=STRAIN"):
            CDP.validate_material(material)

    def test_rejects_kc_at_or_below_official_lower_bound(self):
        for kc in (0.5, 0.25):
            with self.subTest(kc=kc):
                invalid = VALID_INPUT.replace(
                    "36.31, 0.1, 1.16, 0.667, 0.0005",
                    f"36.31, 0.1, 1.16, {kc}, 0.0005",
                )
                self.inp.write_text(invalid, encoding="utf-8")
                material = CDP.choose_material(CDP.parse_materials(self.inp), None)
                with self.assertRaisesRegex(CDP.CDPInputError, "Kc"):
                    CDP.validate_material(material)

    def test_rejects_nonzero_first_table_abscissa(self):
        invalid = VALID_INPUT.replace("2.0e7, 0.0", "2.0e7, 0.0001")
        self.inp.write_text(invalid, encoding="utf-8")
        material = CDP.choose_material(CDP.parse_materials(self.inp), None)
        with self.assertRaisesRegex(CDP.CDPInputError, "首行横坐标必须为 0"):
            CDP.validate_material(material)

    def test_rejects_nonzero_first_damage(self):
        invalid = VALID_INPUT.replace(
            "*Concrete Tension Damage, compression recovery=0.75\n0.0, 0.0",
            "*Concrete Tension Damage, compression recovery=0.75\n0.01, 0.0",
        )
        self.inp.write_text(invalid, encoding="utf-8")
        material = CDP.choose_material(CDP.parse_materials(self.inp), None)
        with self.assertRaisesRegex(CDP.CDPInputError, "首行损伤必须为 0"):
            CDP.validate_material(material)

    def test_requires_material_choice_for_multiple_cdp_materials(self):
        second = VALID_INPUT.replace("name=CONCRETE", "name=CONCRETE-2")
        self.inp.write_text(VALID_INPUT + second, encoding="utf-8")
        with self.assertRaisesRegex(CDP.CDPInputError, "--material"):
            CDP.choose_material(CDP.parse_materials(self.inp), None)

    def test_include_cycle_is_rejected(self):
        other = self.work / "other.inp"
        self.inp.write_text("*Include, input=other.inp\n", encoding="utf-8")
        other.write_text("*Include, input=case.inp\n", encoding="utf-8")
        with self.assertRaisesRegex(CDP.CDPInputError, "循环引用"):
            CDP.parse_materials(self.inp)


if __name__ == "__main__":
    unittest.main()
