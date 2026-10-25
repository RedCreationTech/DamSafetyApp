"""Protocol/static tests only. These do not evaluate CDP mechanics or run MOOSE."""
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

from prepare_point_replay import TABLES, parse_packet, prepare


def encoded(*values):
    return " ".join(struct.pack(">d", x).hex() for x in values)


def synthetic_packet():
    lines = ["AQ_CDP_POINT_V1", "context 19 2 4294967295 0 0 101 " + encoded(1.2113125, .0025, 1., 2., 3.),
             'material "synthetic_protocol_fixture"',
             "local " + encoded(30e9, .2, 36.31, .1, 1.16, .6667) + " 40 " + encoded(1e-9, 1e-7, 1e-6) + " 1",
             "state " + encoded(0., 1., .05, 1e-12), "substeps 256 " + encoded(2.5e-5, 1e-8)]
    for field, count in (("old_strain", 6), ("new_strain", 6), ("old_backbone", 6),
                          ("old_kappa", 2), ("old_viscous", 6), ("old_damage", 2)):
        values = [(i + 1) * 1e-6 for i in range(count)]
        if field == "old_strain":
            values[0] = -0.
        lines.append(field + " " + encoded(*values))
    for table in TABLES:
        lines += [table + " 2", encoded(0., 0.), encoded(1., .1)]
    lines += ['expected_error "synthetic transport error, not a real CDP failure"', "END"]
    return ("\n".join(lines) + "\n").encode()


class PacketProtocolTests(unittest.TestCase):
    def test_complete_state_and_signed_zero_are_preserved(self):
        p = parse_packet(synthetic_packet())
        self.assertEqual(p["step"], 101)
        self.assertEqual(sum(len(p[x]) for x in ("old_backbone", "old_kappa", "old_viscous", "old_damage")), 16)
        self.assertEqual(math.copysign(1., p["old_strain"][0]), -1.)
        self.assertEqual(p["local"]["residual_tolerance"], 1e-9)
        self.assertEqual(p["state_parameters"][2], .05)
        for field, count in (("old_backbone", 6), ("old_kappa", 2),
                             ("old_viscous", 6), ("old_damage", 2)):
            self.assertEqual(p[field], [(i + 1) * 1e-6 for i in range(count)])

    def test_truncated_and_version_mismatch_rejected(self):
        for data in (synthetic_packet()[:100], synthetic_packet().replace(b"V1", b"V2", 1)):
            with self.assertRaises(ValueError):
                parse_packet(data)

    def test_nonfinite_diagnostic_input_requires_manual_review(self):
        data = synthetic_packet().replace(encoded(30e9).encode(), encoded(float("nan")).encode(), 1)
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            parse_packet(data)

    def test_missing_exception_trailing_data_and_oversize_table_rejected(self):
        for data in (synthetic_packet() + b"extra", synthetic_packet().replace(b'"synthetic transport error, not a real CDP failure"', b'""'),
                     synthetic_packet().replace(b"compression_hardening 2", b"compression_hardening 100001")):
            with self.assertRaises(ValueError):
                parse_packet(data)

    def test_bad_binary64_and_unsigned_identifiers_rejected(self):
        for data in (synthetic_packet().replace(b"context 19", b"context -1"),
                     synthetic_packet().replace(encoded(.0025).encode(), b"not_binary64", 1)):
            with self.assertRaises(ValueError):
                parse_packet(data)

    def test_package_spaces_checksums_and_identity_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            packet = root / "point with spaces.txt"
            packet.write_bytes(synthetic_packet())
            output = root / "prepared replay"
            m = prepare(packet, output, "1" * 40, "2" * 64, "job_protocol_fixture")
            self.assertFalse(m["capture_replay_performed"])
            self.assertFalse(m["source_binary_identity_independently_verified"])
            self.assertEqual((output / "capture.txt").read_bytes(), packet.read_bytes())
            import hashlib
            for name, digest in json.loads((output / "checksums.json").read_text()).items():
                self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
            with self.assertRaises(FileExistsError):
                prepare(packet, output, "1" * 40, "2" * 64, "job_protocol_fixture")

    def test_invalid_declared_identity_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            packet = root / "capture.txt"
            packet.write_bytes(synthetic_packet())
            output = root / "out"
            for source, binary, job in (("branch", "2" * 64, "job_x"), ("1" * 40, "unknown", "job_x"),
                                         ("1" * 40, "2" * 64, "not_a_job")):
                with self.assertRaises(ValueError):
                    prepare(packet, output, source, binary, job)
                self.assertFalse(output.exists())


class CandidateDefaultsStaticTests(unittest.TestCase):
    def test_defaults_original_full_integrator_and_exception_preserved(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "src/materials/AbaqusCDPStressUpdate.C").read_text()
        self.assertIn('"enable_exact_result_cache", false', source)
        self.assertIn('"observe_exact_result_cache", false', source)
        self.assertIn('"failure_point_capture_prefix", ""', source)
        self.assertIn('if(!_exact_cache)\n    return _substep_integrator.integrateLinearized', source)
        self.assertIn('throw MooseException(error.what());', source)
        self.assertIn('StressUpdateBase::meshChanged();', source)
        self.assertNotIn('_substep_integrator.integrate(', source)


if __name__ == "__main__":
    unittest.main()
