#!/usr/bin/env python3
"""Validate/package an opt-in CDP point capture; never compile or execute a solver."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import struct

TABLES = ("compression_hardening", "compression_damage", "tension_stiffening", "tension_damage")
HEADERS = ("stress_pa,inelastic_strain", "damage_c,inelastic_strain",
           "stress_pa,cracking_strain", "damage_t,cracking_strain")


def parse_packet(data):
    if len(data) > 32 * 1024 * 1024:
        raise ValueError("Capture exceeds 32 MiB bound")
    lexer = shlex.shlex(data.decode("utf-8"), posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    tokens = iter(lexer)

    def take():
        try:
            return next(tokens)
        except StopIteration:
            raise ValueError("Truncated point packet") from None

    def expect(key):
        if take() != key:
            raise ValueError("Missing/incorrect field: " + key)

    def integer(limit=(1 << 64) - 1):
        token = take()
        if not re.fullmatch(r"[0-9]+", token) or int(token) > limit:
            raise ValueError("Invalid integer")
        return int(token)

    def numbers(count):
        result = []
        for _ in range(count):
            token = take()
            if not re.fullmatch(r"[0-9a-fA-F]{16}", token):
                raise ValueError("Invalid binary64 token")
            value = struct.unpack(">d", bytes.fromhex(token))[0]
            if not math.isfinite(value):
                raise ValueError("Non-finite captured input: preserve raw capture, no automatic replay")
            result.append(value)
        return result

    expect("AQ_CDP_POINT_V1")
    expect("context")
    p = dict(element=integer(), qp=integer((1 << 32)-1), side=integer((1 << 32)-1), rank=integer((1 << 32)-1),
             thread=integer((1 << 32)-1), step=integer())
    context = numbers(5)
    p.update(time=context[0], dt=context[1], position=context[2:])
    expect("material")
    p["material"] = take()
    expect("local")
    local = numbers(6)
    iterations = integer((1 << 32)-1)
    local += numbers(3)
    p["local"] = dict(zip(("E", "nu", "dilation", "eccentricity", "fb_ratio", "Kc",
                            "residual_tolerance", "finite_difference_step", "minimum_line_search"), local))
    p["local"].update(maximum_iterations=iterations, use_AD=bool(integer(1)))
    expect("state")
    p["state_parameters"] = numbers(4)
    expect("substeps")
    p["substep_parameters"] = [integer((1 << 32)-1)] + numbers(2)
    for key, count in (("old_strain", 6), ("new_strain", 6), ("old_backbone", 6),
                       ("old_kappa", 2), ("old_viscous", 6), ("old_damage", 2)):
        expect(key)
        p[key] = numbers(count)
    p["tables"] = {}
    for name in TABLES:
        expect(name)
        size = integer(100000)
        if size < 2:
            raise ValueError("Table must contain at least two rows")
        p["tables"][name] = [numbers(2) for _ in range(size)]
    expect("expected_error")
    p["expected_error"] = take()
    if not p["expected_error"]:
        raise ValueError("Missing original exception; no claimed failure replay")
    expect("END")
    if next(tokens, None) is not None:
        raise ValueError("Trailing point packet data")
    return p


def prepare(packet, output, source_commit, binary_sha256, job_id):
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        raise ValueError("Immutable full source commit required")
    if not re.fullmatch(r"[0-9a-f]{64}", binary_sha256):
        raise ValueError("Captured binary SHA256 required")
    if not re.fullmatch(r"job_[A-Za-z0-9_]+", job_id):
        raise ValueError("Formal capture Job ID required")
    packet, output = Path(packet), Path(output)
    original = packet.read_bytes()
    parsed = parse_packet(original)
    output.mkdir(parents=True, exist_ok=False)
    (output / "capture.txt").write_bytes(original)
    tables = output / "tables"
    tables.mkdir()
    for name, header in zip(TABLES, HEADERS):
        text = header + "\n" + "".join("{:.17g},{:.17g}\n".format(*row) for row in parsed["tables"][name])
        (tables / (name + ".csv")).write_text(text, encoding="ascii")
        # Verify exact binary64 values survive canonical CSV, including signed zero.
        for actual, expected in zip(text.splitlines()[1:], parsed["tables"][name]):
            for got, wanted in zip(map(float, actual.split(",")), expected):
                if struct.pack(">d", got) != struct.pack(">d", wanted):
                    raise ValueError("Derived table lost captured input bits")
    (output / "capture-summary.json").write_text(json.dumps(parsed, indent=2, allow_nan=False) + "\n")
    manifest = dict(schema="AQ-CDP-POINT-REPLAY-1", captured_job_id=job_id,
                    declared_capture_source_commit=source_commit,
                    declared_capture_binary_sha256=binary_sha256,
                    source_binary_identity_independently_verified=False,
                    identity_note="Match declarations with capture Job manifest/node SHA before execution.",
                    capture_sha256=hashlib.sha256(original).hexdigest(),
                    original_exception=parsed["expected_error"],
                    capture_replay_performed=False, new_job_submitted=False,
                    command_arguments=["unit/DamSafetyApp-unit-opt", "--gtest_filter=AbaqusCDPPointReplay.ExternalCapturedFailure"],
                    environment=dict(CDP_REPLAY_PACKET="<absolute package>/capture.txt",
                                     CDP_REPLAY_OUTPUT="<new isolated replay output>"),
                    node_gate="Only after original Job terminal/user notice, pinned environment/build and explicit test authorization; no online solver replacement.")
    (output / "replay-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    files = [p for p in output.rglob("*") if p.is_file()]
    (output / "checksums.json").write_text(json.dumps({str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
                                                      for p in sorted(files)}, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--job-id", required=True)
    args = parser.parse_args()
    try:
        result = prepare(args.packet, args.output, args.source_commit, args.binary_sha256, args.job_id)
    except (OSError, ValueError) as error:
        parser.exit(2, str(error) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
