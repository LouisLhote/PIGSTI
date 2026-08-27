#!/usr/bin/env python3
"""Validate decOM sink list and per-sample .fof files before running decOM."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _parse_fof_line(line: str) -> tuple[str, list[str]]:
    line = line.strip()
    if not line or ":" not in line:
        raise ValueError(f"Invalid .fof line (expected 'sample: path'): {line!r}")
    name, rest = line.split(":", 1)
    name = name.strip()
    paths = [p.strip() for p in rest.split(";") if p.strip()]
    if not name or not paths:
        raise ValueError(f"Invalid .fof line: {line!r}")
    return name, paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-sink", required=True, help="p_sink.txt path")
    parser.add_argument("--p-keys-dir", required=True, help="Directory with {sample}.fof files")
    args = parser.parse_args()

    errors: list[str] = []
    p_sink = Path(args.p_sink)
    keys_dir = Path(args.p_keys_dir)

    if not p_sink.is_file():
        errors.append(f"Missing p_sink file: {p_sink}")
    else:
        sinks = [ln.strip() for ln in p_sink.read_text(encoding="utf-8").splitlines() if ln.strip()]
        if not sinks:
            errors.append(f"Empty p_sink file: {p_sink}")

    fof_files = sorted(keys_dir.glob("*.fof"))
    if not fof_files:
        errors.append(f"No .fof files in {keys_dir}")

    sink_names = set(sinks) if p_sink.is_file() else set()
    fof_names = {p.stem for p in fof_files}
    missing_fof = sorted(sink_names - fof_names)
    extra_fof = sorted(fof_names - sink_names)
    if missing_fof:
        errors.append(f"p_keys missing .fof for sinks: {missing_fof}")
    if extra_fof:
        errors.append(f"Extra .fof files without p_sink entry: {extra_fof}")

    for fof in fof_files:
        for line in fof.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                _name, paths = _parse_fof_line(line)
            except ValueError as exc:
                errors.append(f"{fof}: {exc}")
                continue
            for raw in paths:
                path = Path(raw)
                if not path.is_file():
                    errors.append(f"{fof.name}: FASTQ not found: {path}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(f"decOM inputs OK ({len(fof_files)} .fof files, {len(sink_names)} sinks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
