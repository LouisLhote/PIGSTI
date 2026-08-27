#!/usr/bin/env python3
"""Validate config and samples before running the Snakemake workflow."""

from __future__ import annotations

import argparse
import csv
import os
import sys

import yaml


REQUIRED_SPREADSHEET_COLS = {"Krakenuniq name", "Hops name", "bwa index"}


def _bowtie2_index_ok(prefix: str) -> bool:
    """True if prefix looks like a built Bowtie2 index (.1.bt2 or .1.bt2l)."""
    p = str(prefix).strip()
    if not p:
        return False
    for suf in (".1.bt2", ".1.bt2l"):
        if os.path.isfile(p + suf):
            return True
    return False




def _dir_ok(path: str) -> bool:
    p = str(path).strip()
    return bool(p) and os.path.isdir(p)


def _file_ok(path: str) -> bool:
    p = str(path).strip()
    return bool(p) and os.path.isfile(p)


def _check_index_dict(cfg: dict, key: str, errors: list[str], *, bowtie2: bool) -> None:
    block = cfg.get(key)
    if not isinstance(block, dict) or not block:
        errors.append(f"config key '{key}' must be a non-empty mapping (species -> path)")
        return
    for species, path in sorted(block.items()):
        p = str(path).strip()
        if not p:
            errors.append(f"{key}[{species}] is empty")
            continue
        if bowtie2:
            if not _bowtie2_index_ok(p) and not _file_ok(p):
                errors.append(
                    f"{key}[{species}] Bowtie2 index not found "
                    f"(expected {p}.1.bt2 or {p}.1.bt2l): {p}"
                )
        else:
            if not (_file_ok(p) or os.path.isdir(p)):
                errors.append(f"{key}[{species}] reference not found: {p}")


def _check_bowtie2_with_fasta_fallback(
    cfg: dict,
    bowtie2_key: str,
    fasta_key: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    """Bowtie2 prefixes must exist or a matching FASTA can be indexed at runtime."""
    block = cfg.get(bowtie2_key)
    if not isinstance(block, dict) or not block:
        errors.append(f"config key '{bowtie2_key}' must be a non-empty mapping (species -> path)")
        return
    fasta_block = cfg.get(fasta_key)
    if not isinstance(fasta_block, dict):
        fasta_block = {}
    for species, prefix in sorted(block.items()):
        p = str(prefix).strip()
        if not p:
            errors.append(f"{bowtie2_key}[{species}] is empty")
            continue
        if _bowtie2_index_ok(p):
            continue
        fasta = str(fasta_block.get(species, "")).strip()
        if _file_ok(fasta):
            warnings.append(
                f"{bowtie2_key}[{species}]: Bowtie2 index missing at {p}; "
                f"pipeline will build from {fasta_key} FASTA: {fasta}"
            )
            continue
        if fasta:
            base = fasta.rsplit(".", 1)[0]
            if _bowtie2_index_ok(base):
                warnings.append(
                    f"{bowtie2_key}[{species}]: prefix {p} missing; "
                    f"index found at FASTA-derived prefix {base}"
                )
                continue
        errors.append(
            f"{bowtie2_key}[{species}] Bowtie2 index not found "
            f"(expected {p}.1.bt2 or {p}.1.bt2l) and no {fasta_key} FASTA: "
            f"{fasta or '(missing)'}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--spreadsheet", default="config/Pathogen_spreadsheet.csv")
    parser.add_argument("--output", default="results/workflow/.pigsti_validation_ok")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []

    if not os.path.isfile(args.config):
        errors.append(f"Missing config: {args.config}")
        cfg = {}
    else:
        with open(args.config, encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}

    if not os.path.isfile(args.samples):
        errors.append(f"Missing samples TSV: {args.samples}")
    else:
        with open(args.samples, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames or "sample" not in [h.strip() for h in reader.fieldnames]:
                errors.append("samples.tsv must have a 'sample' column")
            rows = list(reader)
        if not rows:
            errors.append("samples.tsv has no data rows")
        pcr_ids: set[str] = set()
        for row in rows:
            sample = (row.get("sample") or "").strip()
            pcr = (row.get("pcr") or sample).strip()
            if not sample:
                continue
            if pcr in pcr_ids:
                errors.append(f"Duplicate PCR/library id: {pcr}")
            pcr_ids.add(pcr)
            for mate in ("r1", "r2"):
                path = (row.get(mate) or "").strip()
                if path and not os.path.isfile(path):
                    errors.append(f"Missing FASTQ for {pcr} ({mate}): {path}")

    if not os.path.isfile(args.spreadsheet):
        errors.append(f"Missing pathogen spreadsheet: {args.spreadsheet}")
    else:
        import pandas as pd

        df = pd.read_csv(args.spreadsheet)
        missing = REQUIRED_SPREADSHEET_COLS - set(df.columns.str.strip())
        if missing:
            errors.append(f"Pathogen_spreadsheet.csv missing columns: {sorted(missing)}")

    kraken_db = str(cfg.get("kraken_db", "")).strip()
    if not kraken_db:
        errors.append("config key 'kraken_db' is empty")
    elif not _dir_ok(kraken_db):
        errors.append(
            f"config key 'kraken_db' must be an existing KrakenUniq database directory: {kraken_db}"
        )

    host_index = str(cfg.get("host_index", "")).strip()
    if not host_index:
        errors.append("config key 'host_index' is empty (Bowtie2 chimera index PREFIX for bowtie2_unaligned)")
    elif not _bowtie2_index_ok(host_index):
        errors.append(
            "config key 'host_index' — Bowtie2 index not found. "
            f"Expect files {host_index}.1.bt2 or {host_index}.1.bt2l on disk. "
            "Example: .../multi/chimera.1.bt2l -> host_index: .../multi/chimera "
            "(there is no directory named chimera; only the index prefix)."
        )

    host_aligner = str(cfg.get("host_aligner", "bwa")).strip().lower()
    if host_aligner == "bowtie2":
        _check_bowtie2_with_fasta_fallback(
            cfg, "bowtie2_indices", "bwa_indices", errors, warnings
        )
        _check_bowtie2_with_fasta_fallback(
            cfg, "bowtie2_mtDNA_indices", "mtDNA_indices", errors, warnings
        )
    else:
        _check_index_dict(cfg, "bwa_indices", errors, bowtie2=False)
        _check_index_dict(cfg, "mtDNA_indices", errors, bowtie2=False)

    if cfg.get("enable_sexing", True):
        sexing_candidates = (
            os.path.join("scripts", "sexing", "sexing_residual_method.R"),
            os.path.join("scripts", "sexing_residual_method.R"),
        )
        if not any(os.path.isfile(p) for p in sexing_candidates):
            errors.append(
                "enable_sexing=true but missing sexing_residual_method.R "
                f"(tried: {', '.join(sexing_candidates)})"
            )

    if cfg.get("enable_hops"):
        if not cfg.get("hops_malt_index"):
            warnings.append("enable_hops=true but hops_malt_index not set")
        if cfg.get("hops_parallel"):
            jobs = int(cfg.get("hops_parallel_jobs", 2) or 2)
            if jobs < 1:
                errors.append("hops_parallel_jobs must be >= 1 when hops_parallel=true")
            per_job = cfg.get("hops_threads_per_job")
            if per_job is not None and int(per_job) < 1:
                errors.append("hops_threads_per_job must be >= 1 when set")

    if cfg.get("enable_decom"):
        src = cfg.get("decOM_sources", "")
        if not src or not os.path.isdir(str(src)):
            errors.append(f"enable_decom=true but decOM_sources missing or not a directory: {src}")

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        return 1

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("ok\n")
    print(f"Validation passed. Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
