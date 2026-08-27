#!/usr/bin/env python3
"""
Build a browse-friendly view of PIGSTI results:

1. results/final/output_catalog.tsv — every key artefact with analysis / level / path
2. results/browse/ — optional symlinks grouped by analysis (no Snakemake path change)

Run from repo root after a pipeline run, or via Snakemake rule build_results_catalog.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

BY_TYPE_ROOT = "results/browse"
CATALOG_PATH = "results/final/output_catalog.tsv"


def load_samples(samples_tsv: str) -> tuple[list[str], list[str], dict[str, list[str]]]:
    st = pd.read_csv(samples_tsv, sep="\t")
    if "pcr" not in st.columns or "sample" not in st.columns:
        raise ValueError(f"samples.tsv must have 'pcr' and 'sample' columns: {samples_tsv}")
    pcrs = st["pcr"].astype(str).tolist()
    bio = sorted(st["sample"].astype(str).unique().tolist())
    pcr_to_bio = dict(zip(st["pcr"].astype(str), st["sample"].astype(str)))
    bio_to_pcrs: dict[str, list[str]] = {}
    for pcr, bio_id in pcr_to_bio.items():
        bio_to_pcrs.setdefault(bio_id, []).append(pcr)
    return pcrs, bio, bio_to_pcrs


def _add(rows: list[dict], analysis: str, level: str, sample_id: str, path: str) -> None:
    if path and os.path.isfile(path):
        rows.append(
            {
                "analysis": analysis,
                "level": level,
                "sample_id": sample_id,
                "path": path.replace("\\", "/"),
                "size_bytes": os.path.getsize(path),
            }
        )


def _link(target: Path, link_path: Path, use_symlinks: bool) -> None:
    link_path.parent.mkdir(parents=True, exist_ok=True)
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    if not target.is_file():
        return
    if use_symlinks:
        try:
            link_path.symlink_to(target.resolve())
            return
        except OSError:
            pass
    # Windows / OneDrive fallback: write a pointer file
    link_path.write_text(str(target.resolve()), encoding="utf-8")


def collect_catalog(
    pcrs: list[str],
    bio_samples: list[str],
    bio_to_pcrs: dict[str, list[str]],
    enable_hops: bool,
) -> list[dict]:
    rows: list[dict] = []

    for pcr in pcrs:
        _add(rows, "adapter_removal", "pcr", pcr, f"results/{pcr}/adapter_removal/{pcr}.collapsed.gz")
        _add(rows, "prinseq", "pcr", pcr, f"results/{pcr}/prinseq/{pcr}-passed.fq.gz")
        _add(rows, "fastq_screen", "pcr", pcr, f"results/{pcr}/fastq_screen/{pcr}_best_species.txt")
        _add(rows, "host_mapping", "pcr", pcr, f"results/libraries/{pcr}/host_mapping/{pcr}.dedup.bam")
        _add(rows, "mtdna_mapping", "pcr", pcr, f"results/libraries/{pcr}/mtdna_mapping/{pcr}.dedup.bam")
        _add(rows, "unaligned_fastq", "pcr", pcr, f"results/{pcr}/unaligned_fastq/{pcr}_unaligned.fastq.gz")

    for bio in bio_samples:
        _add(
            rows,
            "krakenuniq",
            "bio",
            bio,
            f"results/metagenomics/krakenuniq/{bio}/{bio}_kraken-report.txt",
        )
        _add(
            rows,
            "krakenuniq",
            "bio",
            bio,
            f"results/metagenomics/krakenuniq/{bio}/{bio}_output.txt",
        )
        _add(rows, "evalue_pathogen", "pathogen", bio, f"results/pathogen/{bio}/evalue/pathogen/{bio}_pathogen.csv")
        _add(rows, "evalue_genus", "pathogen", bio, f"results/pathogen/{bio}/evalue/genus/{bio}_genus.csv")
        _add(rows, "evalue_species", "pathogen", bio, f"results/pathogen/{bio}/evalue/species/{bio}_species.csv")
        _add(rows, "qualimap_host", "sample", bio, f"results/samples/{bio}/qualimap/genome_results.txt")
        _add(rows, "qualimap_mtdna", "sample", bio, f"results/samples/{bio}/qualimap_mtdna/genome_results.txt")
        _add(rows, "sexing", "sample", bio, f"results/samples/{bio}/sexing/{bio}_sexing.pdf")
        _add(rows, "sexing", "sample", bio, f"results/samples/{bio}/sexing/{bio}_sexing.pdf")
        _add(rows, "sexing", "sample", bio, f"results/samples/{bio}/sexing/{bio}_sexing.tsv")
        _add(
            rows,
            "sample_unaligned_merged",
            "bio",
            bio,
            f"results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz",
        )
        _add(
            rows,
            "sample_host_merged_bam",
            "host",
            bio,
            f"results/host/host_mapping/{bio}.dedup.merged.bam",
        )
        _add(
            rows,
            "sample_mtdna_merged_bam",
            "host",
            bio,
            f"results/host/mtdna_mapping/{bio}.dedup.merged.bam",
        )
        _add(rows, "sample_report", "bio", bio, f"results/pathogen/{bio}/summary/{bio}_sample_report.pdf")
        if enable_hops:
            _add(rows, "comparison", "bio", bio, f"results/pathogen/{bio}/comparison/{bio}_comparison.tsv")

    for bio in bio_samples:
        bwa_dir = Path(f"results/pathogen/{bio}/pathogen_mapping")
        if bwa_dir.is_dir():
            for bam in sorted(bwa_dir.glob(f"{bio}_*.dedup.bam")):
                _add(rows, "pathogen_mapping", "bio", bio, str(bam).replace("\\", "/"))

    final_files = [
        ("final_pathogen_summary", "final", "all", "results/final/pathogen_summary_all_samples.xlsx"),
        ("final_host_mtdna", "final", "all", "results/final/host_mtdna_summary_all_samples.xlsx"),
        ("final_comprehensive", "final", "all", "results/final/comprehensive_summary_all_samples.xlsx"),
        ("run_manifest", "final", "all", "results/final/run_manifest.json"),
        ("pathogen_targets", "workflow", "all", "results/workflow/pathogen_targets.txt"),
        ("pathogen_targets_manifest", "workflow", "all", "results/workflow/pathogen_targets.manifest.json"),
        ("kraken_abundance", "metagenomics", "all", "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv"),
    ]
    if enable_hops:
        final_files.append(
            ("hops_heatmap", "metagenomics", "all", "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv")
        )
    for analysis, level, sid, path in final_files:
        _add(rows, analysis, level, sid, path)

    return rows


def build_symlinks(rows: list[dict], use_symlinks: bool, by_type_root: str) -> int:
    n = 0
    for row in rows:
        src = Path(row["path"])
        if not src.is_file():
            continue
        analysis = row["analysis"]
        sid = row["sample_id"]
        dest = Path(by_type_root) / analysis / sid / src.name
        _link(src, dest, use_symlinks)
        n += 1
    return n


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--catalog", default=CATALOG_PATH)
    parser.add_argument("--by-type-root", default=BY_TYPE_ROOT)
    parser.add_argument("--symlinks", action="store_true", help="Create results/browse symlinks")
    parser.add_argument("--no-symlinks", action="store_true")
    args = parser.parse_args()

    cfg = {}
    if os.path.isfile(args.config):
        with open(args.config, encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}

    samples_tsv = cfg.get("samples", "config/samples.tsv")
    enable_hops = bool(cfg.get("enable_hops", False))
    pcrs, bio, bio_to_pcrs = load_samples(samples_tsv)

    rows = collect_catalog(pcrs, bio, bio_to_pcrs, enable_hops)
    if not rows:
        print("No result files found to catalog.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows).sort_values(["analysis", "level", "sample_id", "path"])
    os.makedirs(os.path.dirname(args.catalog) or ".", exist_ok=True)
    df.to_csv(args.catalog, sep="\t", index=False)
    print(f"Wrote {len(df)} rows → {args.catalog}")

    do_links = args.symlinks and not args.no_symlinks
    if do_links:
        n = build_symlinks(rows, use_symlinks=True, by_type_root=args.by_type_root)
        print(f"Created {n} entries under {args.by_type_root}/")


if __name__ == "__main__":
    main()
