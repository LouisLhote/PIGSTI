#!/usr/bin/env python3
"""
Build (bio_sample, pathogen) pairs for pathogen BWA/Bowtie2 mapping.

Used by Snakemake checkpoint ``generate_pathogen_targets``. When ``--strict-cohort``
is set, every biological sample must appear in ``--evalue`` inputs (Snakemake
guarantees those files came from completed jobs, not stale partial runs).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pigsti_naming import safe_pathogen_name as safe_name
from pigsti_paths import sample_from_evalue_pathogen_csv as sample_from_escore_pathogen_csv


def load_detection_cfg(config_path: str) -> dict:
    defaults = {
        "map_all_escore_pathogens": True,
        "use_evalue_for_detection": True,
        "escore_threshold": 5,
        "reads_threshold": 50,
        "guellil_evalue_threshold": 0.001,
    }
    cfg = {}
    if config_path and os.path.isfile(config_path):
        with open(config_path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
    det = cfg.get("pathogen_detection_criteria") or {}
    if isinstance(det, dict):
        defaults.update(det)
    legacy = os.path.join(os.path.dirname(config_path or "."), "pathogen_detection_criteria.yaml")
    if os.path.isfile(legacy):
        with open(legacy, "r", encoding="utf-8") as f:
            leg = yaml.safe_load(f) or {}
        if isinstance(leg, dict):
            for k, v in leg.items():
                defaults.setdefault(k, v)
    return defaults


def read_hops_dataframe(hops_path: str | None) -> pd.DataFrame | None:
    if not hops_path or not os.path.isfile(hops_path):
        return None
    hops_df = pd.read_csv(hops_path, sep="\t")
    hops_df.columns = hops_df.columns.str.replace('"', "").str.strip()
    hops_df.rename(columns={"node": "Species"}, inplace=True)
    hops_df["Species"] = hops_df["Species"].str.replace('"', "").str.strip()
    return hops_df


def pairs_for_sample(
    sample: str,
    escore_path: str,
    hops_df: pd.DataFrame | None,
    spreadsheet_df: pd.DataFrame,
    detection_cfg: dict,
) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    escore_set: set[str] = set()
    hops_set: set[str] = set()

    escore_df = pd.read_csv(escore_path)
    escore_df.columns = escore_df.columns.str.strip()
    if "Guellil_et_al_Evalue" in escore_df.columns:
        escore_df["Guellil_et_al_Evalue"] = pd.to_numeric(
            escore_df["Guellil_et_al_Evalue"], errors="coerce"
        )
    if "taxonomy" not in escore_df.columns:
        return pairs

    escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()

    if bool(detection_cfg.get("map_all_escore_pathogens", True)):
        for p in escore_df["taxonomy"].dropna().astype(str):
            p = p.strip()
            if p:
                escore_set.add(p)
    else:
        reads_col = None
        if "reads" in escore_df.columns:
            reads_col = "reads"
        elif "taxReads" in escore_df.columns:
            reads_col = "taxReads"
        elif "# of reads" in escore_df.columns:
            reads_col = "# of reads"

        use_evalue = bool(detection_cfg.get("use_evalue_for_detection", False))
        escore_min_default = detection_cfg.get("escore_threshold", 5)
        reads_min_default = detection_cfg.get("reads_threshold", 50)
        evalue_max_default = detection_cfg.get("guellil_evalue_threshold", 0.001)

        for _, row in escore_df.iterrows():
            pathogen_name = row.get("taxonomy")
            if not isinstance(pathogen_name, str):
                continue
            pathogen_name = pathogen_name.strip()

            row_cfg = spreadsheet_df[spreadsheet_df["Krakenuniq name"] == pathogen_name]
            min_escore = escore_min_default
            min_reads = reads_min_default
            if not row_cfg.empty:
                if "min_escore" in row_cfg.columns:
                    val = row_cfg.iloc[0]["min_escore"]
                    if pd.notna(val):
                        try:
                            min_escore = float(val)
                        except Exception:
                            pass
                if "min_reads" in row_cfg.columns:
                    val = row_cfg.iloc[0]["min_reads"]
                    if pd.notna(val):
                        try:
                            min_reads = int(val)
                        except Exception:
                            pass

            max_evalue = evalue_max_default
            if not row_cfg.empty:
                for colname in [
                    "Guellil_et_al_Evalue_threshold",
                    "max_evalue",
                    "evalue_threshold",
                ]:
                    if colname in row_cfg.columns:
                        val = row_cfg.iloc[0][colname]
                        if pd.notna(val):
                            try:
                                max_evalue = float(val)
                            except Exception:
                                pass
                        break

            reads_ok = True
            if reads_col is not None and reads_col in escore_df.columns:
                try:
                    reads_val = int(row.get(reads_col, 0))
                except Exception:
                    reads_val = 0
                reads_ok = reads_val >= min_reads

            keep_pathogen = False
            if use_evalue:
                if "Guellil_et_al_Evalue" in escore_df.columns:
                    try:
                        evalue_raw = row.get("Guellil_et_al_Evalue")
                        if pd.notna(evalue_raw):
                            evalue_val = float(evalue_raw)
                            keep_pathogen = evalue_val > max_evalue and reads_ok
                    except (ValueError, TypeError, KeyError):
                        keep_pathogen = False
            else:
                try:
                    escore_raw = row.get("Escore")
                    if pd.notna(escore_raw):
                        escore_val = float(escore_raw)
                        keep_pathogen = escore_val >= min_escore and reads_ok
                except (ValueError, TypeError, KeyError):
                    keep_pathogen = False

            if keep_pathogen:
                escore_set.add(pathogen_name)

    if hops_df is not None:
        sample_col = f"{sample}_unaligned.rma6"
        if sample_col in hops_df.columns:
            hops_species = hops_df.loc[hops_df[sample_col] > 1, "Species"]
            for hops_name in hops_species:
                row = spreadsheet_df[
                    spreadsheet_df["Hops name"].str.replace(" ", "_") == hops_name
                ]
                if not row.empty:
                    hops_set.add(row.iloc[0]["Krakenuniq name"])

    union_pathogens = escore_set | hops_set
    for pathogen in sorted(union_pathogens):
        if pathogen in spreadsheet_df["Krakenuniq name"].values:
            pairs.append((sample, pathogen))
    return pairs


def compute_pairs(
    evalue_paths: list[str],
    bio_samples: list[str],
    spreadsheet_df: pd.DataFrame,
    detection_cfg: dict,
    hops_path: str | None = None,
    strict_cohort: bool = True,
) -> list[tuple[str, str]]:
    escore_by_sample: dict[str, str] = {}
    for p in evalue_paths:
        escore_by_sample[sample_from_escore_pathogen_csv(p)] = p

    if strict_cohort:
        missing = sorted(set(bio_samples) - set(escore_by_sample.keys()))
        if missing:
            raise ValueError(
                "Pathogen mapping cohort not ready: missing E-value pathogen CSV for "
                f"biological samples {missing}. Ensure all KrakenUniq/E-score jobs finished "
                "(rule metagenomics_screening_cohort_ready) before generate_pathogen_targets."
            )

    hops_df = read_hops_dataframe(hops_path)
    all_pairs: list[tuple[str, str]] = []
    for sample in bio_samples:
        if sample not in escore_by_sample:
            continue
        all_pairs.extend(
            pairs_for_sample(
                sample,
                escore_by_sample[sample],
                hops_df,
                spreadsheet_df,
                detection_cfg,
            )
        )
    return all_pairs


def write_targets_and_manifest(
    pairs: list[tuple[str, str]],
    bio_samples: list[str],
    evalue_paths: list[str],
    hops_path: str | None,
    targets_path: str,
    manifest_path: str,
) -> None:
    os.makedirs(os.path.dirname(targets_path) or ".", exist_ok=True)
    with open(targets_path, "w", encoding="utf-8") as f:
        for sample, pathogen in pairs:
            ps = safe_name(pathogen)
            f.write(f"results/pathogen/{sample}/pathogen_mapping/{sample}_{ps}.dedup.bam\n")

    manifest = {
        "status": "ok",
        "cohort_samples": sorted(bio_samples),
        "n_cohort_samples": len(bio_samples),
        "n_mapping_pairs": len(pairs),
        "pairs": [[s, p] for s, p in pairs],
        "escore_inputs": sorted(evalue_paths),
        "hops_input": hops_path,
        "targets_file": targets_path,
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--spreadsheet", default="config/Pathogen_spreadsheet.csv")
    parser.add_argument("--bio-samples", nargs="+", required=True)
    parser.add_argument(
        "--evalue",
        "--escore",
        nargs="+",
        required=True,
        dest="evalue",
        metavar="PATH",
        help="E-value pathogen CSV paths (results/pathogen/{bio}/evalue/pathogen/...)",
    )
    parser.add_argument("--hops", default="", help="HOPS heatmap TSV (optional)")
    parser.add_argument("--output-targets", required=True)
    parser.add_argument("--output-manifest", required=True)
    parser.add_argument(
        "--strict-cohort",
        action="store_true",
        default=True,
        help="Require every --bio-samples to have a matching --evalue file",
    )
    parser.add_argument(
        "--no-strict-cohort",
        action="store_false",
        dest="strict_cohort",
    )
    args = parser.parse_args()

    detection_cfg = load_detection_cfg(args.config)
    spreadsheet_df = pd.read_csv(args.spreadsheet)
    hops_path = args.hops.strip() or None

    pairs = compute_pairs(
        list(args.evalue),
        list(args.bio_samples),
        spreadsheet_df,
        detection_cfg,
        hops_path=hops_path,
        strict_cohort=args.strict_cohort,
    )
    write_targets_and_manifest(
        pairs,
        list(args.bio_samples),
        list(args.evalue),
        hops_path,
        args.output_targets,
        args.output_manifest,
    )
    print(
        f"Wrote {len(pairs)} mapping targets for {len(args.bio_samples)} samples "
        f"→ {args.output_targets}"
    )


if __name__ == "__main__":
    main()
