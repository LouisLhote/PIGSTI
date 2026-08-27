#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: str) -> dict:
    return {"sha256": sha256_file(path), "size_bytes": os.path.getsize(path)}


def git_revision() -> dict | None:
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return {
            "commit": commit,
            "branch": branch,
            "dirty": bool(dirty),
        }
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def load_effective_config(path: str) -> dict:
    try:
        import yaml
    except ImportError:
        return {}
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    # Keys that affect scientific outputs (paths omitted — use file_hashes for configs).
    keys = [
        "pathogen_mapping_mode",
        "pathogen_aligner",
        "host_aligner",
        "enable_hops",
        "enable_decom",
        "enable_diamond",
        "strict_inputs",
        "dedup_tool",
        "merge_tool",
        "defaults",
        "pathogen_detection_criteria",
        "edit_distance_damage_split",
        "edit_distance_damage_window_size",
        "cleanup_intermediates",
    ]
    return {k: cfg[k] for k in keys if k in cfg}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/run_manifest.json")
    parser.add_argument(
        "--snakefile",
        default="Snakefile",
        help="Path to Snakefile (default: Snakefile in CWD).",
    )
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Main config YAML path.",
    )
    args = parser.parse_args()

    files_to_record = [
        args.snakefile,
        "PIGSTI_snakemake.yaml",
        args.config,
        "config/samples.tsv",
        "config/Pathogen_spreadsheet.csv",
        "config/pathogen_detection_criteria.yaml",
    ]

    env_dir = "workflow/envs"
    try:
        for name in sorted(os.listdir(env_dir)):
            if name.endswith(".yaml") or name.endswith(".yml"):
                files_to_record.append(os.path.join(env_dir, name))
    except FileNotFoundError:
        pass

    seen = set()
    files_to_record_dedup = []
    for p in files_to_record:
        if p not in seen:
            seen.add(p)
            files_to_record_dedup.append(p)

    inputs_hashes = {}
    for p in files_to_record_dedup:
        if os.path.exists(p):
            inputs_hashes[p] = file_record(p)

    manifest = {
        "pipeline": "PIGSTI",
        "manifest_version": "1.1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "snakemake_command": os.environ.get("SNAKEMAKE_COMMAND"),
        "git": git_revision(),
        "effective_config": load_effective_config(args.config),
        "file_hashes": inputs_hashes,
        "checkpoint_outputs": {
            "pathogen_targets": file_record("results/workflow/pathogen_targets.txt")
            if os.path.exists("results/workflow/pathogen_targets.txt")
            else None
        },
        "reproducibility_notes": {
            "scientific_outputs": (
                "Re-run with the same commit, config file hashes, reference DBs, "
                "and FASTQ inputs. Compare sha256/size_bytes in file_hashes and "
                "checkpoint_outputs to verify a match."
            ),
            "non_scientific_fields": (
                "created_at_utc and monitoring HTML timestamps may differ between runs."
            ),
        },
    }

    out_path = args.output
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
