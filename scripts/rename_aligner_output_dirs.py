#!/usr/bin/env python3
"""
LEGACY ONLY — one-off upgrade for old result trees that still use bwa_host/ / bwa_mtdna/.

The current pipeline (Snakefile + pigsti_paths.py) already uses host_mapping/ and
mtdna_mapping/. Normal runs do not need this script; use a fresh results/ directory.

--move-results renames directories under an existing results/ tree from a pre-v3 naming run.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

# Order: longer / more specific tokens first
TEXT_REPLACEMENTS = [
    ("pathogen_bwa_complete", "pathogen_mapping_complete"),
    ("bwa_mtdna", "mtdna_mapping"),
    ("bwa_host", "host_mapping"),
]

GLOBS = [
    "Snakefile",
    "README.md",
    "config/config.example.yaml",
    "scripts/*.py",
    "scripts/*.R",
    "docs/*.md",
]

SKIP = {"rename_aligner_output_dirs.py", "migrate_snakefile_paths_v2.py", "reorganize_results_v3.py"}

RESULT_DIR_MOVES = [
    ("results/libraries", "bwa_host", "host_mapping"),
    ("results/libraries", "bwa_mtdna", "mtdna_mapping"),
    ("results/pools", "bwa_host", "host_mapping"),
    ("results/pools", "bwa_mtdna", "mtdna_mapping"),
    ("results/browse", "bwa_host", "host_mapping"),
    ("results/browse", "bwa_mtdna", "mtdna_mapping"),
    ("results/workflow", "pathogen_bwa_complete.txt", "pathogen_mapping_complete.txt"),
    ("logs", "bwa_host", "host_mapping"),
    ("logs", "bwa_mtdna", "mtdna_mapping"),
]


def migrate_text(text: str) -> str:
    for old, new in TEXT_REPLACEMENTS:
        text = text.replace(old, new)
    return text


def update_repo(root: Path) -> int:
    n = 0
    for pat in GLOBS:
        for path in root.glob(pat):
            if path.name in SKIP:
                continue
            try:
                original = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            updated = migrate_text(original)
            if updated != original:
                path.write_text(updated, encoding="utf-8")
                n += 1
                print(f"updated {path.relative_to(root)}")
    return n


def move_results(root: Path) -> None:
    for parent, old_name, new_name in RESULT_DIR_MOVES:
        base = root / parent
        if not base.is_dir():
            continue
        if old_name.endswith(".txt"):
            for old in base.rglob(old_name):
                new = old.with_name(new_name)
                if old.exists() and not new.exists():
                    old.rename(new)
                    print(f"moved {old} -> {new}")
            continue
        for old_dir in base.rglob(old_name):
            if not old_dir.is_dir():
                continue
            new_dir = old_dir.parent / new_name
            if new_dir.exists():
                print(f"skip (exists): {new_dir}")
                continue
            old_dir.rename(new_dir)
            print(f"moved {old_dir} -> {new_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="PIGSTI repo root (default: cwd)",
    )
    parser.add_argument(
        "--move-results",
        action="store_true",
        help="Rename directories/files under results/, logs/ on disk",
    )
    args = parser.parse_args()
    root = Path(args.root).resolve()
    n = update_repo(root)
    print(f"Updated {n} files under {root}")
    if args.move_results:
        move_results(root)


if __name__ == "__main__":
    main()
