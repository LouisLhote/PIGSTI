#!/usr/bin/env python3
"""Rename v2 layout: cohort->final, integrations->metagenomics, bwa_pathogen->pathogen_mapping."""
from __future__ import annotations

import sys
from pathlib import Path

SKIP = {
    "rename_layout_v2_labels.py",
    "migrate_snakefile_paths_v2.py",
    "fix_wildcards_v2.py",
}


def migrate_text(text: str) -> str:
    text = text.replace("results/cohort/", "results/final/")
    text = text.replace("results/integrations/", "results/metagenomics/")
    text = text.replace("/bwa_pathogen/", "/pathogen_mapping/")
    text = text.replace("/bwa_pathogen", "/pathogen_mapping")
    text = text.replace("bwa_pathogen/", "pathogen_mapping/")
    text = text.replace('"bwa_pathogen"', '"pathogen_mapping"')
    text = text.replace("'bwa_pathogen'", "'pathogen_mapping'")
    text = text.replace("pools/bwa_pathogen", "pools/pathogen_mapping")
    text = text.replace("sample_bwa_pathogen", "sample_pathogen_mapping")
    # prose in docs (avoid touching workflow json basename with 'cohort' in the middle)
    text = text.replace("`cohort/`", "`final/`")
    text = text.replace("`integrations/`", "`metagenomics/`")
    text = text.replace("results/cohort", "results/final")
    text = text.replace("results/integrations", "results/metagenomics")
    return text


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    globs = ["Snakefile", "scripts/*.py", "scripts/*.R", "docs/*.md", "README.md", "config/config.example.yaml"]
    changed = []
    for pat in globs:
        for path in root.glob(pat):
            if path.name in SKIP:
                continue
            text = path.read_text(encoding="utf-8")
            new = migrate_text(text)
            if new != text:
                path.write_text(new, encoding="utf-8")
                changed.append(path)
    print(f"Updated {len(changed)} files under {root.resolve()}")


if __name__ == "__main__":
    main()
