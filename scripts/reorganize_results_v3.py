#!/usr/bin/env python3
"""v3 layout: kraken -> metagenomics/, pathogen outputs -> pathogen/, samples/ = merged host/mtDNA QC."""
from __future__ import annotations

import sys
from pathlib import Path

SKIP = {"reorganize_results_v3.py", "migrate_snakefile_paths_v2.py", "fix_wildcards_v2.py", "rename_layout_v2_labels.py"}

REPLACEMENTS = [
    # Kraken under metagenomics
    ("results/samples/{sample}/krakenuniq/", "results/metagenomics/krakenuniq/{sample}/"),
    ("results/samples/{wildcards.sample}/krakenuniq/", "results/metagenomics/krakenuniq/{wildcards.sample}/"),
    # Pathogen branch
    ("results/samples/{sample}/evalue/", "results/pathogen/{sample}/evalue/"),
    ("results/samples/{sample}/pathogen_mapping/", "results/pathogen/{sample}/pathogen_mapping/"),
    ("results/samples/{sample}/summary/", "results/pathogen/{sample}/summary/"),
    ("results/samples/{sample}/comparison/", "results/pathogen/{sample}/comparison/"),
    ("results/samples/{wildcards.sample}/evalue/", "results/pathogen/{wildcards.sample}/evalue/"),
    ("results/samples/{wildcards.sample}/pathogen_mapping/", "results/pathogen/{wildcards.sample}/pathogen_mapping/"),
    ("results/samples/{wildcards.sample}/summary/", "results/pathogen/{wildcards.sample}/summary/"),
    ("results/samples/{wc.sample}/pathogen_mapping/", "results/pathogen/{wc.sample}/pathogen_mapping/"),
    ("results/samples/{wc.sample}/evalue/", "results/pathogen/{wc.sample}/evalue/"),
    # Merged host/mtDNA QC -> samples (from pools)
    ("results/pools/qualimap_mtdna/{sample}/", "results/samples/{sample}/qualimap_mtdna/"),
    ("results/pools/qualimap/{sample}/", "results/samples/{sample}/qualimap/"),
    ("results/pools/qualimap_mtdna/{wildcards.sample}", "results/samples/{wildcards.sample}/qualimap_mtdna"),
    ("results/pools/qualimap/{wildcards.sample}", "results/samples/{wildcards.sample}/qualimap"),
    ("mkdir -p results/pools/qualimap_mtdna/", "mkdir -p results/samples/{wildcards.sample}/qualimap_mtdna"),
    ("mkdir -p results/pools/qualimap/", "mkdir -p results/samples/{wildcards.sample}/qualimap"),
    ("-outdir results/pools/qualimap_mtdna/", "-outdir results/samples/{wildcards.sample}/qualimap_mtdna/"),
    ("-outdir results/pools/qualimap/", "-outdir results/samples/{wildcards.sample}/qualimap/"),
    # Scripts / prose
    ("results/samples/*/pathogen_mapping/", "results/pathogen/*/pathogen_mapping/"),
    ("results/samples/*/bwa_pathogen/", "results/pathogen/*/pathogen_mapping/"),
    ("results/samples/{bio}/evalue/", "results/pathogen/{bio}/evalue/"),
    ("results/samples/{bio}/pathogen_mapping/", "results/pathogen/{bio}/pathogen_mapping/"),
    ("results/samples/{bio}/summary/", "results/pathogen/{bio}/summary/"),
    ("results/samples/{bio}/comparison/", "results/pathogen/{bio}/comparison/"),
    ("results/samples/{bio}/krakenuniq/", "results/metagenomics/krakenuniq/{bio}/"),
    ("results/samples/{sample_name}/comparison/", "results/pathogen/{sample_name}/comparison/"),
    ("results/samples/{sample_name}/", "results/pathogen/{sample_name}/"),
    # parse regex
    ("results/samples/([^/]+)/evalue/pathogen/", "results/pathogen/([^/]+)/evalue/pathogen/"),
    ('f"results/samples/{s}/evalue/pathogen/{s}_pathogen.csv"', 'f"results/pathogen/{s}/evalue/pathogen/{s}_pathogen.csv"'),
    ("results/samples/{s}/evalue/pathogen/", "results/pathogen/{s}/evalue/pathogen/"),
    # R report base
    ('file.path("results", "samples", sample)', 'file.path("results", "pathogen", sample)'),
    # Canonical path comment
    (
        "results/libraries/, results/samples/, results/pools/, results/final/, results/metagenomics/, results/workflow/.",
        "results/libraries/, results/samples/, results/pools/, results/pathogen/, results/metagenomics/, results/final/, results/workflow/.",
    ),
    (
        "samples/{bio}/       # krakenuniq, evalue, pathogen_mapping, summary, comparison",
        "samples/{bio}/       # merged host/mtDNA qualimap (and related sample QC)",
    ),
    (
        "samples/{bio}/            # krakenuniq, evalue, pathogen_mapping, summary, comparison",
        "samples/{bio}/            # merged host/mtDNA qualimap",
    ),
    (
        "pathogen/{bio}/           # evalue, pathogen_mapping, summary, comparison",
        "pathogen/{bio}/           # evalue, pathogen_mapping, summaries",
    ),
]

GLOBS = ["Snakefile", "scripts/*.py", "scripts/*.R", "docs/*.md", "README.md", "config/config.example.yaml"]


def migrate_text(text: str) -> str:
    for old, new in REPLACEMENTS:
        text = text.replace(old, new)
    return text


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    n = 0
    for pat in GLOBS:
        for path in root.glob(pat):
            if path.name in SKIP:
                continue
            t = path.read_text(encoding="utf-8")
            new = migrate_text(t)
            if new != t:
                path.write_text(new, encoding="utf-8")
                n += 1
    print(f"Updated {n} files in {root.resolve()}")


if __name__ == "__main__":
    main()
