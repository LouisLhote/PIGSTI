#!/usr/bin/env python3
"""Apply v2 path layout + Evalue directory renames to text files. Run from repo root."""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Phase 1: structural paths (order: longest / most specific first)
STRUCTURAL = [
    ("results/sample_qualimap_mtdna/", "results/pools/qualimap_mtdna/"),
    ("results/sample_qualimap/", "results/pools/qualimap/"),
    ("results/sample_unaligned_fastq/", "results/pools/unaligned_fastq/"),
    ("results/sample_bwa_pathogen/", "results/pools/bwa_pathogen/"),
    ("results/sample_bwa_mtdna/", "results/pools/bwa_mtdna/"),
    ("results/sample_bwa_host/", "results/pools/bwa_host/"),
    ("results/KRAKENUNIQ_ABUNDANCE_MATRIX/", "results/metagenomics/kraken_abundance/"),
    ("results/pathogen_detection/", "results/integrations/pathogen_detection/"),
    ("results/hops/", "results/integrations/hops/"),
    ('directory("results/decOM")', 'directory("results/integrations/decOM")'),
    ("results/decOM", "results/integrations/decOM"),
    ("results/p_keys/", "results/integrations/decOM/p_keys/"),
    ("results/p_sink.txt", "results/integrations/decOM/p_sink.txt"),
    ("results/final/", "results/cohort/"),
    ("results/publication_figures/", "results/cohort/publication_figures/"),
    ("results/multi_qc_dashboard.html", "results/cohort/multi_qc_dashboard.html"),
    ("results/pipeline_", "results/cohort/pipeline_"),
    ("results/pathogen_reports_complete.txt", "results/workflow/pathogen_reports_complete.txt"),
    ("results/pathogen_targets.manifest.json", "results/workflow/pathogen_targets.manifest.json"),
    ("results/pathogen_targets.txt", "results/workflow/pathogen_targets.txt"),
    ("results/metagenomics_screening_cohort_ready.json", "results/workflow/metagenomics_screening_cohort_ready.json"),
    ("results/pathogen_bwa_complete.txt", "results/workflow/pathogen_bwa_complete.txt"),
    ("results/.pigsti_validation_ok", "results/workflow/.pigsti_validation_ok"),
    ("results/index_status/", "results/workflow/index_status/"),
    ("results/failures/", "results/workflow/failures/"),
    ("results/fastq_screen/fastq_screen.conf", "results/workflow/fastq_screen.conf"),
    ('BY_TYPE_ROOT = "results_by_type"', 'BY_TYPE_ROOT = "results/browse"'),
    ("results_by_type/", "results/browse/"),
    # Bio-sample subtrees under results/{sample}/
    ("results/{sample}/Escore/", "results/samples/{sample}/evalue/"),
    ("results/{sample}/krakenuniq/", "results/samples/{sample}/krakenuniq/"),
    ("results/{sample}/bwa_pathogen/", "results/samples/{sample}/bwa_pathogen/"),
    ("results/{sample}/comparison/", "results/samples/{sample}/comparison/"),
    ("results/{sample}/summary/", "results/samples/{sample}/summary/"),
    ("results/{wildcards.sample}/bwa_pathogen/", "results/samples/{wildcards.sample}/bwa_pathogen/"),
    ("results/{wildcards.sample}/krakenuniq/", "results/samples/{wildcards.sample}/krakenuniq/"),
    ("results/{wildcards.sample}/summary/", "results/samples/{wildcards.sample}/summary/"),
    ("results/{wildcards.pcr}/bwa_pathogen/", "results/libraries/{wildcards.pcr}/bwa_pathogen/"),
    ("results/{wc.pcr}/bwa_pathogen/", "results/libraries/{wc.pcr}/bwa_pathogen/"),
    ("results/{wc.sample}/bwa_pathogen/", "results/samples/{wc.sample}/bwa_pathogen/"),
    # PCR / library subtrees
    ("results/{sample}/adapter_removal/", "results/libraries/{sample}/adapter_removal/"),
    ("results/{sample}/prinseq/", "results/libraries/{sample}/prinseq/"),
    ("results/{sample}/fastq_screen/", "results/libraries/{sample}/fastq_screen/"),
    ("results/{sample}/bowtie2/", "results/libraries/{sample}/bowtie2/"),
    ("results/{sample}/unaligned_fastq/", "results/libraries/{sample}/unaligned_fastq/"),
    ("results/{sample}/bwa_host/", "results/libraries/{sample}/bwa_host/"),
    ("results/{sample}/bwa_mtdna/", "results/libraries/{sample}/bwa_mtdna/"),
    ("results/{sample}/qualimap_mtdna/", "results/libraries/{sample}/qualimap_mtdna/"),
    ("results/{sample}/qualimap/", "results/libraries/{sample}/qualimap/"),
    ("results/{sample}/damageprofiler_", "results/libraries/{sample}/damageprofiler_"),
    ("results/{sample}/damageprofiler/", "results/libraries/{sample}/damageprofiler/"),
    ("results/{sample}/qc/", "results/libraries/{sample}/qc/"),
    ("results/{sample}/flags/", "results/libraries/{sample}/flags/"),
    ("results/{wildcards.sample}/fastq_screen/", "results/libraries/{wildcards.sample}/fastq_screen/"),
    ("results/{wildcards.sample}/adapter_removal/", "results/libraries/{wildcards.sample}/adapter_removal/"),
    ("results/{wildcards.sample}/prinseq/", "results/libraries/{wildcards.sample}/prinseq/"),
    ("results/{wildcards.sample}/bowtie2/", "results/libraries/{wildcards.sample}/bowtie2/"),
    # Browse symlinks: legacy flat names -> v2 browse layout
    (f"{BY_TYPE_PLACEHOLDER}/sample_unaligned_fastq/", f"{BY_TYPE_PLACEHOLDER}/pools_unaligned_fastq/")
    if False else ("{BY_TYPE}/sample_unaligned_fastq/", "{BY_TYPE}/pools_unaligned_fastq/"),
]

# Fix BY_TYPE placeholder - use literal after BY_TYPE_ROOT change
STRUCTURAL = [t if not isinstance(t, tuple) or "BY_TYPE_PLACEHOLDER" not in str(t) else t for t in STRUCTURAL]
STRUCTURAL.extend([
    ("{BY_TYPE_ROOT}/Escore/", "{BY_TYPE_ROOT}/evalue/"),
    ("{BY_TYPE_ROOT}/sample_unaligned_fastq/", "{BY_TYPE_ROOT}/pools_unaligned_fastq/"),
    ("/Escore/", "/evalue/"),
])

# Phase 2: Snakemake rule / I/O identifiers (not config thresholds or script filenames)
IDENTIFIER_REPLACEMENTS = [
    ("_sample_from_escore_pathogen_csv", "_sample_from_evalue_pathogen_csv"),
    ("rule escore_by_type:", "rule evalue_by_type:"),
    ("rule escore:", "rule evalue:"),
    ("workflow/envs/escore.yaml", "workflow/envs/evalue.yaml"),
    ("logs/escore/", "logs/evalue/"),
    ("input.escore_pathogen", "input.evalue_pathogen"),
    ("input.escore_genus", "input.evalue_genus"),
    ("input.escore_species", "input.evalue_species"),
    ("input.escore_files", "input.evalue_files"),
    ("escore_pathogen=", "evalue_pathogen="),
    ("escore_genus=", "evalue_genus="),
    ("escore_species=", "evalue_species="),
    ("escore_files=", "evalue_files="),
    ("--escore ", "--evalue "),
    ("escore_paths", "evalue_paths"),
    ("``input.escore_files``", "``input.evalue_files``"),
    ("some escore jobs", "some evalue jobs"),
    ("completed E-score CSV", "completed E-value CSV"),
    ("E-score path", "E-value path"),
    ("from E-score path", "from E-value path"),
    ("E-score CSV", "E-value CSV"),
    ("# Load pathogen detection criteria (for Escore / reads", "# Load pathogen detection criteria (for E-value / reads"),
]

# Regex for legacy get_sample_ref_pairs default paths
LEGACY_EVALUE_CSV = re.compile(
    r'results/\{s\}/Escore/pathogen/\{s\}_pathogen\.csv'
)
LEGACY_EVALUE_CSV_V2 = 'results/samples/{s}/evalue/pathogen/{s}_pathogen.csv'

LEGACY_PARSE = re.compile(
    r'results/\(\[\^/\]\+\)/Escore/pathogen/'
)
LEGACY_PARSE_V2 = r'results/samples/([^/]+)/evalue/pathogen/'


def migrate_text(text: str) -> str:
    for old, new in STRUCTURAL:
        text = text.replace(old, new)
    for old, new in IDENTIFIER_REPLACEMENTS:
        text = text.replace(old, new)
    text = LEGACY_EVALUE_CSV.sub(
        lambda _: LEGACY_EVALUE_CSV_V2.replace("{s}", "{s}"), text
    )
    text = text.replace(
        "f\"results/{s}/Escore/pathogen/{s}_pathogen.csv\"",
        'f"results/samples/{s}/evalue/pathogen/{s}_pathogen.csv"',
    )
    text = text.replace(
        'if os.path.isfile(f"results/{s}/Escore/pathogen/{s}_pathogen.csv")',
        'if os.path.isfile(f"results/samples/{s}/evalue/pathogen/{s}_pathogen.csv")',
    )
    text = re.sub(
        r'results/\(\[\^/\]\+\)/Escore/pathogen/',
        LEGACY_PARSE_V2,
        text,
    )
    text = text.replace(
        'm = re.search(r"results/([^/]+)/Escore/pathogen/", path)',
        'm = re.search(r"results/samples/([^/]+)/evalue/pathogen/", path)',
    )
    text = text.replace(
        "Cannot parse biological sample from E-score path",
        "Cannot parse biological sample from E-value path",
    )
    return text


def migrate_file(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    new = migrate_text(text)
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    patterns = [
        "Snakefile",
        "scripts/*.py",
        "docs/*.md",
        "README.md",
        "config/config.example.yaml",
    ]
    changed = []
    for pat in patterns:
        for path in root.glob(pat):
            if path.name == "migrate_snakefile_paths_v2.py":
                continue
            if migrate_file(path):
                changed.append(path)
    print(f"Updated {len(changed)} files:")
    for p in changed:
        print(f"  {p}")


if __name__ == "__main__":
    main()
