#!/usr/bin/env python3
"""Second-pass fixes for Snakefile wildcards + remaining v1 paths."""
import re
from pathlib import Path

p = Path("Snakefile")
t = p.read_text(encoding="utf-8")
lib = [
    "qc", "adapter_removal", "fastq_screen", "host_mapping", "qualimap",
    "mtdna_mapping", "qualimap_mtdna", "flags", "prinseq", "bowtie2",
    "unaligned_fastq", "damageprofiler_host", "damageprofiler_mtdna",
]
for sub in lib:
    t = t.replace(
        f"results/{{wildcards.sample}}/{sub}",
        f"results/libraries/{{wildcards.sample}}/{sub}",
    )
t = t.replace(
    "results/{wildcards.sample}/bwa_pathogen",
    "results/samples/{wildcards.sample}/bwa_pathogen",
)
t = t.replace(
    "results/{wildcards.pcr}/bwa_pathogen",
    "results/libraries/{wildcards.pcr}/bwa_pathogen",
)
for old, new in [
    ("results/sample_unaligned_fastq", "results/pools/unaligned_fastq"),
    ("results/sample_host_mapping", "results/host/host_mapping"),
    ("results/sample_mtdna_mapping", "results/host/mtdna_mapping"),
    ("results/pools/host_mapping", "results/host/host_mapping"),
    ("results/pools/mtdna_mapping", "results/host/mtdna_mapping"),
    ("results/sample_bwa_pathogen", "results/pools/bwa_pathogen"),
]:
    t = t.replace(old, new)
t = re.sub(
    r"(\s+)escore(\s*=\s*\"results/samples)",
    r"\1evalue\2",
    t,
)
t = re.sub(r"input\.escore\b", "input.evalue", t)
t = t.replace(
    "Canonical Snakemake paths under results/{sample}/ are unchanged.",
    "Canonical paths: results/libraries/, results/samples/, results/pools/, "
    "results/cohort/, results/integrations/, results/workflow/.",
)
p.write_text(t, encoding="utf-8")
print("fixed", p)
