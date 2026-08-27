# PIGSTI

**P**athogen an**I**mal **G**enome **S**equence **T**oolk**I**t

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org)
[![Conda](https://img.shields.io/badge/conda%20%2F%20mamba-enabled-green.svg)](https://conda.io)

**PIGSTI: a modular, reproducible pipeline for detecting species identity, pathogens, and microbes from animal palaeogenomic data.**

A [Snakemake](https://snakemake.github.io) workflow for shotgun libraries from ancient animal remains: competitive host identification, host/mtDNA mapping, metagenomic screening, pathogen reference mapping, and multi-criteria authentication.

<p align="center">
  <img src="docs/Fig1_v2.1.png" alt="PIGSTI pipeline overview (Figure 1)" width="100%">
</p>

<p align="center">
  <sub><b>Figure 1.</b> Schematic overview of the PIGSTI analysis workflow. Subworkflows: preprocessing, host identification and mapping, metagenomic screening, and pathogen detection and authentication. Line colours match each subworkflow (see legend); dashed lines are optional steps.</sub>
</p>

---

## Table of contents

- [What does PIGSTI do?](#what-does-pigsti-do)
- [What do you need?](#what-do-you-need)
- [Prepare the sample sheet](#prepare-the-sample-sheet)
- [Configure the pipeline](#configure-the-pipeline)
- [Quick start](#quick-start)
- [Outputs](#outputs)
- [Pathogen authentication](#pathogen-authentication)
- [Optional modules](#optional-modules)
- [Core tools](#core-tools)
- [Documentation](#documentation)
- [License](#license)

---

## What does PIGSTI do?

PIGSTI (**P**athogen an**I**mal **G**enome **S**equence **T**oolk**I**t) is a modular Snakemake pipeline for **species identity**, **pathogen screening**, and **microbial authentication** in shotgun data from ancient animal remains. After adapter trimming, reads follow **two parallel routes** — host mapping and metagenomic screening — before pathogen candidates are mapped and scored.

At a high level (see **Figure 1**):

1. **Preprocessing** — adapter trimming and quality filtering (**AdapterRemoval** for paired-end with collapse; **cutadapt** for single-end).
2. **Host route** — competitive host identification on collapsed reads (**FastQ Screen** / Bowtie2): species with the highest total mapped reads and the highest *one-hit-one-genome* count; optional full-dataset rescreen when endogenous content is very low (`fastq_screen_full_dataset_rescreen`). PIGSTI flags **human contamination**. Nuclear and mitochondrial mapping (**BWA** or **Bowtie2**), PCR replicate merging (**samtools**), 4 bp terminal soft-clipping (`softclip_mod.py`), **DamageProfiler**, **Qualimap**, endogenous DNA estimates, and optional genetic sexing (cattle, goat, sheep, dog).
3. **Metagenomics route** — in parallel, exact-duplicate removal and complexity filtering (**PRINSEQ++**), mammalian read depletion against a multi-host chimera index (**Bowtie2**), pooling of libraries, optional **decOM** source tracking, then **KrakenUniq** classification (optional **MALT/HOPS**).
4. **Pathogen route** — candidates from Guellil **E-value** screening (default) **∪ optional HOPS**, with per-pathogen thresholds in the spreadsheet; reference mapping (**BWA** or **Bowtie2**); **authentication** with a composite score out of 10 (or 13 with HOPS): KrakenUniq reads, E-value, relative entropy, edit-distance decay (damaged / non-damaged), 5′ C→T damage, ANI, breadth ratio, mapping ratio, and genus rank.
5. **Reports** — cohort Excel tables, per-pathogen PDFs, and a run provenance manifest.

**Two IDs matter throughout**

| ID | Meaning |
|----|---------|
| `pcr` | Sequencing library (one FASTQ pair / single-end file) |
| `sample` | Biological individual — PCR libraries are merged per sample after host mapping; pooled reads are merged at sample level for KrakenUniq, HOPS, and pathogen mapping |

**Host identification (FastQ Screen)** — By default, screening runs on a read subset; if `#One_hit_one_genome` is below 50, PIGSTI can rescreen the full collapsed FASTQ (`fastq_screen_full_dataset_min_one_hit: 50`). Host species is assigned from competitive mapping against your reference panel (see [`docs/CONFIG.md`](docs/CONFIG.md)).

---

## What do you need?

| Requirement | Details |
|-------------|---------|
| **OS** | Linux (recommended) or macOS |
| **Software** | [Miniconda](https://docs.conda.io/en/latest/miniconda.html) / [Mamba](https://mamba.readthedocs.io/), Snakemake ≥ 7.32 |
| **Config file** | `config/config.yaml` — easiest via the interactive HTML helper (below) |
| **Sample sheet** | Tab-separated `config/samples.tsv` (see [Prepare the sample sheet](#prepare-the-sample-sheet)) |
| **Pathogen table** | `config/Pathogen_spreadsheet.csv` — one reference FASTA per row (`bwa index` column); optional per-pathogen E-value / read thresholds |
| **Pathogen panel builder** | Optional: `python create_pigsti_pathogen_database.py` downloads references and builds indices |
| **KrakenUniq database** | NCBI NT–style DB used by **aMETA** — [doi:10.17044/scilifelab.20205504](https://doi.org/10.17044/scilifelab.20205504) · set `kraken_db:` to that directory |
| **Chimera / mammal index** | Bowtie2 prefix for host-read depletion (`host_index`) |
| **Host & mtDNA references** | FASTA paths (BWA) or Bowtie2 prefixes per species — depends on `host_aligner` |
| **Pathogen references** | One FASTA (or index) per spreadsheet row (`bwa index` column); build a panel with [`create_pigsti_pathogen_database.py`](create_pigsti_pathogen_database.py) |

Optional: HOPS MALT index + Resources; decOM sources.

> **KrakenUniq tip:** Point `kraken_db` at the unpacked aMETA NT database root (the folder that contains the KrakenUniq index files). You do not need to rebuild it for a standard PIGSTI run.

---

## Prepare the sample sheet

Copy the template and edit paths on your machine:

```bash
cp config/samples.example.tsv config/samples.tsv
```

**Required columns** (tab-separated):

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | yes | Biological sample ID |
| `pcr` | recommended | Library ID (defaults to `sample` if omitted) |
| `r1` | yes | Absolute or repo-relative path to R1 FASTQ (`.fastq.gz`) |
| `r2` | for PE | R2 path; leave empty for single-end |
| `RGLB` | recommended | Read-group library ID |
| `sequencing_run` | optional | Run / batch label |

**Example**

```tsv
sample	pcr	r1	r2	RGLB	sequencing_run
Pig01	Pig01_PCR1	/path/to/Pig01_PCR1_R1.fastq.gz	/path/to/Pig01_PCR1_R2.fastq.gz	LIB01	Run1
Pig01	Pig01_PCR2	/path/to/Pig01_PCR2_R1.fastq.gz	/path/to/Pig01_PCR2_R2.fastq.gz	LIB02	Run1
Sample02	Sample02_SE	/path/to/Sample02_R1.fastq.gz		LIB03	Run1
```

- Same `sample`, different `pcr` → libraries are merged at the biological-sample level after host mapping and filtering, and for KrakenUniq / pathogen mapping.
- Paths must exist before the run (startup validation checks them).

Also prepare the pathogen spreadsheet:

```bash
cp config/Pathogen_spreadsheet.example.csv config/Pathogen_spreadsheet.csv
```

Fill `Krakenuniq name`, `Hops name`, and `bwa index` (pathogen FASTA path). Optional per-pathogen overrides: `Guellil_et_al_Evalue_threshold`, `min_reads`, `min_escore` (e.g. stricter filtering for commensals such as *E. coli*).

---

## Configure the pipeline

### Recommended — HTML config facilitator

Open this file in any browser, fill the form, then **download / copy** `config.yaml`:

**[`config/pigsti_config_generator.html`](config/pigsti_config_generator.html)**

It covers manifests, adapters, aligners, host references, Kraken / chimera paths, optional HOPS (including parallel MALT), detection thresholds (E-value vs E-score), and authentication criteria.

### Or edit YAML by hand

```bash
cp config/config.example.yaml config/config.yaml
```

Minimum keys to set:

```yaml
samples: "config/samples.tsv"
pathogen_spreadsheet: "config/Pathogen_spreadsheet.csv"

kraken_db: "/path/to/aMETA_NT_krakenuniq_db"   # from doi:10.17044/scilifelab.20205504
host_index: "/path/to/refs/chimera"            # Bowtie2 prefix

host_aligner: bwa
pathogen_aligner: bwa
pathogen_mapping_mode: default                 # recommended; or super_careful

bwa_indices:
  Pig: "/path/to/sus_scrofa.fa"
mtDNA_indices:
  Pig: "/path/to/pig_mito.fa"

pathogen_detection_criteria:
  use_evalue_for_detection: true
  guellil_evalue_threshold: 0.001
  reads_threshold: 50
```

When `host_aligner: bowtie2`, use `bowtie2_indices` and `bowtie2_mtDNA_indices` instead of `bwa_indices` / `mtDNA_indices`.

Full key list: [`docs/CONFIG.md`](docs/CONFIG.md).

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/LouisLhote/PIGSTI.git
cd PIGSTI

# 2. Driver environment (per-rule conda envs are created on first run)
conda env create -n pigsti-snake -f PIGSTI_snakemake.yaml
conda activate pigsti-snake

# 3. Config + manifests
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/Pathogen_spreadsheet.example.csv config/Pathogen_spreadsheet.csv
# → edit paths, or use config/pigsti_config_generator.html

# 4. Dry-run (validation runs at Snakefile load)
snakemake -n -p --use-conda --conda-frontend mamba \
  --conda-prefix .snakemake/conda --cores 32

# 5. Run
snakemake --use-conda --conda-frontend mamba \
  --conda-prefix .snakemake/conda --cores 32 --rerun-incomplete
```

Write results to another disk:

```bash
snakemake --use-conda --cores 32 --config results_root=/path/to/output
```

---

## Outputs

```
results/
├── libraries/{pcr}/       adapter_removal, fastq_screen, host_mapping/, mtdna_mapping/
├── samples/{bio}/         merged Qualimap, optional sexing
├── pools/                 unaligned FASTQs, merged BAMs
├── metagenomics/          krakenuniq/, hops/, decOM/
├── pathogen/{bio}/        evalue/, pathogen_mapping/, summary/
├── final/                 cohort Excel, heatmaps, run_manifest.json
└── workflow/              checkpoints, validation stamp
```

| Deliverable | Path |
|-------------|------|
| Pathogen cohort table | `results/final/pathogen_summary_all_samples.xlsx` |
| Host + mtDNA cohort table | `results/final/host_mtdna_summary_all_samples.xlsx` |
| Comprehensive table | `results/final/comprehensive_summary_all_samples.xlsx` |
| Per-pathogen PDF | `results/pathogen/{bio}/summary/{bio}_{pathogen}_pathogen_report.pdf` |
| E-value / Kraken table | `results/pathogen/{bio}/evalue/pathogen/{bio}_pathogen.csv` |
| Provenance | `results/final/run_manifest.json` |

---

## Pathogen authentication

Screening uses KrakenUniq clade read count (default **≥ 50**) and Guellil **E-value** (default **E > 0.001**), with optional **HOPS** union. Mapped candidates are scored **out of 10** (or **13** with HOPS). Full definitions: [`docs/METRICS_DEFINITIONS.md`](docs/METRICS_DEFINITIONS.md).

| Criterion | Default | Role |
|-----------|---------|------|
| KrakenUniq clade reads | ≥ 50 | Screening abundance |
| Guellil E-value *E = (K/R) × C* | > 0.001 | Screening confidence |
| Relative entropy | ≥ 0.9 (virus ≥ 0.7) | Evenness of read starts (100 bp / 1 kb windows) |
| Edit-distance decay (damaged / non-damaged) | ≥ 0.65 / ≥ 0.55 | Terminal-deamination split; composite decay score |
| 5′ C→T damage | ≥ 0.01 | Postmortem deamination (DamageProfiler) |
| ANI | > 96.5% | Similarity to the reference |
| Breadth ratio | ≥ 0.8 | Observed / expected Poisson breadth |
| Mapping ratio | ≥ 0.5 | Mapped vs KrakenUniq-assigned reads |
| Genus ranking | = 1 | Dominant genus-level KrakenUniq hit |
| HOPS / MaltExtract (optional) | +3 criteria | Edit-distance decline, terminal damage, damaged ED = 0 |

---

## Optional modules

| Flag | Effect |
|------|--------|
| `enable_hops: true` | HOPS (MALT + MaltExtract); candidates = E-value hits ∪ HOPS hits. Supports `hops_parallel` / `hops_malt_mmap` |
| `enable_decom: true` | decOM source tracking |
| `enable_sexing: true` | Residual sexing (Cow, Goat, Sheep, Dog) |
| `pathogen_screening_only: true` | Skip host/mtDNA mapping |
| `cleanup_intermediates: true` | Drop large intermediates after finals |

---

## Core tools

| Stage | Function | Tool |
|-------|----------|------|
| Preprocessing | Adapter / quality trim | **AdapterRemoval** (PE) · **cutadapt** (SE) |
| Host ID | Competitive mapping | **FastQ Screen** (Bowtie2) |
| Host / mtDNA | Alignment · dedup | **BWA aln** or **Bowtie2** · **samtools** / **Picard** |
| Host QC | Damage · mapping QC · soft-clip | **DamageProfiler** · **Qualimap2** · `softclip_mod.py` (4 bp) |
| Host QC | Genetic sex (optional) | Residual method (cattle, goat, sheep, dog) |
| Metagenomics | Dedup · host removal | **PRINSEQ++** · **Bowtie2** (chimera index) |
| Metagenomics | Classification · optional | **KrakenUniq** · **HOPS/MALT** · **decOM** |
| Pathogen | Detection · mapping · auth | Guellil **E-value** · **BWA/Bowtie2** · composite score |

---

## Documentation

| Document | Contents |
|----------|----------|
| [`docs/CONFIG.md`](docs/CONFIG.md) | Config keys, HOPS parallel, `results_root` |
| [`docs/METRICS_DEFINITIONS.md`](docs/METRICS_DEFINITIONS.md) | Authentication metrics (score out of 10 / 13) |
| [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md) | Result directory layout |

---

## License

[MIT](LICENSE) — © 2026 Louis Lhote  

**Issues:** [github.com/LouisLhote/PIGSTI/issues](https://github.com/LouisLhote/PIGSTI/issues)
