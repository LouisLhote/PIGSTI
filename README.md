
# 🐖 PIGSTI – Pathogen Identification and Genomic Screening of Targeted Individuals

**PIGSTI** is a Snakemake pipeline for screening ancient or modern sequencing libraries for pathogen DNA, combining authenticity checks, metagenomic classification, and reference-based analysis.

---

## 🚀 Features

- Supports paired-end and single-end reads
- Automatically links samples to candidate pathogens via Escore + spreadsheet match
- Generates:
  - Authenticity profiles (DamageProfiler, edit distance, read length)
  - aDNA-specific BAM plots and entropy maps
  - HOPS reports (if applicable)
  - Summary PDFs per sample-pathogen pair
- Compatible with high-throughput batch processing via Snakemake

---

## 📁 Directory Structure

```bash
.
├── config/
│   ├── config.yaml                # Global settings
│   ├── samples.tsv                # Sample metadata
│   └── Pathogen_spreadsheet.csv  # List of candidate pathogens and metadata
├── results/                       # Output directory (auto-generated)
├── Snakefile                      # Workflow definition
└── README.md                      # You're here
```

---

## 🔧 Installation

### 🐍 Install Miniconda (Linux)

```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts and restart your shell or source your .bashrc
source ~/.bashrc
```

### 📦 Install Snakemake and dependencies

```bash
conda create -n pigsti -c bioconda snakemake
conda activate pigsti
```

Optional: set up tools like KrakenUniq, DamageProfiler, HOPS, samtools, bwa, etc.

---

## 📋 Input Files

### `config/samples.tsv`

| sample | pcr_id | r1 | r2 | RGLB | sequencing_run |
|--------|--------|----|----|------|----------------|
| S1     | P1     | S1_R1.fastq.gz | S1_R2.fastq.gz | Lib1 | SeqRun1 |

### `config/Pathogen_spreadsheet.csv`

- Must contain a column named **`Krakenuniq name`**, which matches taxonomic names from Escore output.
- Optional columns for HOPS support or plotting metadata.

### `config/config.yaml`

Basic YAML file for any shared parameters (e.g. reference directory, thread counts, etc.).

---

## 🧬 How It Works

1. Reads sample sheet and builds a dictionary of paired reads and read groups.
2. Reads Escore output to find high-confidence pathogen matches.
3. Cross-references with the master spreadsheet to generate a list of `(sample, pathogen)` targets.
4. Runs rules for:
   - Adapter removal
   - Alignment
   - Authenticity analysis (DamageProfiler)
   - KrakenUniq classification
   - BAM plotting and coverage entropy maps
   - Final reporting

---

## 🏁 Running the Pipeline

To run locally:

```bash
snakemake --use-conda --cores 8
```

To run on a cluster or HPC:

```bash
snakemake --profile your-cluster-profile
```

---

## 📤 Output Summary

For each `(sample, pathogen)` pair:
- Escore CSV
- HOPS PDF report (if detected)
- DamageProfiler plots (damage, read length, edit distance)
- aDNA-BAMPlotter figures
- Coverage/ANI plots
- Entropy maps (1kb sliding windows)

A final combined **multi-panel PDF report** is generated per sample-pathogen pair.

---

## 🔐 Notes & Tips

- Escore outputs must be located at:  
  `results/{sample}/Escore/pathogen/{sample}_pathogen.csv`
- HOPS summary PDFs must be found at:  
  `results/hops/maltExtract/pdf_candidate_profiles/{Hops name}/{prefix}_summary.pdf`
- If Escore results are missing for a sample, it is skipped.
- Only pathogens present in both Escore *and* the spreadsheet are processed.

---

## 📜 Citation

If you use this pipeline, please cite the underlying tools:

- [Snakemake](https://snakemake.readthedocs.io/)
- [DamageProfiler](https://github.com/Integrative-Transcriptomics/DamageProfiler)
- [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
- [HOPS](https://github.com/rhuebler/HOPS)

---

## 🧑‍💻 Author

**Louis Lhote**  
GitHub: [@LouisLhote](https://github.com/LouisLhote)
