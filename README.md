# 🐖 PIGSTI                                                                                                                                                                                                                                                                (Pathogen anImal Genome Sequence ToolkIt)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org)
[![Conda](https://img.shields.io/badge/conda-enabled-green.svg)](https://conda.io)

A comprehensive bioinformatics pipeline for pathogen identification, genomic surveillance, and ancient DNA analysis using metagenomic sequencing data.

## 🧬 Overview

PIGSTI is a state-of-the-art bioinformatics pipeline designed for:
- **Pathogen Detection**: Automated identification of pathogens from metagenomic data
- **Ancient DNA Analysis**: Specialized tools for aDNA damage pattern analysis
- **Host Genome Analysis**: Comprehensive host and mitochondrial DNA mapping
- **Quality Control**: Extensive QC metrics and contamination screening
- **Comparative Analysis**: Integration of multiple detection methods (KrakenUniq, HOPS, BWA)

## ✨ Key Features

### 🔍 **Multi-Method Pathogen Detection**
- **KrakenUniq**: Taxonomic classification with unique k-mer counting
- **HOPS**: Hierarchical Oligonucleotide Primer Signature analysis
- **E-Score**: Statistical scoring for pathogen detection confidence
- **BWA Alignment**: Targeted alignment to pathogen references

### 🧬 **Ancient DNA Specialized Analysis**
- **DamageProfiler**: DNA damage pattern analysis
- **aDNA-BAMPlotter**: Visualization of ancient DNA characteristics
- **Soft-clipping Analysis**: Detection of DNA damage patterns
- **Entropy Profiling**: Sequence complexity analysis

### 📊 **Comprehensive Monitoring & Visualization**
- **Real-time Pipeline Monitoring**: Interactive execution heatmaps
- **Step Completion Tracking**: Visual progress monitoring
- **Execution Timing Analysis**: Performance metrics and bottleneck identification
- **Workflow Visualization**: Interactive pipeline diagram
- **Comprehensive Logging**: Detailed logs for every pipeline step

### 🎯 **Smart Sample Processing**
- **Adaptive Adapter Removal**: Automatic single-end vs paired-end detection
- **Dynamic Pathogen Selection**: Post-Escore targeted analysis
- **Quality-based Filtering**: Multi-level quality control
- **Contamination Screening**: Host and contaminant detection

## 🚀 Quick Start

### Prerequisites

- **Snakemake** ≥ 7.0.0
- **Conda** or **Mamba**
- **Python** 3.9+
- **Linux/macOS** environment

### Installation

1. **Clone the repository**:
```bash
git clone https://github.com/yourusername/PIGSTI.git
cd PIGSTI
```

2. **Set up the environment**:
```bash
# Create conda environment
conda env create -f environment.yml

# Or use mamba for faster installation
mamba env create -f environment.yml
```

3. **Configure the pipeline**:
```bash
# Edit configuration files
cp config/config.yaml.example config/config.yaml
cp config/samples.tsv.example config/samples.tsv

# Update paths in config/config.yaml
# Add your samples to config/samples.tsv
```

### Basic Usage

```bash
# Activate environment
conda activate pigsti

# Run the pipeline
snakemake --cores 40 --use-conda --rerun-incomplete

# Run with specific targets
snakemake --cores 40 --use-conda results/pipeline_execution_report.html
```

## 📁 Project Structure

```
PIGSTI/
├── config/                          # Configuration files
│   ├── config.yaml                  # Main pipeline configuration
│   ├── samples.tsv                  # Sample information
│   ├── Pathogen_spreadsheet.csv     # Pathogen reference database
│   └── config_hops2.0.txt          # HOPS configuration
├── scripts/                         # Pipeline scripts
│   ├── generate_pipeline_report.py  # Monitoring and reporting
│   ├── calculate_detection_scores.py # Pathogen scoring
│   └── ...                         # Other analysis scripts
├── workflow/
│   └── envs/                       # Conda environment definitions
├── results/                         # Pipeline outputs
│   ├── pipeline_execution_report.html # Interactive monitoring report
│   ├── pipeline_workflow_diagram.png # Workflow visualization
│   └── ...                         # Analysis results
├── logs/                           # Execution logs
└── Snakefile                       # Main pipeline definition
```

## 📋 Configuration

### Sample Configuration (`config/samples.tsv`)

```tsv
sample	r1	r2	RGLB	sequencing_run
Sample1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz	Sample1-LIB	Run1
Sample2	/path/to/Sample2_R1.fastq.gz		Sample2-LIB	Run2
```

**Note**: Leave `r2` empty for single-end samples.

### Main Configuration (`config/config.yaml`)

```yaml
# Database paths
kraken_db: "/path/to/krakenuniq/database"
host_index: "/path/to/host/bwa/index"

# BWA indices for different species
bwa_indices:
  "Homo sapiens": "/path/to/human/index"
  "Sus scrofa": "/path/to/pig/index"

mtDNA_indices:
  "Homo sapiens": "/path/to/human/mtdna/index"
  "Sus scrofa": "/path/to/pig/mtdna/index"

# Analysis parameters
min_read_length: 30
quality_threshold: 30
```

## 🔄 Pipeline Workflow

```mermaid
graph TD
    A[Raw FASTQ] --> B[Adapter Removal]
    B --> C[FastQ Screen]
    B --> D[Prinseq]
    D --> E[Bowtie2 Host]
    E --> F[KrakenUniq]
    F --> G[E-Score Analysis]
    C --> H[Host BWA]
    C --> I[mtDNA BWA]
    G --> J[Pathogen BWA]
    H --> K[DamageProfiler]
    I --> K
    J --> K
    H --> L[QualiMap]
    I --> L
    J --> L
    K --> M[Pathogen Detection]
    L --> M
    M --> N[Final Reports]
```

## 📊 Output Files

### Core Analysis Results
- **Pathogen Detection**: `results/{sample}/Escore/pathogen/{sample}_pathogen.csv`
- **Host Alignment**: `results/{sample}/bwa_host/{sample}.dedup.bam`
- **mtDNA Alignment**: `results/{sample}/bwa_mtdna/{sample}.dedup.bam`
- **Pathogen Alignments**: `results/{sample}/bwa_pathogen/{sample}_{pathogen}.dedup.bam`

### Quality Control Reports
- **FastQ Screen**: `results/{sample}/fastq_screen/{sample}.collapsed_screen.html`
- **QualiMap**: `results/{sample}/qualimap/genome_results.txt`
- **DamageProfiler**: `results/{sample}/damageprofiler_host/`

### Monitoring & Visualization
- **Pipeline Report**: `results/pipeline_execution_report.html`
- **Workflow Diagram**: `results/pipeline_workflow_diagram.png`
- **Timing Data**: `results/pipeline_timing_data.csv`
- **Detection Scores**: `results/pathogen_detection/detection_scores_matrix.csv`

## 🎛️ Advanced Usage

### Running Specific Steps

```bash
# Run only adapter removal
snakemake --cores 8 results/*/adapter_removal/*.collapsed.gz

# Run only pathogen detection
snakemake --cores 16 results/*/Escore/pathogen/*.csv

# Generate monitoring report
snakemake --cores 4 results/pipeline_execution_report.html
```

### Custom Configuration

```bash
# Use custom config file
snakemake --cores 40 --configfile custom_config.yaml

# Override specific parameters
snakemake --cores 40 --config min_read_length=50 quality_threshold=20
```

### Parallel Execution

```bash
# Run with cluster support
snakemake --cores 40 --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}"

# Dry run to check workflow
snakemake --cores 40 --dry-run
```

## 📈 Monitoring & Troubleshooting

### Real-time Monitoring

The pipeline includes comprehensive monitoring capabilities:

1. **Interactive Report**: Open `results/pipeline_execution_report.html` in your browser
2. **Execution Heatmap**: Visual progress tracking for all samples and steps
3. **Timing Analysis**: Performance metrics and bottleneck identification
4. **Log Files**: Detailed logs in `logs/` directory

### Common Issues

**Issue**: `KeyError: 'Sample' not found in samples.tsv`
**Solution**: Ensure all sample names in your analysis match exactly with `config/samples.tsv`

**Issue**: `ValueError: no index available for pileup`
**Solution**: The pipeline automatically creates BAM indices, but ensure sufficient disk space

**Issue**: Conda environment creation fails
**Solution**: Use `mamba` instead of `conda` for faster and more reliable package resolution

### Performance Optimization

- **Memory**: Use `--resources mem_mb=32000` for memory-intensive steps
- **Storage**: Ensure sufficient disk space (recommend 100GB+ per sample)
- **Cores**: Optimal performance with 20-40 cores depending on sample size

## 🧪 Testing

```bash
# Run test with example data
snakemake --cores 4 --use-conda test_data/

# Validate configuration
snakemake --cores 1 --dry-run --configfile config/config.yaml
```


## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.



---

**Made with ❤️ for the bioinformatics community**
