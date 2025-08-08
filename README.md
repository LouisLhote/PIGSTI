# 🐖 PIGSTI

**PIGSTI** (Pathogen anImal Genome Sequence ToolkIt) is a Snakemake-based pipeline designed to **detect pathogens** and **screen animal endogenous sequences** in ancient DNA (aDNA) sequencing datasets.

---

## 🚀 Features

- **Pathogen detection** – Automated pipeline to identify microbial or pathogen DNA in aDNA samples.  
- **Animal endogenous screening** – Checks for and quantifies:
  - Species-specific endogenous DNA  
  - Human contamination  
  - Mitochondrial genome coverage  
- **Modular workflow** – Structured via Snakemake for reproducibility, scalability, and easy maintenance.

---

## 🔧 Installation

# 1️⃣ Install Miniconda (Linux)

If not already installed:

### Download Miniconda installer (Linux, Python 3)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
### Install Miniconda
```
bash Miniconda3-latest-Linux-x86_64.sh
```

### Follow prompts, then reload shell
```
source ~/.bashrc
```

# 2️⃣ Download PIGSTI and Install Snakemake

### Clone the repository
```
git clone https://github.com/LouisLhote/PIGSTI.git
cd PIGSTI
```
### Create the conda environment from the provided YAML
```
conda env create -n pigsti-snake -f PIGSTI_snakemake.yaml
```

### Activate the environment
```
conda activate pigsti-snake
```
---

# ▶️ Running the Pipeline

Once the environment is active, run the workflow with:
```
snakemake --cores <N> --use-conda
```
Replace `<N>` with the number of CPU cores you want to use.

---


