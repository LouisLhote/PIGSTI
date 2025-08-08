# PIGSTI

**PIGSTI** (Pathogen anImal Genome Sequence ToolkIt) is a Snakemake-based pipeline designed to **detect pathogens** and **screen animal endogenous sequences** in ancient DNA (aDNA) sequencing datasets.

---

## 🚀 Features

- **Pathogen detection**: Automated pipeline to identify microbial or pathogen DNA in aDNA samples.
- **Animal endogenous screening**: Checks for and quantifies endogenous animal DNA contamination or presence.
- **Modular workflow**: Structured via Snakemake for reproducibility, scalability, and manageability.
- **Configurable steps**: Easily adjustable through configuration files (`config/` directory).

---

## 🔧 Installation

### 1️⃣ Install Miniconda (Linux)

If not already installed:

```bash
# Download Miniconda installer (Linux, Python 3)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts, then reload shell or source .bashrc
source ~/.bashrc
