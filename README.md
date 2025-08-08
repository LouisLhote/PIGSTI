## 🔧 Installation

### 1️⃣ Install Miniconda (Linux)

If not already installed:

```bash
# Download Miniconda installer (Linux, Python 3)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install Miniconda
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts, then reload shell
source ~/.bashrc
2️⃣ Download PIGSTI and Install Snakemake
bash
Copy
Edit
# Clone the repository
git clone https://github.com/LouisLhote/PIGSTI.git
cd PIGSTI

# Create the conda environment from the provided YAML
conda env create -n pigsti-snake -f PIGSTI_snakemake.yaml

# Activate the environment
conda activate pigsti-snake
▶️ Running the Pipeline
Once the environment is active, run the workflow with:

bash
Copy
Edit
snakemake --cores <N> --use-conda
