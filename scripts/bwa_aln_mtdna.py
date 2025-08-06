import os

if "snakemake" in globals():
    reads_file = snakemake.input.reads
    species_file = snakemake.input.species
    output_file = snakemake.output.bam
    threads = snakemake.threads
    config = snakemake.config
    sample = snakemake.wildcards.sample
else:
    raise RuntimeError("This script must be run via Snakemake")

# Read species
with open(species_file) as f:
    species = f.read().strip()

# Get reference index path
bwa_indices = config.get("mtDNA_indices", {})
index_prefix = bwa_indices.get(species)
if index_prefix is None:
    fallback = "Sheep"
    print(f"Warning: Species '{species}' not found, using fallback '{fallback}'.")
    index_prefix = bwa_indices.get(fallback)
    if index_prefix is None:
        raise ValueError(f"No BWA index for species '{species}' or fallback '{fallback}' found.")

# Output .sai and reference path
sai_file = output_file.replace("_F4.cram", ".sai")
ref_fasta = index_prefix  # assuming .fa path is same as index prefix

# Ensure output dir exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

sai_file = output_file.replace(".bam", ".sai")
bam_file = output_file  # final BAM output

cmd = (
    f"bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {index_prefix} {reads_file} > {sai_file} && "
    f"bwa samse -r '@RG\\tID:{sample}_host\\tSM:{sample}\\tPL:ILLUMINA' {index_prefix} {sai_file} {reads_file} | "
    f"samtools view -@ {threads} -F 4 -b -o {bam_file} -"
)

print(f"Running command:\n{cmd}")
exit_code = os.system(cmd)
if exit_code != 0:
    raise RuntimeError(f"BWA mapping failed with exit code {exit_code}")

