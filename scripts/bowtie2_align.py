import os

if "snakemake" in globals():
    reads_file = snakemake.input.reads
    species_file = snakemake.input.species
    output_file = snakemake.output.bam
    threads = snakemake.threads
    config = snakemake.config
else:
    raise RuntimeError("This script must be run via Snakemake")

# Read species exactly as-is
with open(species_file) as f:
    species = f.read().strip()

# Lookup bowtie2 index prefix in config
bt2_indices = config.get("bowtie2_indices", {})
index_prefix = bt2_indices.get(species)
if index_prefix is None:
    fallback = "Sheep"
    print(f"Warning: Species '{species}' not found in bowtie2_indices config, falling back to '{fallback}'.")
    index_prefix = bt2_indices.get(fallback)
    if index_prefix is None:
        raise ValueError(f"No bowtie2 index for species '{species}' or fallback '{fallback}' found in config.")

# Construct bowtie2 + samtools command
cmd = (
    f"bowtie2 --very-sensitive-local -x {index_prefix} -U {reads_file} -p {threads} "
    f"| samtools view -Sb - -F4 > {output_file}"
)

print(f"Running command:\n{cmd}")
exit_code = os.system(cmd)

if exit_code != 0:
    raise RuntimeError(f"bowtie2 alignment failed with exit code {exit_code}")
