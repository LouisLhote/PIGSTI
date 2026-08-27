import os
import subprocess
import sys
import time

if "snakemake" in globals():
    reads_file = snakemake.input.reads
    species_file = snakemake.input.species
    host_mapped_bam = snakemake.output.bam
    host_unmapped_bam = snakemake.output.unmapped_bam
    host_unmapped_fastq = snakemake.output.unmapped_fastq
    threads = snakemake.threads
    config = snakemake.config
    sample = snakemake.wildcards.sample
    lb = snakemake.params.lb
    rgid = snakemake.params.rgid
else:
    raise RuntimeError("This script must be run via Snakemake")


def _ensure_aligner_utils_path():
    here = os.path.dirname(os.path.abspath(__file__))
    cur = here
    for _ in range(12):
        sd = os.path.join(cur, "scripts")
        if os.path.isfile(os.path.join(sd, "aligner_ref_utils.py")):
            if sd not in sys.path:
                sys.path.insert(0, sd)
            return
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    raise RuntimeError("Could not locate scripts/aligner_ref_utils.py (run via Snakemake from repo root).")


_ensure_aligner_utils_path()
from aligner_ref_utils import resolve_bowtie2_prefix_and_fasta, resolve_bwa_database_prefix

# Read species
with open(species_file) as f:
    species = f.read().strip()

# Get reference index path(s)
bwa_indices = config.get("bwa_indices", {})
bowtie2_indices = config.get("bowtie2_indices", {})

host_aligner = config.get("host_aligner", "bwa")

# Make read-group values robust: avoid passing empty rg-id/LB to bowtie2.
lb = (str(lb).strip() if lb is not None else "")
rgid = (str(rgid).strip() if rgid is not None else "")
if not rgid:
    rgid = sample + "_host"

fallback = "Sheep"

if host_aligner == "bowtie2":
    bt2_cfg = bowtie2_indices.get(species) or bowtie2_indices.get(fallback)
    fa_cfg = bwa_indices.get(species) or bwa_indices.get(fallback)
    if not bt2_cfg and not fa_cfg:
        raise ValueError(
            f"No bowtie2_indices or bwa_indices entry for species '{species}' "
            f"(or fallback '{fallback}') with host_aligner=bowtie2."
        )
    index_prefix, ref_fasta = resolve_bowtie2_prefix_and_fasta(bt2_cfg, fa_cfg)
else:
    fa_cfg = bwa_indices.get(species) or bwa_indices.get(fallback)
    if not fa_cfg:
        raise ValueError(
            f"No bwa_indices entry for species '{species}' (or fallback '{fallback}') with host_aligner=bwa."
        )
    index_prefix = resolve_bwa_database_prefix(fa_cfg)
    ref_fasta = None

def ensure_bowtie2_index(index_prefix: str, ref_fasta: str, timeout_sec: int = 3600):
    """
    Ensure that a bowtie2 index exists and is complete at the given prefix.
    Uses a simple lock file and a 'done' marker to avoid race conditions
    when many samples try to build the same index concurrently.
    """
    needed_files = [
        index_prefix + ".1.bt2",
        index_prefix + ".2.bt2",
        index_prefix + ".3.bt2",
        index_prefix + ".4.bt2",
        index_prefix + ".rev.1.bt2",
        index_prefix + ".rev.2.bt2",
    ]
    done_marker = index_prefix + ".bt2_build_done"
    lock_file = index_prefix + ".bt2_build_lock"

    def index_complete() -> bool:
        # Consider the index complete if all .bt2 files are present, even if
        # our local done-marker is missing (e.g. prebuilt system index).
        return all(os.path.exists(p) for p in needed_files)

    # Fast path: already complete
    if index_complete():
        return

    # Try to acquire build lock
    start = time.time()
    while True:
        try:
            fd = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            # We hold the lock: (re)build the index
            os.makedirs(os.path.dirname(index_prefix), exist_ok=True)
            build_cmd = f"bowtie2-build {ref_fasta} {index_prefix}"
            print(f"[bwa_aln_host] Bowtie2 index not complete at {index_prefix}, building with:\n{build_cmd}")
            rc = os.system(build_cmd)
            if rc != 0:
                # Best effort cleanup
                if os.path.exists(lock_file):
                    os.remove(lock_file)
                raise RuntimeError(f"bowtie2-build failed with exit code {rc}")
            # Mark as done (ignore errors if directory is not writable, since the index itself is what matters)
            try:
                open(done_marker, "w").close()
            except OSError:
                pass
            if os.path.exists(lock_file):
                os.remove(lock_file)
            return
        except FileExistsError:
            # Someone else is building; wait for completion
            if index_complete():
                return
            if time.time() - start > timeout_sec:
                raise RuntimeError(
                    f"Timeout while waiting for bowtie2 index at {index_prefix} "
                    f"to be built (lock file '{lock_file}' still present)."
                )
            time.sleep(10)


# If using bowtie2, auto-build index if missing / incomplete (with concurrency safety)
if host_aligner == "bowtie2":
    ensure_bowtie2_index(index_prefix, ref_fasta)

# Determine pathogen mapping mode (default vs super_careful)
pathogen_mode = config.get("pathogen_mapping_mode", "super_careful")

# Ensure output dir exists
os.makedirs(os.path.dirname(host_mapped_bam), exist_ok=True)


def run_bowtie2_pipe_to_bam(
    bowtie2_cmd: list,
    log_file: str,
    sam_threads: int,
    out_bam: str,
) -> None:
    """
    Run bowtie2 | samtools view without a shell pipeline so both exit codes are checked.
    """
    os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
    with open(log_file, "wb") as err_f:
        p_bt2 = subprocess.Popen(
            bowtie2_cmd,
            stdout=subprocess.PIPE,
            stderr=err_f,
        )
        try:
            p_st = subprocess.Popen(
                [
                    "samtools",
                    "view",
                    f"-@{sam_threads}",
                    "-F",
                    "4",
                    "-b",
                    "-o",
                    out_bam,
                    "-",
                ],
                stdin=p_bt2.stdout,
            )
        finally:
            if p_bt2.stdout is not None:
                p_bt2.stdout.close()
        st_rc = p_st.wait()
        bt2_rc = p_bt2.wait()
    if bt2_rc != 0:
        raise RuntimeError(
            f"bowtie2 failed with exit code {bt2_rc} (see {log_file}); "
            f"samtools exit code was {st_rc}"
        )
    if st_rc != 0:
        raise RuntimeError(f"samtools view failed with exit code {st_rc}")


if host_aligner == "bowtie2":
    log_file = host_mapped_bam.replace(".bam", "_bowtie2.log")
    bowtie2_cmd = [
        "bowtie2",
        "--end-to-end",
        "--sensitive",
        "-x",
        index_prefix,
        "-U",
        reads_file,
        "-p",
        str(threads),
        "--rg-id",
        rgid,
        "--rg",
        f"SM:{sample}",
        "--rg",
        "PL:ILLUMINA",
    ]
    if lb:
        bowtie2_cmd.extend(["--rg", f"LB:{lb}"])
    if pathogen_mode != "default":
        bowtie2_cmd.extend(["--un-gz", host_unmapped_fastq])
    print(f"Running bowtie2 -> samtools (log: {log_file})")
    run_bowtie2_pipe_to_bam(bowtie2_cmd, log_file, threads, host_mapped_bam)
    if pathogen_mode == "default":
        subprocess.run(
            [
                "bash",
                "-c",
                f': > "{host_unmapped_bam}" && printf "" | pigz > "{host_unmapped_fastq}"',
            ],
            check=True,
        )
    else:
        subprocess.run(["bash", "-c", f': > "{host_unmapped_bam}"'], check=True)
    cmd = None
else:
    # BWA host mapping (original behaviour)
    # We keep a single SAI file and reuse it
    sai_file = host_mapped_bam.replace(".bam", ".sai")
    if pathogen_mode == "default":
        # Fast mode: only produce host-mapped BAM; create empty placeholders for unaligned outputs
        cmd = (
            f"bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {index_prefix} {reads_file} > {sai_file} && "
            f"bwa samse -r '@RG\\tID:{sample}_host\\tSM:{sample}\\tPL:ILLUMINA' "
            f"{index_prefix} {sai_file} {reads_file} | "
            f"samtools view -@ {threads} -F 4 -b -o {host_mapped_bam} && "
            f': > "{host_unmapped_bam}" && '
            f'printf "" | pigz > "{host_unmapped_fastq}"'
        )
    else:
        # Super-careful mode: also produce host-unaligned BAM + FASTQ for pathogen mapping
        cmd = (
            f"bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {index_prefix} {reads_file} > {sai_file} && "
            f"bwa samse -r '@RG\\tID:{sample}_host\\tSM:{sample}\\tPL:ILLUMINA' "
            f"{index_prefix} {sai_file} {reads_file} | "
            f"samtools view -@ {threads} -F 4 -b -o {host_mapped_bam} && "
            f"bwa samse -r '@RG\\tID:{sample}_host\\tSM:{sample}\\tPL:ILLUMINA' "
            f"{index_prefix} {sai_file} {reads_file} | "
            f"samtools view -@ {threads} -f 4 -b -o {host_unmapped_bam} && "
            f"samtools bam2fq {host_unmapped_bam} | pigz > {host_unmapped_fastq}"
        )

if cmd is not None:
    print(f"Running command:\n{cmd}")
    try:
        subprocess.run(
            ["bash", "-c", f"set -euo pipefail; {cmd}"],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"BWA host alignment (mapped{' + unaligned' if pathogen_mode != 'default' else ''}) "
            f"failed with exit code {e.returncode}"
        ) from e

