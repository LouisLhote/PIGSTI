import os
import re
import subprocess
import sys
import time

if "snakemake" in globals():
    reads_file = snakemake.input.reads
    species_file = snakemake.input.species
    output_file = snakemake.output.bam
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
bwa_indices = config.get("mtDNA_indices", {})
bowtie2_indices = config.get("bowtie2_mtDNA_indices", {})

host_aligner = config.get("host_aligner", "bwa")

# Make read-group values robust: avoid passing empty rg-id/LB to bowtie2.
lb = (str(lb).strip() if lb is not None else "")
rgid = (str(rgid).strip() if rgid is not None else "")
if not rgid:
    rgid = sample + "_mtdna"
# Bowtie2 / SAM RG ID: avoid whitespace and odd characters that can confuse parsers.
rgid = re.sub(r"[^\w.\-]+", "_", rgid).strip("_") or (sample + "_mtdna")


def ensure_bowtie2_index(index_prefix: str, ref_fasta: str, timeout_sec: int = 3600):
    """
    Ensure a complete Bowtie2 index exists at index_prefix (all six .bt2 files).
    Uses a lock file to avoid concurrent bowtie2-build on the same prefix.
    """
    needed_files = [
        index_prefix + ".1.bt2",
        index_prefix + ".2.bt2",
        index_prefix + ".3.bt2",
        index_prefix + ".4.bt2",
        index_prefix + ".rev.1.bt2",
        index_prefix + ".rev.2.bt2",
    ]
    lock_file = index_prefix + ".bt2_build_lock"

    def index_complete() -> bool:
        return all(os.path.exists(p) for p in needed_files)

    if index_complete():
        return

    start = time.time()
    while True:
        try:
            fd = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            os.makedirs(os.path.dirname(index_prefix) or ".", exist_ok=True)
            build_cmd = ["bowtie2-build", ref_fasta, index_prefix]
            print(f"[bwa_aln_mtdna] Building Bowtie2 mtDNA index:\n{' '.join(build_cmd)}")
            r = subprocess.run(build_cmd)
            if r.returncode != 0:
                if os.path.exists(lock_file):
                    os.remove(lock_file)
                raise RuntimeError(f"bowtie2-build (mtDNA) failed with exit code {r.returncode}")
            if os.path.exists(lock_file):
                os.remove(lock_file)
            if not index_complete():
                raise RuntimeError(
                    f"bowtie2-build finished but mtDNA index still incomplete at prefix {index_prefix!r}"
                )
            return
        except FileExistsError:
            if index_complete():
                return
            if time.time() - start > timeout_sec:
                raise RuntimeError(
                    f"Timeout waiting for Bowtie2 mtDNA index at {index_prefix} (lock {lock_file})"
                )
            time.sleep(10)


fallback = "Sheep"

if host_aligner == "bowtie2":
    bt2_cfg = bowtie2_indices.get(species) or bowtie2_indices.get(fallback)
    fa_cfg = bwa_indices.get(species) or bwa_indices.get(fallback)
    if not bt2_cfg and not fa_cfg:
        raise ValueError(
            f"No bowtie2_mtDNA_indices or mtDNA_indices entry for species '{species}' "
            f"(or fallback '{fallback}') with host_aligner=bowtie2."
        )
    index_prefix, ref_fasta = resolve_bowtie2_prefix_and_fasta(bt2_cfg, fa_cfg)
else:
    fa_cfg = bwa_indices.get(species) or bwa_indices.get(fallback)
    if not fa_cfg:
        raise ValueError(
            f"No mtDNA_indices entry for species '{species}' (or fallback '{fallback}') with host_aligner=bwa."
        )
    index_prefix = resolve_bwa_database_prefix(fa_cfg)
    ref_fasta = None

# If using bowtie2, ensure index is complete (not just .1.bt2 present)
if host_aligner == "bowtie2":
    ensure_bowtie2_index(index_prefix, ref_fasta)
    inspect = subprocess.run(
        ["bowtie2-inspect", "-s", index_prefix],
        capture_output=True,
        text=True,
    )
    if inspect.returncode != 0:
        raise RuntimeError(
            f"bowtie2-inspect failed for mtDNA index prefix {index_prefix!r} "
            f"(index may be corrupt or built with another Bowtie2 version). "
            f"stderr: {inspect.stderr.strip()}"
        )

# Ensure output dir exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)


def run_bowtie2_pipe_to_bam(
    bowtie2_cmd: list,
    log_file: str,
    sam_threads: int,
    out_bam: str,
) -> None:
    """
    Run bowtie2 | samtools view without a shell pipeline so both exit codes are checked.
    (Default /bin/sh does not use pipefail; bowtie2 can fail while samtools sees EOF.)
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
    # Bowtie2 mtDNA mapping with read group
    log_file = output_file.replace(".bam", "_bowtie2.log")
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
    print(f"Running bowtie2 -> samtools (log: {log_file})")
    run_bowtie2_pipe_to_bam(bowtie2_cmd, log_file, threads, output_file)
    cmd = None
else:
    # BWA mtDNA mapping (original behaviour)
    sai_file = output_file.replace(".bam", ".sai")
    cmd = (
        f"bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {index_prefix} {reads_file} > {sai_file} && "
        f"bwa samse -r '@RG\\tID:{sample}_host\\tSM:{sample}\\tPL:ILLUMINA' {index_prefix} {sai_file} {reads_file} | "
        f"samtools view -@ {threads} -F 4 -b -o {output_file} -"
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
            f"Host/mtDNA mapping failed with exit code {e.returncode}"
        ) from e

