import os
import sys
import subprocess
import yaml


def main():
    if "snakemake" not in globals():
        print("This script must be run by Snakemake.", file=sys.stderr)
        sys.exit(2)

    bam = snakemake.input.get("bam")
    species_file = snakemake.input.get("species_file")
    cfg_path = snakemake.input.get("config") or snakemake.params.get("config")
    out_dir = snakemake.output.get("dir") if isinstance(snakemake.output, dict) else snakemake.output[0]
    log_path = snakemake.log[0] if snakemake.log else None
    index_key = snakemake.params.get("index_key")

    os.makedirs(out_dir, exist_ok=True)

    # Read species robustly
    try:
        with open(species_file, "r", encoding="utf-8", errors="ignore") as f:
            species = f.read().strip()
    except Exception as e:
        msg = f"Failed to read species file '{species_file}': {e}"
        print(msg, file=sys.stderr)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
        sys.exit(1)

    # Load config and resolve reference
    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f)
    ref_map = cfg.get(index_key, {}) or {}
    ref_path = ref_map.get(species, "")

    if not ref_path:
        msg = f"Reference for species '{species}' not found under '{index_key}' in {cfg_path}"
        print(msg, file=sys.stderr)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
        sys.exit(1)

    # Sanity checks
    if not os.path.exists(bam) or os.path.getsize(bam) == 0:
        msg = f"Missing or empty BAM: {bam}"
        print(msg, file=sys.stderr)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
        sys.exit(1)
    if not os.path.exists(ref_path) or os.path.getsize(ref_path) == 0:
        msg = f"Missing or empty reference: {ref_path}"
        print(msg, file=sys.stderr)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
        sys.exit(1)

    # Ensure FASTA index exists for tools that may expect it
    fai_path = ref_path + ".fai"
    if not os.path.exists(fai_path):
        try:
            subprocess.run(["samtools", "faidx", ref_path], check=True)
        except Exception as e:
            msg = f"Failed to index reference with samtools faidx: {e}"
            print(msg, file=sys.stderr)
            if log_path:
                with open(log_path, "a") as lf:
                    lf.write(msg + "\n")
            sys.exit(1)

    if log_path:
        with open(log_path, "a") as lf:
            lf.write(f"[run_damageprofiler] species={species} ref={ref_path}\n")

    # If BAM has zero mapped reads, create placeholder outputs and exit successfully
    try:
        mapped_proc = subprocess.run(["samtools", "view", "-c", "-F", "4", bam], capture_output=True, text=True, check=True)
        mapped_count = int(mapped_proc.stdout.strip() or 0)
    except Exception as e:
        msg = f"Failed to count mapped reads in BAM: {e}"
        print(msg, file=sys.stderr)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write(msg + "\n")
        # Be conservative: if we cannot count, continue to try DP
        mapped_count = None

    if mapped_count == 0:
        note_path = os.path.join(out_dir, "NO_MAPPED_READS.txt")
        try:
            with open(note_path, "w") as nf:
                nf.write(f"No mapped reads found in {bam}. DamageProfiler skipped.\n")
            # Create placeholder file expected by downstream, if any
            placeholder = os.path.join(out_dir, "misincorporation.txt")
            with open(placeholder, "w") as pf:
                pf.write("")
        except Exception as e:
            print(f"Failed to write placeholder outputs: {e}", file=sys.stderr)
            sys.exit(1)
        if log_path:
            with open(log_path, "a") as lf:
                lf.write("No mapped reads. Skipping DamageProfiler and creating placeholders.\n")
        sys.exit(0)

    cmd = [
        "damageprofiler",
        "-Xmx12G",
        "-i", bam,
        "-o", out_dir,
        "-r", ref_path,
    ]

    # Run and stream output to log if provided
    with open(log_path, "a") if log_path else open(os.devnull, "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=lf)
    if proc.returncode != 0:
        print(f"DamageProfiler failed with exit code {proc.returncode}", file=sys.stderr)
        sys.exit(proc.returncode)


if __name__ == "__main__":
    main()


