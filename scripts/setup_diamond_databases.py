#!/usr/bin/env python3
"""
Standalone DIAMOND DB builder for PIGSTI.

RefSeq FTP note:
NCBI's RefSeq FTP no longer provides single monolithic files like
`bacteria.protein.faa.gz`. Instead, protein FASTAs are split into many files
matching `<category>*.protein.faa.gz`. This script supports that layout by
concatenating all matching split files found in `--refseq-dir`.

Modes:
  - targeted: keeps only pathogen/ruminant-relevant proteins (smaller DB)
  - full_refseq: combines complete RefSeq categories (large DB)
"""

import argparse
import gzip
import glob
import os
import shutil
import subprocess

PATHOGEN_KEYWORDS = [
    "bacillus anthracis",
    "yersinia pestis",
    "mycobacterium tuberculosis",
    "mycobacterium bovis",
    "brucella",
    "salmonella",
    "clostridium",
    "echinococcus",
    "toxoplasma gondii",
]

RUMINANT_KEYWORDS = [
    "bos taurus",
    "bos primigenius",
    "ovis aries",
    "capra hircus",
    "cervus",
]

MT_PROTEIN_KEYWORDS = [
    "cytochrome c oxidase subunit",
    "cytochrome b",
    "nad",
    "nadh dehydrogenase",
    "atp synthase",
]


def header_matches(header: str, patterns):
    h = header.lower()
    return any(p in h for p in patterns)


def get_split_protein_fastas(refseq_dir: str, category: str):
    """
    Return sorted list of split protein FASTA archives for a RefSeq category.
    Matches files like: <category>*.protein.faa.gz (including wp_protein parts).
    """
    pat = os.path.join(refseq_dir, f"{category}*.protein.faa.gz")
    return sorted(glob.glob(pat))


def download_splits_with_wget(refseq_dir: str, category: str):
    """
    Download all <category>*.protein.faa.gz files from the RefSeq release directory.
    This uses wget to avoid Python SSL issues on locked-down HPC nodes.
    """
    wget = shutil.which("wget")
    if not wget:
        raise RuntimeError("wget not found in PATH; either install wget or use --skip-download.")

    base = f"https://ftp.ncbi.nlm.nih.gov/refseq/release/{category}/"
    pat = f"{category}*.protein.faa.gz"

    print(f"[download] {category} protein splits via wget pattern: {pat}")
    cwd = os.getcwd()
    try:
        os.makedirs(refseq_dir, exist_ok=True)
        os.chdir(refseq_dir)
        # -r/--recursive is needed because the server is a directory listing with many files.
        # --no-parent prevents traversing upwards.
        subprocess.check_call(
            [
                "wget",
                "-c",
                "-r",
                "-np",
                "-nd",
                "-A",
                pat,
                "--reject",
                "*gpff.gz",
                base,
            ]
        )
    finally:
        os.chdir(cwd)


def filter_fasta_gz_to_handle(fin_path: str, out_handle, include_patterns):
    """Stream-parse a gz FASTA and write only matching sequences to an open file handle."""
    kept = 0
    with gzip.open(fin_path, "rt", encoding="utf-8", errors="ignore") as fin:
        keep = False
        for line in fin:
            if line.startswith(">"):
                keep = header_matches(line, include_patterns)
                if keep:
                    kept += 1
                    out_handle.write(line)
            elif keep:
                out_handle.write(line)
    return kept


def append_labeled_fasta_gz_to_out(in_gz: str, out_handle, label: str):
    """Append labeled FASTA sequences from a gz file to an open output handle."""
    count = 0
    with gzip.open(in_gz, "rt", encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            if line.startswith(">"):
                count += 1
                line = f">{label}|{line[1:]}"
            out_handle.write(line)
    return count


def build_targeted(refseq_dir, outdir):
    bacteria_parts = get_split_protein_fastas(refseq_dir, "bacteria")
    viral_parts = get_split_protein_fastas(refseq_dir, "viral")
    mito_parts = get_split_protein_fastas(refseq_dir, "mitochondrion")
    if not bacteria_parts or not viral_parts or not mito_parts:
        raise FileNotFoundError(
            "Missing split protein FASTA parts. Ensure --refseq-dir contains "
            "<category>*.protein.faa.gz files (download them with the bash script)."
        )

    pathogen_bac_faa = os.path.join(outdir, "pathogen_bacteria_refseq.faa")
    pathogen_viral_faa = os.path.join(outdir, "pathogen_viral_refseq.faa")
    host_mito_faa = os.path.join(outdir, "ruminant_mito_refseq.faa")
    pathogen_combined_faa = os.path.join(outdir, "pathogen_refseq.faa")
    combined_faa = os.path.join(outdir, "combined_ancient.faa")

    os.makedirs(outdir, exist_ok=True)

    kept_bac = 0
    with open(pathogen_bac_faa, "w", encoding="utf-8") as out_bac:
        for part in bacteria_parts:
            kept_bac += filter_fasta_gz_to_handle(part, out_bac, PATHOGEN_KEYWORDS)

    kept_vir = 0
    with open(pathogen_viral_faa, "w", encoding="utf-8") as out_vir:
        for part in viral_parts:
            kept_vir += filter_fasta_gz_to_handle(part, out_vir, PATHOGEN_KEYWORDS)

    kept_mito = 0
    with open(host_mito_faa, "w", encoding="utf-8") as out_mito:
        for part in mito_parts:
            kept_mito += filter_fasta_gz_to_handle(
                part, out_mito, RUMINANT_KEYWORDS + MT_PROTEIN_KEYWORDS
            )

    with open(pathogen_combined_faa, "w", encoding="utf-8") as out_pc:
        for p in [pathogen_bac_faa, pathogen_viral_faa]:
            with open(p, "r", encoding="utf-8") as inp:
                out_pc.write(inp.read())

    with open(combined_faa, "w", encoding="utf-8") as out_all:
        for p in [pathogen_combined_faa, host_mito_faa]:
            with open(p, "r", encoding="utf-8") as inp:
                out_all.write(inp.read())

    return {
        "combined_faa": combined_faa,
        "manifest_lines": [
            f"pathogen_bacteria_sequences\t{kept_bac}",
            f"pathogen_viral_sequences\t{kept_vir}",
            f"ruminant_mito_sequences\t{kept_mito}",
        ],
    }


def build_full_refseq(refseq_dir, outdir):
    categories = ["bacteria", "viral", "mitochondrion", "fungi", "archaea"]
    combined_faa = os.path.join(outdir, "full_refseq_combined.faa")
    manifest_lines = []

    os.makedirs(outdir, exist_ok=True)
    with open(combined_faa, "w", encoding="utf-8") as out_all:
        for cat in categories:
            parts = get_split_protein_fastas(refseq_dir, cat)
            if not parts:
                raise FileNotFoundError(f"Missing split protein FASTA parts for '{cat}' in {refseq_dir}")

            n_total = 0
            for part in parts:
                n_total += append_labeled_fasta_gz_to_out(part, out_all, label=cat.upper())
            manifest_lines.append(f"{cat}_sequences\t{n_total}")

    return {"combined_faa": combined_faa, "manifest_lines": manifest_lines}


def run(cmd: str):
    print(f"[run] {cmd}")
    subprocess.check_call(cmd, shell=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--db-name", default="combined_ancient_db")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--mode", choices=["targeted", "full_refseq"], default="full_refseq")
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.refseq_dir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)

    required = ["bacteria", "viral", "mitochondrion"]
    if args.mode == "full_refseq":
        required += ["fungi", "archaea"]

    if not args.skip_download:
        # Downloads split <category>*.protein.faa.gz files into --refseq-dir.
        for cat in required:
            download_splits_with_wget(args.refseq_dir, cat)

    # Ensure inputs exist (at least one split part per category).
    for cat in required:
        parts = get_split_protein_fastas(args.refseq_dir, cat)
        if not parts:
            raise FileNotFoundError(
                f"Missing required split protein FASTA parts for '{cat}' in {args.refseq_dir}. "
                f"Expected files like: {cat}*.protein.faa.gz"
            )

    built = build_targeted(args.refseq_dir, args.outdir) if args.mode == "targeted" else build_full_refseq(args.refseq_dir, args.outdir)

    db_prefix = os.path.join(args.outdir, args.db_name)
    run(
        f'diamond makedb --in "{built["combined_faa"]}" --db "{db_prefix}" --threads {args.threads}'
    )

    manifest = os.path.join(args.outdir, "diamond_db_manifest.txt")
    with open(manifest, "w", encoding="utf-8") as f:
        f.write(f"mode\t{args.mode}\n")
        f.write(f"combined_faa\t{built['combined_faa']}\n")
        f.write(f"diamond_db\t{db_prefix}.dmnd\n")
        for line in built["manifest_lines"]:
            f.write(line + "\n")

    print(f"[done] DB built: {db_prefix}.dmnd")
    print(f"[done] Manifest: {manifest}")


if __name__ == "__main__":
    main()
