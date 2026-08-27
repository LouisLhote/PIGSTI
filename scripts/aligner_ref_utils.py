"""
Resolve reference paths for host / mtDNA alignment.

Bowtie2 mode
  - `bowtie2_*` config may be either:
    - a complete Bowtie2 index prefix (six .bt2 files), or
    - a path to a FASTA (.fa / .fasta / .fna); prefix = basename without extension.
  - A companion FASTA from `bwa_indices` / `mtDNA_indices` is still used when the
    Bowtie2 prefix has no sibling FASTA and the index must be built.

BWA mode
  - Config value should point to a FASTA (or resolvable to one). We ensure
    `bwa index` has been run (standard `ref.fa.bwt` naming).
"""

from __future__ import annotations

import os
import subprocess


def bowtie2_index_complete(prefix: str) -> bool:
    if not prefix:
        return False
    return all(
        os.path.exists(prefix + suf)
        for suf in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    )


def pick_fasta_near_path(p: str) -> str:
    """Return an existing non-empty FASTA path: `p` if it is a FASTA file, else try p+.fa etc."""
    if not p:
        return ""
    p = p.strip()
    if os.path.isfile(p) and os.path.getsize(p) > 0:
        low = p.lower()
        if low.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
            return p
    base = p.rstrip("/")
    for suf in (".fa", ".fasta", ".fna"):
        c = base + suf
        if os.path.isfile(c) and os.path.getsize(c) > 0:
            return c
    return ""


def default_bowtie2_prefix_from_fasta(fa_path: str) -> str:
    low = fa_path.lower()
    for ext in (".fa.gz", ".fasta.gz", ".fna.gz", ".fasta", ".fna", ".fa"):
        if low.endswith(ext):
            return fa_path[: -len(ext)]
    return os.path.splitext(fa_path)[0]


def resolve_bowtie2_prefix_and_fasta(
    bowtie2_configured: str | None,
    fasta_fallback: str | None,
) -> tuple[str, str]:
    """
    Returns (bowtie2_index_prefix, ref_fasta) for `bowtie2 -x` and optional `bowtie2-build`.

    ref_fasta is always a concrete FASTA path when returned (needed for builds / sanity).
    """
    b = (bowtie2_configured or "").strip()
    r = (fasta_fallback or "").strip()

    if not b and not r:
        raise ValueError("resolve_bowtie2_prefix_and_fasta: no bowtie2 or FASTA path configured.")

    # 1) Explicit Bowtie2 prefix with a complete index
    if b and bowtie2_index_complete(b):
        fa = pick_fasta_near_path(r) or pick_fasta_near_path(b)
        if not fa:
            raise ValueError(
                f"Bowtie2 index is complete at {b!r}, but no FASTA found under bowtie2 path or "
                f"bwa_indices/mtDNA_indices (need .fa sibling or explicit FASTA in config)."
            )
        return b, fa

    # 2) bowtie2 column is a FASTA file → prefix next to / derived from FASTA
    if b and os.path.isfile(b) and os.path.getsize(b) > 0:
        low = b.lower()
        if low.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
            return default_bowtie2_prefix_from_fasta(b), b

    # 3) Incomplete prefix in bowtie2_indices; need FASTA from mtDNA/host indices
    fa = pick_fasta_near_path(r) or pick_fasta_near_path(b)
    if b and fa:
        return b, fa

    # 4) No bowtie2 path: derive prefix from FASTA-only config (single-column style)
    if not b and r:
        fa = pick_fasta_near_path(r)
        if fa:
            return default_bowtie2_prefix_from_fasta(fa), fa

    raise ValueError(
        f"Cannot resolve Bowtie2 prefix + FASTA from bowtie2={bowtie2_configured!r}, "
        f"fasta_fallback={fasta_fallback!r}"
    )


def bwa_index_present(fa_path: str) -> bool:
    # `bwa index ref.fa` creates ref.fa.bwt (among others)
    bwt = fa_path + ".bwt"
    return os.path.isfile(bwt) and os.path.getsize(bwt) > 0


def ensure_bwa_index(fa_path: str) -> None:
    if bwa_index_present(fa_path):
        return
    subprocess.run(["bwa", "index", fa_path], check=True)


def resolve_bwa_database_prefix(configured: str) -> str:
    """
    Return the BWA database prefix path to pass to `bwa aln` / `bwa samse`
    (same convention as `bwa index ref.fa` → use `ref.fa` as prefix).
    """
    p = (configured or "").strip()
    if not p:
        raise ValueError("resolve_bwa_database_prefix: empty path.")

    if os.path.isfile(p) and os.path.getsize(p) > 0:
        low = p.lower()
        if low.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
            ensure_bwa_index(p)
            return p

    fa = pick_fasta_near_path(p)
    if fa:
        ensure_bwa_index(fa)
        return fa

    raise ValueError(
        f"BWA mode: cannot resolve {p!r} to a FASTA file (expected .fa/.fasta/.fna or prefix with sibling .fa)."
    )
