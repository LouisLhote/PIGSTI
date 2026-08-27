#!/usr/bin/env python3
"""
Build .idx from host BAM idxstats and run residual-method sexing (R) for supported species.

Supported FastQ Screen best-species labels: Cow, Goat, Sheep, Dog.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

from sexing_utils import resolve_bio_host_species

# FastQ Screen species -> optional expected autosome count (QC hint only).
# X is detected from reference names in idxstats; R script infers N from .idx rows.
SEXING_KITS = {
    "Cow": {"expected_n_auto": 29, "r_script": "sexing_residual_method.R"},
    "Goat": {"expected_n_auto": 29, "r_script": "sexing_residual_method.R"},
    "Sheep": {"expected_n_auto": 26, "r_script": "sexing_residual_method.R"},
    "Dog": {"expected_n_auto": 38, "r_script": "sexing_residual_method.R"},
}

SUPPORTED_SPECIES = set(SEXING_KITS.keys())


def _parse_chr_ref(ref_name: str) -> tuple[str, int | None]:
    """Return ('auto', n) or ('x', None) or ('skip', None)."""
    base = ref_name.split()[0].strip()
    if base == "*" or not base:
        return "skip", None
    low = base.lower()
    if low.startswith("chr"):
        base = base[3:]
    if base.upper() in ("X", "CHRX") or re.fullmatch(r"[xX]", base):
        return "x", None
    if base.upper() in ("Y", "CHRY") or re.fullmatch(r"[yY]", base):
        return "y", None
    m = re.fullmatch(r"0*([0-9]+)", base)
    if m:
        return "auto", int(m.group(1))
    return "skip", None


def _run_idxstats(bam_path: str) -> list[tuple[str, int, int]]:
    proc = subprocess.run(
        ["samtools", "idxstats", bam_path],
        check=True,
        capture_output=True,
        text=True,
    )
    rows: list[tuple[str, int, int]] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        ref, length_s, mapped_s, _unmapped_s = parts[0], parts[1], parts[2], parts[3]
        try:
            length = int(length_s)
            mapped = int(mapped_s)
        except ValueError:
            continue
        rows.append((ref, length, mapped))
    return rows


def build_idx_matrix(
    idxstats_rows: list[tuple[str, int, int]],
    expected_n_auto: int | None,
    sample_id: str,
) -> tuple[list[str], list[int], list[int], str | None]:
    """
    Build .idx rows: autosomes 1..N (from reference) then X.

    Finds X by name (chrX, X, …). Y is ignored. Returns optional QC warning string.
    """
    autosomes: dict[int, tuple[int, int]] = {}
    x_candidates: list[tuple[int, int, str]] = []

    for ref, length, mapped in idxstats_rows:
        kind, num = _parse_chr_ref(ref)
        if kind == "skip" or kind == "y":
            continue
        if kind == "auto" and num is not None:
            prev = autosomes.get(num)
            if prev is None or mapped > prev[1]:
                autosomes[num] = (length, mapped)
        elif kind == "x":
            x_candidates.append((length, mapped, ref.split()[0]))

    if not autosomes:
        raise ValueError(
            "No numbered autosomes in idxstats (expected contigs named 1, 2, … or chr1, chr2, …)"
        )
    if not x_candidates:
        raise ValueError(
            "No X chromosome in idxstats (expected contig X or chrX; Y is ignored)"
        )

    x_length, x_mapped, x_label = max(x_candidates, key=lambda t: t[1])
    auto_nums = sorted(autosomes.keys())
    n_auto = len(auto_nums)

    # Prefer contiguous 1..max (typical reference builds)
    max_auto = auto_nums[-1]
    missing = [n for n in range(1, max_auto + 1) if n not in autosomes]
    if missing:
        raise ValueError(
            f"Non-contiguous autosome numbering in reference (missing 1..{max_auto}): "
            f"{missing[:12]}{'...' if len(missing) > 12 else ''}"
        )
    auto_nums = list(range(1, max_auto + 1))
    n_auto = len(auto_nums)

    warning = None
    if expected_n_auto is not None and n_auto != expected_n_auto:
        warning = (
            f"autosome_count_mismatch: expected {expected_n_auto} for species, "
            f"found {n_auto} numbered autosomes (1..{max_auto})"
        )

    chrs: list[str] = []
    lengths: list[int] = []
    reads: list[int] = []
    for n in auto_nums:
        length, mapped = autosomes[n]
        chrs.append(str(n))
        lengths.append(length)
        reads.append(mapped)
    chrs.append(x_label if x_label else "X")
    lengths.append(x_length)
    reads.append(x_mapped)
    return chrs, lengths, reads, warning


def write_idx(path: str, chrs: list[str], lengths: list[int], reads: list[int], sample_id: str) -> None:
    lines = ["chr\tlength\t" + sample_id]
    for c, ln, rd in zip(chrs, lengths, reads):
        lines.append(f"{c}\t{ln}\t{rd}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_skip_outputs(
    idx_path: str,
    pdf_path: str,
    tsv_path: str,
    sample_id: str,
    species: str,
    reason: str,
) -> None:
    write_idx(idx_path, ["1"], [0], [0], sample_id)
    Path(tsv_path).write_text(
        "sample\tspecies\tstatus\tcall\tfemale_prob\tlikelihood_ratio\tnote\n"
        f"{sample_id}\t{species}\tskipped\tNA\tNA\tNA\t{reason}\n",
        encoding="utf-8",
    )
    # Minimal valid PDF
    subprocess.run(
        [
            "Rscript",
            "-e",
            f'pdf("{pdf_path}"); plot.new(); title(main="Sexing skipped: {reason}"); dev.off()',
        ],
        check=True,
    )


def main() -> int:
    if "snakemake" not in globals():
        print("run_pigsti_sexing.py must be run via Snakemake", file=sys.stderr)
        return 1

    bam = snakemake.input.bam
    species_files = list(getattr(snakemake.input, "species_files", []) or [])
    mismatch_warning = getattr(snakemake.input, "species_mismatch_warning", None)
    idx_out = snakemake.output.idx
    pdf_out = snakemake.output.pdf
    tsv_out = snakemake.output.tsv
    sample = snakemake.wildcards.sample  # biological sample ID
    enabled = int(snakemake.params.enabled)
    host_ref_map = json.loads(snakemake.params.host_ref_map)
    r_script = Path(snakemake.input.r_script).resolve()

    os.makedirs(os.path.dirname(idx_out), exist_ok=True)

    species, species_note = resolve_bio_host_species(species_files, mismatch_warning)

    if not enabled:
        write_skip_outputs(idx_out, pdf_out, tsv_out, sample, species, "enable_sexing=false")
        return 0

    if species_note:
        write_skip_outputs(idx_out, pdf_out, tsv_out, sample, species, species_note)
        return 0

    if species not in SUPPORTED_SPECIES:
        write_skip_outputs(
            idx_out,
            pdf_out,
            tsv_out,
            sample,
            species,
            f"species_not_supported ({species})",
        )
        return 0

    if not host_ref_map.get(species):
        write_skip_outputs(
            idx_out,
            pdf_out,
            tsv_out,
            sample,
            species,
            "no_host_reference_in_config",
        )
        return 0

    if not os.path.isfile(bam) or os.path.getsize(bam) == 0:
        write_skip_outputs(idx_out, pdf_out, tsv_out, sample, species, "empty_host_bam")
        return 0

    kit = SEXING_KITS[species]
    try:
        stats = _run_idxstats(bam)
        chrs, lengths, reads, idx_warning = build_idx_matrix(
            stats, kit.get("expected_n_auto"), sample
        )
        write_idx(idx_out, chrs, lengths, reads, sample)
        if idx_warning:
            print(f"[sexing] {sample}: {idx_warning}", file=sys.stderr)
    except (subprocess.CalledProcessError, ValueError) as exc:
        write_skip_outputs(idx_out, pdf_out, tsv_out, sample, species, str(exc))
        return 0

    if not r_script.is_file():
        write_skip_outputs(
            idx_out,
            pdf_out,
            tsv_out,
            sample,
            species,
            f"missing_r_script:{r_script} (expected scripts/sexing_residual_method.R in repo)",
        )
        return 0

    try:
        subprocess.run(
            ["Rscript", str(r_script), idx_out, pdf_out, tsv_out],
            check=True,
        )
        # Annotate TSV
        lines = Path(tsv_out).read_text(encoding="utf-8").strip().splitlines()
        if lines and "status" not in lines[0]:
            out_lines = [lines[0] + "\tspecies\tstatus"]
            for row in lines[1:]:
                out_lines.append(row + f"\t{species}\tok")
            Path(tsv_out).write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    except subprocess.CalledProcessError as exc:
        write_skip_outputs(idx_out, pdf_out, tsv_out, sample, species, f"r_script_failed:{exc}")
    return 0


if "snakemake" in globals():
    raise SystemExit(main())
