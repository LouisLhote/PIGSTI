"""
PIGSTI v2/v3 canonical results paths.

Layout:
  results/libraries/{pcr}/       PCR-level QC, host_mapping/, mtdna_mapping/
  results/samples/{bio}/         Merged host/mtDNA sample QC (qualimap, etc.)
  results/pools/                 Merged unaligned FASTQs (Kraken/decOM input)
  results/host/                  Merged host and mtDNA BAMs per biological sample
  results/metagenomics/          KrakenUniq per sample, HOPS, decOM, abundance matrices
  results/pathogen/{bio}/        E-value tables, pathogen mapping, summaries, comparison
  results/final/                 Cohort Excel, heatmaps, manifest, catalog
  results/workflow/              Checkpoints and barriers
"""

from __future__ import annotations

from pigsti_naming import safe_pathogen_name

# ---------------------------------------------------------------------------
# PCR / library level
# ---------------------------------------------------------------------------

def library_root(pcr: str) -> str:
    return f"results/libraries/{pcr}"


def library_collapsed(pcr: str) -> str:
    return f"results/libraries/{pcr}/adapter_removal/{pcr}.collapsed.gz"


def library_prinseq_passed(pcr: str) -> str:
    return f"results/libraries/{pcr}/prinseq/{pcr}-passed.fq.gz"


def sexing_r_script() -> str:
    """Repo-relative path to unified sexing R script (scripts/ or scripts/sexing/)."""
    import os

    for rel in (
        "scripts/sexing/sexing_residual_method.R",
        "scripts/sexing_residual_method.R",
    ):
        if os.path.isfile(rel):
            return rel
    return "scripts/sexing_residual_method.R"


def library_fastq_screen_html(pcr: str) -> str:
    return f"results/libraries/{pcr}/fastq_screen/{pcr}.collapsed_screen.html"


def library_host_dedup_bam(pcr: str) -> str:
    return f"results/libraries/{pcr}/host_mapping/{pcr}.dedup.bam"


def library_mtdna_dedup_bam(pcr: str) -> str:
    return f"results/libraries/{pcr}/mtdna_mapping/{pcr}.dedup.bam"


def library_qualimap_host(pcr: str) -> str:
    return f"results/libraries/{pcr}/qualimap/genome_results.txt"


def library_qualimap_mtdna(pcr: str) -> str:
    return f"results/libraries/{pcr}/qualimap_mtdna/genome_results.txt"


def library_damage_host(pcr: str) -> str:
    return f"results/libraries/{pcr}/damageprofiler_host"


def library_damage_mtdna(pcr: str) -> str:
    return f"results/libraries/{pcr}/damageprofiler_mtdna"


def library_pathogen_mapping_dir(pcr: str) -> str:
    return f"results/libraries/{pcr}/pathogen_mapping"


# ---------------------------------------------------------------------------
# Biological sample — merged host / mtDNA QC
# ---------------------------------------------------------------------------

def sample_root(bio: str) -> str:
    return f"results/samples/{bio}"


def sample_qualimap_host(bio: str) -> str:
    return f"results/samples/{bio}/qualimap/genome_results.txt"


def sample_qualimap_mtdna(bio: str) -> str:
    return f"results/samples/{bio}/qualimap_mtdna/genome_results.txt"


# ---------------------------------------------------------------------------
# Metagenomics (per bio sample + cohort)
# ---------------------------------------------------------------------------

def metagenomics_kraken_report(bio: str) -> str:
    return f"results/metagenomics/krakenuniq/{bio}/{bio}_kraken-report.txt"


def metagenomics_kraken_output(bio: str) -> str:
    return f"results/metagenomics/krakenuniq/{bio}/{bio}_output.txt"


def metagenomics_hops_heatmap() -> str:
    return "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"


def metagenomics_decom_dir() -> str:
    """decOM tool output directory (inputs live in results/metagenomics/decOM/p_keys/)."""
    return "results/metagenomics/decOM/decOM_out"


def metagenomics_kraken_matrix_abs() -> str:
    return "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv"


def metagenomics_kraken_matrix_norm() -> str:
    return "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_normalized.csv"


def metagenomics_pathogen_detection_matrix() -> str:
    return "results/metagenomics/pathogen_detection/detection_scores_matrix.csv"


# ---------------------------------------------------------------------------
# Pathogen (per bio sample)
# ---------------------------------------------------------------------------

def pathogen_root(bio: str) -> str:
    return f"results/pathogen/{bio}"


def pathogen_evalue_genus(bio: str) -> str:
    return f"results/pathogen/{bio}/evalue/genus/{bio}_genus.csv"


def pathogen_evalue_species(bio: str) -> str:
    return f"results/pathogen/{bio}/evalue/species/{bio}_species.csv"


def pathogen_evalue_pathogen(bio: str) -> str:
    return f"results/pathogen/{bio}/evalue/pathogen/{bio}_pathogen.csv"


def pathogen_mapping_dir(bio: str) -> str:
    return f"results/pathogen/{bio}/pathogen_mapping"


def pathogen_mapping_bam(bio: str, pathogen: str) -> str:
    ps = safe_pathogen_name(pathogen)
    return f"results/pathogen/{bio}/pathogen_mapping/{bio}_{ps}.dedup.bam"


def pathogen_mapping_bam_safe(bio: str, pathogen_safe: str) -> str:
    return f"results/pathogen/{bio}/pathogen_mapping/{bio}_{pathogen_safe}.dedup.bam"


def pathogen_qualimap(bio: str, pathogen_safe: str) -> str:
    return f"results/pathogen/{bio}/pathogen_mapping/qualimap_{pathogen_safe}"


def pathogen_damage(bio: str, pathogen_safe: str) -> str:
    return f"results/pathogen/{bio}/pathogen_mapping/damageprofiler_{pathogen_safe}"


def pathogen_metric(bio: str, pathogen_safe: str, suffix: str) -> str:
    return f"results/pathogen/{bio}/pathogen_mapping/{bio}_{pathogen_safe}.{suffix}"


def pathogen_entropy_100bp(bio: str, pathogen_safe: str) -> str:
    return pathogen_metric(bio, pathogen_safe, "entropy_100bp.txt")


def pathogen_entropy_1000bp(bio: str, pathogen_safe: str) -> str:
    return pathogen_metric(bio, pathogen_safe, "entropy_1000bp.txt")


def pathogen_report_pdf(bio: str, pathogen_safe: str) -> str:
    return f"results/pathogen/{bio}/summary/{bio}_{pathogen_safe}_pathogen_report.pdf"


def pathogen_summary_csv(bio: str) -> str:
    return f"results/pathogen/{bio}/summary/{bio}_pathogen_summary.csv"


def pathogen_sample_report_pdf(bio: str) -> str:
    return f"results/pathogen/{bio}/summary/{bio}_sample_report.pdf"


def pathogen_comparison_tsv(bio: str) -> str:
    return f"results/pathogen/{bio}/comparison/{bio}_comparison.tsv"


def pathogen_comparison_html(bio: str) -> str:
    return f"results/pathogen/{bio}/comparison/{bio}_heatmap.html"


# ---------------------------------------------------------------------------
# Pools (unaligned metagenomics input) + host lane (merged host/mtDNA)
# ---------------------------------------------------------------------------

def pools_unaligned_merged(bio: str) -> str:
    return f"results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz"


def host_host_merged_bam(bio: str) -> str:
    return f"results/host/host_mapping/{bio}.dedup.merged.bam"


def host_host_merged_metrics(bio: str) -> str:
    return f"results/host/host_mapping/{bio}.dedup.merged.metrics.txt"


def host_host_species_warning(bio: str) -> str:
    return f"results/host/host_mapping/{bio}.species_mismatch_warning.txt"


def host_mtdna_merged_bam(bio: str) -> str:
    return f"results/host/mtdna_mapping/{bio}.dedup.merged.bam"


def host_mtdna_merged_metrics(bio: str) -> str:
    return f"results/host/mtdna_mapping/{bio}.dedup.merged.metrics.txt"


def host_mtdna_species_warning(bio: str) -> str:
    return f"results/host/mtdna_mapping/{bio}.species_mismatch_warning.txt"


# Aliases (older names)
def pools_host_merged_bam(bio: str) -> str:
    return host_host_merged_bam(bio)


def pools_mtdna_merged_bam(bio: str) -> str:
    return host_mtdna_merged_bam(bio)


def pools_pathogen_mapping_merged_bam(bio: str, pathogen_safe: str) -> str:
    return f"results/pools/pathogen_mapping/{bio}_{pathogen_safe}.dedup.merged.bam"


# ---------------------------------------------------------------------------
# Final / workflow
# ---------------------------------------------------------------------------

def final_pathogen_summary_xlsx() -> str:
    return "results/final/pathogen_summary_all_samples.xlsx"


def final_host_mtdna_xlsx() -> str:
    return "results/final/host_mtdna_summary_all_samples.xlsx"


def final_comprehensive_xlsx() -> str:
    return "results/final/comprehensive_summary_all_samples.xlsx"


def final_output_catalog() -> str:
    return "results/final/output_catalog.tsv"


def workflow_pathogen_targets_txt() -> str:
    return "results/workflow/pathogen_targets.txt"


def workflow_run_manifest() -> str:
    return "results/final/run_manifest.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sample_from_evalue_pathogen_csv(path: str) -> str:
    path = str(path).replace("\\", "/")
    import re
    m = re.search(r"results/pathogen/([^/]+)/evalue/pathogen/", path)
    if not m:
        m = re.search(r"results/[^/]+/([^/]+)/evalue/pathogen/", path)
    if not m:
        raise ValueError(f"Cannot parse biological sample from E-value path: {path}")
    return m.group(1)


def parse_pathogen_targets_bam_line(line: str) -> tuple[str, str] | None:
    """
    Parse one line from ``results/workflow/pathogen_targets.txt``.

    v3: ``results/pathogen/{bio}/pathogen_mapping/{bio}_{pathogen_safe}.dedup.bam``
    legacy: ``results/{bio}/pathogen_mapping/{bio}_{pathogen_safe}.dedup.bam``
    """
    line = str(line).strip().replace("\\", "/")
    if not line or not line.endswith(".dedup.bam"):
        return None
    parts = line.split("/")
    if len(parts) >= 5 and parts[0] == "results" and parts[1] == "pathogen":
        sample = parts[2]
    elif len(parts) >= 4 and parts[0] == "results":
        sample = parts[1]
    else:
        return None
    bam_name = parts[-1]
    prefix = f"{sample}_"
    if not bam_name.startswith(prefix):
        return None
    pathogen_safe = bam_name[len(prefix) :]
    if pathogen_safe.endswith(".dedup.bam"):
        pathogen_safe = pathogen_safe[: -len(".dedup.bam")]
    return sample, pathogen_safe


def kraken_input_fastq_path(bio: str) -> str:
    return pools_unaligned_merged(bio)


# Aliases (older names)
def sample_kraken_report(bio: str) -> str:
    return metagenomics_kraken_report(bio)


def sample_evalue_pathogen(bio: str) -> str:
    return pathogen_evalue_pathogen(bio)


def sample_pathogen_mapping_bam(bio: str, pathogen: str) -> str:
    return pathogen_mapping_bam(bio, pathogen)


def kraken_report(bio_sample: str) -> str:
    return metagenomics_kraken_report(bio_sample)


def escore_pathogen_csv(bio_sample: str) -> str:
    return pathogen_evalue_pathogen(bio_sample)


def pathogen_bam(sample: str, pathogen: str) -> str:
    return pathogen_mapping_bam(sample, pathogen)
