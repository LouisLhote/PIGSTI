configfile: "config/config.yaml"

ruleorder: adapter_removal_pe > adapter_removal_se

global_wildcard_constraints = {
    "krakenuniq_name": r"[A-Za-z0-9_.-]+",
    "sample": r"[A-Za-z0-9_.-]+",
    "ref_name_safe": r"[A-Za-z0-9_]+"
}

import os
import sys
import pandas as pd
import re
import csv
from pathlib import Path
import json
import yaml
import collections
from snakemake.io import temp

# Load config once globally (needed by rules' params)
with open("config/config.yaml") as f:
    CFG = yaml.safe_load(f)


def _deep_merge_dicts(base, overlay):
    """
    Deep merge nested dicts for keys like `bwa_indices` / `mtDNA_indices`.
    - For dict values: merge recursively
    - For non-dict values: overlay wins
    """
    for k, v in (overlay or {}).items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            _deep_merge_dicts(base[k], v)
        else:
            base[k] = v
    return base


# Optional overlay config for "user friendly" profiles.
# You can set either:
# - config key: `config_overlay: <path-to-yaml>`
# - environment variable: `PIGSTI_CONFIG_OVERLAY=<path-to-yaml>`
_overlay_path = CFG.get("config_overlay") or os.environ.get("PIGSTI_CONFIG_OVERLAY")
if _overlay_path:
    if not os.path.exists(_overlay_path):
        raise FileNotFoundError(f"[Snakefile] config_overlay not found: {_overlay_path}")
    with open(_overlay_path) as _of:
        _overlay_cfg = yaml.safe_load(_of) or {}
    CFG = _deep_merge_dicts(CFG, _overlay_cfg)
    try:
        _deep_merge_dicts(config, _overlay_cfg)  # also update Snakemake's `config` dict
    except Exception:
        # If `config` isn't available (shouldn't happen in Snakemake), we just proceed with CFG.
        pass

# Snakemake CLI/profile overrides (--config key=value) update `config`, not CFG.
try:
    CFG = _deep_merge_dicts(CFG, dict(config))
except NameError:
    pass


def setup_results_root(cfg: dict) -> str:
    """
    Redirect pipeline output away from ./results when requested.

    All workflow paths still use the literal prefix ``results/`` in the Snakefile;
    at load time we replace ``./results`` with a symlink to ``results_root`` so
    existing rules and scripts keep working.

    Set via config.yaml, Snakemake CLI, or environment variable::

        snakemake --config results_root=/path/to/output
        export PIGSTI_RESULTS_ROOT=/path/to/output

    Relative paths are resolved from the directory where you launch Snakemake.
    """
    raw = (
        cfg.get("results_root")
        or cfg.get("results_dir")
        or os.environ.get("PIGSTI_RESULTS_ROOT")
        or "results"
    )
    target = os.path.abspath(os.path.expanduser(str(raw).strip()))
    link_name = os.path.abspath("results")

    if target == link_name:
        os.makedirs(target, exist_ok=True)
        resolved = target
    elif os.path.lexists(link_name):
        if os.path.islink(link_name):
            current = os.path.abspath(os.readlink(link_name))
            if current == target:
                resolved = target
            else:
                os.unlink(link_name)
                resolved = None
        elif os.listdir(link_name):
            raise ValueError(
                f"[PIGSTI] results_root is {target!r} but non-empty directory {link_name!r} "
                "already exists in the pipeline folder. Rename or move ./results, then retry:\n"
                f"  mv results results_local_backup && mkdir -p {target!r}\n"
                f"  # optional: rsync -a results_local_backup/ {target!r}/"
            )
        else:
            os.rmdir(link_name)
            resolved = None
    else:
        resolved = None

    if resolved is None:
        os.makedirs(os.path.dirname(target) or ".", exist_ok=True)
        os.makedirs(target, exist_ok=True)
        os.symlink(target, link_name, target_is_directory=True)
        print(f"[PIGSTI] results -> {target}", file=sys.stderr)
        resolved = target

    os.makedirs(".pigsti", exist_ok=True)
    Path(".pigsti/results_root").write_text(resolved + "\n", encoding="utf-8")
    actual = os.path.realpath(link_name)
    if actual != resolved:
        print(
            f"[PIGSTI] WARNING: results resolves to {actual!r}, expected {resolved!r}",
            file=sys.stderr,
        )
    return resolved


RESULTS_ROOT = setup_results_root(CFG)

# Duplicate / merge tool configuration
# - dedup_tool:  "picard" (default) or "samtools"
# - merge_tool:  "picard" (default) or "samtools"
# - markdup_threads: threads to use for duplicate marking (default: 6 for samtools)
DEDUP_TOOL = CFG.get("dedup_tool", "picard")
MERGE_TOOL = CFG.get("merge_tool", "picard")
MARKDUP_THREADS = int(CFG.get("markdup_threads", 6))

# FastQ Screen / Picard always come from the rule conda env (PATH), not config paths.
FASTQ_SCREEN_EXE = "fastq_screen"
# Re-run FastQ Screen on full FASTQ (--subset 0) when best #One_hit_one_genome is low.
FASTQ_SCREEN_FULL_RESCREEN = bool(CFG.get("fastq_screen_full_dataset_rescreen", True))
FASTQ_SCREEN_MIN_ONE_HIT = int(CFG.get("fastq_screen_full_dataset_min_one_hit", 50))
# Genetic sexing (residual method) for Cow, Goat, Sheep, Dog after host mapping.
ENABLE_SEXING = bool(CFG.get("enable_sexing", True))


def _resolve_sexing_r_script() -> str:
    """Unified R script under scripts/ or scripts/sexing/ (both layouts supported)."""
    for rel in (
        "scripts/sexing/sexing_residual_method.R",
        "scripts/sexing_residual_method.R",
    ):
        if os.path.isfile(rel):
            return rel
    return "scripts/sexing_residual_method.R"


SEXING_R_SCRIPT = _resolve_sexing_r_script()
PICARD_CMD = "picard"


# -------------------- Load sample & PCR info --------------------
# Backwards-compatible:
# - If a column "pcr" exists in the samples table, it is used as the PCR/library ID.
# - Otherwise, each row is treated as its own PCR, with pcr_id == sample.

SAMPLES_DICT = {}  # per-PCR r1/r2 lists (keyed by pcr_id)
PCRS = []          # list of PCR IDs
PCR_INFO = {}      # pcr_id -> metadata dict
SAMPLE_TO_PCRS = collections.defaultdict(list)  # biological sample -> list of PCR IDs

SAMPLES_TSV = CFG.get("samples", "config/samples.tsv")

with open(SAMPLES_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    # Normalize header names to avoid issues with BOMs or stray whitespace
    if reader.fieldnames is not None:
        reader.fieldnames = [fn.strip().lstrip("\ufeff") for fn in reader.fieldnames]
    has_pcr_col = "pcr" in reader.fieldnames if reader.fieldnames is not None else False

    for row in reader:
        # Skip completely empty / malformed rows
        if not row or "sample" not in row or row["sample"] is None or not str(row["sample"]).strip():
            continue

        # Normalize values from TSV to avoid hidden whitespace/CRLF issues
        sample = str(row["sample"]).strip()
        pcr_id = str(row["pcr"]).strip() if has_pcr_col and row.get("pcr") else sample
        r1 = str(row["r1"]).strip() if row.get("r1") is not None else ""
        r2 = str(row["r2"]).strip() if row.get("r2") is not None else ""

        # Per-PCR read lists
        if pcr_id not in SAMPLES_DICT:
            SAMPLES_DICT[pcr_id] = {"sample": sample, "r1": [], "r2": []}
        SAMPLES_DICT[pcr_id]["r1"].append(r1)
        if r2:
            SAMPLES_DICT[pcr_id]["r2"].append(r2)

        # PCR metadata
        if pcr_id not in PCR_INFO:
            source = row.get("source", "").strip() if row.get("source") is not None else ""
            if not source:
                # Default: treat rows without an explicit source as LOCAL
                # ENA-derived rows written by ena_to_pigsti_samplesheet.py
                # will have source == "ENA".
                source = "LOCAL"
            PCR_INFO[pcr_id] = {
                "sample": sample,
                "r1": r1,
                "r2": r2,
                "RGLB": row.get("RGLB", ""),
                "sequencing_run": row.get("sequencing_run", ""),
                "source": source,
            }
            PCRS.append(pcr_id)
            SAMPLE_TO_PCRS[sample].append(pcr_id)

# For now, treat each PCR as a "sample" in the workflow for PCR-level rules.
SAMPLES = list(SAMPLES_DICT.keys())

# Biological samples (keys of SAMPLE_TO_PCRS) for sample-level rules (KrakenUniq, HOPS, summaries, etc.)
BIO_SAMPLES = list(SAMPLE_TO_PCRS.keys())

if not SAMPLES:
    raise ValueError(
        f"[PIGSTI] No PCR/libraries loaded from {SAMPLES_TSV}. "
        "Check the TSV has a header row and at least one valid sample row."
    )
if not BIO_SAMPLES:
    raise ValueError(
        f"[PIGSTI] No biological samples after parsing {SAMPLES_TSV}."
    )

def _clean_path_field(path_value):
    """Normalize TSV path cells (trim whitespace and wrapping quotes)."""
    if path_value is None:
        return ""
    p = str(path_value).strip()
    if len(p) >= 2 and ((p[0] == '"' and p[-1] == '"') or (p[0] == "'" and p[-1] == "'")):
        p = p[1:-1].strip()
    return p

def _resolve_existing_fastq_path(path_value, sample_id, mate_label):
    """
    Return a normalized path and verify it exists.
    If the provided extension is .fastq.gz/.fq.gz and the file is missing,
    try the alternative extension before failing.
    """
    p = _clean_path_field(path_value)
    if not p:
        return p
    if os.path.exists(p):
        return p

    alt = None
    if p.endswith(".fastq.gz"):
        alt = p[:-9] + ".fq.gz"
    elif p.endswith(".fq.gz"):
        alt = p[:-6] + ".fastq.gz"
    if alt and os.path.exists(alt):
        return alt

    raise FileNotFoundError(
        f"[samples.tsv] Missing {mate_label} FASTQ for sample/PCR '{sample_id}'. "
        f"Normalized path: {repr(p)} "
        f"(length={len(p)})."
    )

BY_TYPE_ROOT = "results/browse"

# Pathogen mapping mode:
# - "default":      use Bowtie2-unmapped per-PCR reads for pathogen BWA mapping (faster)
# - "super_careful": use BWA host-unmapped per-PCR reads for pathogen BWA mapping (more stringent)
PATHOGEN_MODE = CFG.get("pathogen_mapping_mode", "super_careful")

# Pathogen aligner for mapping to pathogen references:
# - "bwa"      (current behaviour)
# - "bowtie2"  (use bowtie2 --end-to-end --sensitive)
PATHOGEN_ALIGNER = CFG.get("pathogen_aligner", "bwa")

# Optional HOPS step (for users who want to skip HOPS entirely)
ENABLE_HOPS = CFG.get("enable_hops", False)
_hops_params = CFG.get("hops_parameters") if isinstance(CFG.get("hops_parameters"), dict) else {}
HOPS_PARALLEL = bool(CFG.get("hops_parallel", False))
HOPS_MALT_MMAP = bool(CFG.get("hops_malt_mmap", False))
HOPS_PARALLEL_JOBS = max(1, int(CFG.get("hops_parallel_jobs", 2)))
_default_malt_threads = int(_hops_params.get("threadsMalt", 15))
HOPS_THREADS_PER_JOB = max(1, int(CFG.get("hops_threads_per_job", _default_malt_threads // HOPS_PARALLEL_JOBS)))
HOPS_HEAP_GB = max(1, int(CFG.get("hops_heap_gb", 800)))
HOPS_MALTEX_THREADS = max(1, int(_hops_params.get("threadsMaltEx", 15)))
HOPS_MALT_ROOT = "results/metagenomics/hops/malt"
HOPS_RMA_ROOT = "results/metagenomics/hops/rma"
HOPS_HEATMAP = "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"
# Limit concurrent MALT jobs at runtime: snakemake --resources hops_jobs=N
# (N defaults from hops_parallel_jobs in config.yaml; hops_malt_per_sample uses hops_jobs=1)


def hops_malt_done_path(bio: str) -> str:
    return f"{HOPS_MALT_ROOT}/{bio}/.malt_done"


def hops_staged_rma_path(bio: str) -> str:
    return f"{HOPS_RMA_ROOT}/{bio}_unaligned.rma6"


# Host/mtDNA summary parallelism (avoid threads: 999 on shared clusters)
SUMMARY_THREADS = int(CFG.get("summary_threads", 8))

# decOM: default False matches legacy ``decOM ... || true`` (pipeline continues on tool failure).
DECOM_FAIL_ON_ERROR = bool(CFG.get("decom_fail_on_error", False))
# decOM allocates approximately (decom_memory × decom_threads) RAM (legacy PIGSTI: 64GB × 8).
DECOM_MEMORY = str(CFG.get("decom_memory", "64GB"))
DECOM_THREADS = int(CFG.get("decom_threads", 8))

# Edit-distance damage-vs-no-damage split (HOPS-like) — always enabled.
EDIT_DISTANCE_DAMAGE_SPLIT = True
EDIT_DISTANCE_DAMAGE_WINDOW = int(CFG.get("edit_distance_damage_window_size", 5))

# Multi-QC HTML dashboard (per biological sample + pathogen info)
ENABLE_MULTI_QC_DASHBOARD = CFG.get("enable_multi_qc_dashboard", False)
MULTI_QC_TOP_PATHOGENS = int(CFG.get("multi_qc_top_pathogens", 5))

# Failure policy:
# - strict_inputs=True (default): fail fast on missing/empty critical inputs.
# - strict_inputs=False: preserve legacy fail-open behavior for some rules.
STRICT_INPUTS = bool(CFG.get("strict_inputs", True))

# End-of-run TSV catalog + optional results/browse/ symlinks (canonical paths unchanged).
BUILD_RESULTS_CATALOG = bool(CFG.get("build_results_catalog", True))
BUILD_RESULTS_CATALOG_SYMLINKS = bool(CFG.get("build_results_catalog_symlinks", True))

# Optional decOM step (off by default; opt in via config.yaml)
ENABLE_DECOM = CFG.get("enable_decom", False)

# Pathogen screening only (skip host/mtDNA alignment, qualimap, damage, and host summaries)
PATHOGEN_SCREENING_ONLY = CFG.get("pathogen_screening_only", False)
HOST_MTDNA_ANALYSIS_ENABLED = not PATHOGEN_SCREENING_ONLY

# In screening-only mode we rely on Bowtie2-based "host unaligned FASTQ" generation
# (via rules like `bowtie2_unaligned` -> `sample_unaligned_fastq_by_type`).
# Force pathogen_mapping_mode to `default` so we don't require `host_mapping/*host_unaligned.fastq.gz`.
if PATHOGEN_SCREENING_ONLY:
    PATHOGEN_MODE = "default"

# Host / mtDNA aligner: "bwa" or "bowtie2"
HOST_ALIGNER = CFG.get("host_aligner", "bwa")

# Adapter sequences (Illumina TruSeq defaults; override in config.yaml)
# AdapterRemoval (PE): adapter_removal_adapter1 / adapter_removal_adapter2
# cutadapt (SE):       cutadapt_adapter
# Legacy keys adapter1/adapter2 are used as fallbacks.
_ILLUMINA_ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
_ILLUMINA_ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
AR_ADAPTER1 = CFG.get("adapter_removal_adapter1", CFG.get("adapter1", _ILLUMINA_ADAPTER1))
AR_ADAPTER2 = CFG.get("adapter_removal_adapter2", CFG.get("adapter2", _ILLUMINA_ADAPTER2))
CUTADAPT_ADAPTER = CFG.get("cutadapt_adapter", CFG.get("adapter1", _ILLUMINA_ADAPTER1))

# FastQ Screen configuration file is generated by the workflow from config indices.
FASTQ_SCREEN_CONF = "results/workflow/fastq_screen.conf"

# decOM tool output (separate from gather_sinks p_sink / p_keys inputs)
DECOM_OUTPUT_DIR = "results/metagenomics/decOM/decOM_out"

# Aligner-neutral result dirs (tool = config host_aligner / pathogen_aligner: bwa | bowtie2)
HOST_MAPPING_DIR = "host_mapping"
MTDNA_MAPPING_DIR = "mtdna_mapping"
PATHOGEN_MAPPING_COMPLETE = "results/workflow/pathogen_mapping_complete.txt"

# Gate file written before the DAG is built (see ensure_pigsti_validated below).
PIGSTI_VALIDATION_OK = "results/workflow/.pigsti_validation_ok"


def ensure_pigsti_validated() -> None:
    """
    Fail fast on config/samples/paths when Snakemake loads this Snakefile.
    Not a separate workflow job — validation runs once before any rule is scheduled.
    """
    import subprocess
    import sys

    config_path = "config/config.yaml"
    spreadsheet = "config/Pathogen_spreadsheet.csv"
    inputs = [config_path, SAMPLES_TSV, spreadsheet]
    if os.path.isfile(PIGSTI_VALIDATION_OK):
        try:
            ok_mtime = os.path.getmtime(PIGSTI_VALIDATION_OK)
            if all(
                os.path.getmtime(p) <= ok_mtime
                for p in inputs
                if os.path.isfile(p)
            ):
                return
        except OSError:
            pass
    os.makedirs(os.path.dirname(PIGSTI_VALIDATION_OK) or ".", exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            "scripts/validate_pigsti_setup.py",
            "--config",
            config_path,
            "--samples",
            SAMPLES_TSV,
            "--spreadsheet",
            spreadsheet,
            "--output",
            PIGSTI_VALIDATION_OK,
        ],
        check=True,
    )


ensure_pigsti_validated()

# Cleanup intermediate files (set to false to keep all intermediate files)
CLEANUP_ENABLED = CFG.get("cleanup_intermediates", True)

# Helper wrapper so we can conditionally mark outputs as temp() based on config.
def maybe_temp(path):
    return temp(path) if CLEANUP_ENABLED else path

# Read-group mapping per PCR (SM is still the biological sample)
READ_GROUPS = {}
with open(SAMPLES_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    # Normalize header names here as well
    if reader.fieldnames is not None:
        reader.fieldnames = [fn.strip().lstrip("\ufeff") for fn in reader.fieldnames]
    has_pcr_col = "pcr" in reader.fieldnames if reader.fieldnames is not None else False
    for row in reader:
        # Skip empty / malformed rows
        if not row or "sample" not in row or row["sample"] is None or not str(row["sample"]).strip():
            continue
        sample = str(row["sample"]).strip()
        pcr_id = row["pcr"] if has_pcr_col and row.get("pcr") else sample
        rglb = str(row.get("RGLB", "")).strip()
        sequencing_run = str(row.get("sequencing_run", "")).strip()
        rg_id = f"{rglb}_{sequencing_run}"
        read_group = f"@RG\\tID:{rg_id}\\tPL:ILLUMINA\\tLB:{rglb}\\tSM:{sample}"
        READ_GROUPS[pcr_id] = read_group


def get_read_group(wc):
    """
    Return the appropriate @RG line for a given wildcard context.
    The wildcard 'sample' may be a PCR ID or a biological sample ID;
    in the latter case we map to the first PCR for that sample.
    """
    sid = wc.sample
    if sid not in READ_GROUPS and sid in SAMPLE_TO_PCRS:
        sid = SAMPLE_TO_PCRS[sid][0]
    return READ_GROUPS[sid]


def pcrs_for_sample(sample):
    """Return list of PCR IDs for a given biological sample."""
    return SAMPLE_TO_PCRS.get(sample, [])


def kraken_input_fastq_path(bio: str) -> str:
    """FASTQ path for KrakenUniq / HOPS / decOM (merged unaligned per biological sample)."""
    return f"results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz"


def kraken_input_fastq(wildcards):
    return kraken_input_fastq_path(wildcards.sample)


# Load spreadsheet once globally and strip columns once
spreadsheet_df = pd.read_csv("config/Pathogen_spreadsheet.csv")
spreadsheet_df.columns = spreadsheet_df.columns.str.strip()

# Load pathogen detection criteria (for E-value / reads / E-value thresholds)
# Precedence:
#  1) `pathogen_detection_criteria` section in config/config.yaml
#  2) `config/pathogen_detection_criteria.yaml` (backwards compatible fallback)
#  3) built-in defaults
_criteria_defaults = {
    "escore_threshold": 5,
    "reads_threshold": 50,
    "guellil_evalue_threshold": 0.001,  # Guellil E-value metric threshold (larger is better)
    "use_evalue_for_detection": True,
}
DETECTION_CFG = dict(_criteria_defaults)

_criteria_from_main = CFG.get("pathogen_detection_criteria")
if isinstance(_criteria_from_main, dict):
    DETECTION_CFG.update(_criteria_from_main)
else:
    _criteria_path = "config/pathogen_detection_criteria.yaml"
    if os.path.exists(_criteria_path):
        with open(_criteria_path) as _f:
            _legacy = yaml.safe_load(_f) or {}
        if isinstance(_legacy, dict):
            DETECTION_CFG.update(_legacy)

def _scripts_dir_on_path():
    sd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")
    if sd not in sys.path:
        sys.path.insert(0, sd)


def _sample_from_evalue_pathogen_csv(path):
    path = str(path).replace("\\", "/")
    m = re.search(r"results/pathogen/([^/]+)/evalue/pathogen/", path)
    if not m:
        m = re.search(r"results/[^/]+/([^/]+)/evalue/pathogen/", path)
    if not m:
        raise ValueError(f"Cannot parse biological sample from E-value path: {path}")
    return m.group(1)


def get_sample_ref_pairs(evalue_paths=None, hops_path=None, strict_cohort=False):
    """
    (sample, pathogen) pairs for mapping.

    For the pathogen_targets checkpoint, pass Snakemake ``input.evalue_files`` and
    ``strict_cohort=True`` so every BIO_SAMPLES has a completed E-value CSV (never
    scan disk with os.path.exists during --rerun-incomplete).

    When ``evalue_paths`` is None, only existing on-disk CSVs are used (legacy /
    cleanup markers — not for scheduling mapping).
    """
    _scripts_dir_on_path()
    from build_pathogen_target_pairs import compute_pairs, load_detection_cfg

    detection_cfg = load_detection_cfg("config/config.yaml")
    if evalue_paths is None:
        evalue_paths = [
            f"results/pathogen/{s}/evalue/pathogen/{s}_pathogen.csv"
            for s in BIO_SAMPLES
            if os.path.isfile(f"results/pathogen/{s}/evalue/pathogen/{s}_pathogen.csv")
        ]
        if hops_path is None and ENABLE_HOPS:
            hp = "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"
            hops_path = hp if os.path.isfile(hp) else None
        strict_cohort = False
    else:
        evalue_paths = list(evalue_paths)

    return compute_pairs(
        evalue_paths,
        list(BIO_SAMPLES),
        spreadsheet_df,
        detection_cfg,
        hops_path=hops_path,
        strict_cohort=strict_cohort,
    )

def get_single_end_samples():
    """Get samples that have no r2 reads (single-end only)"""
    return [sample for sample in SAMPLES if not SAMPLES_DICT[sample]["r2"]]

def get_paired_end_samples():
    """Get samples that have r2 reads (paired-end)"""
    return [sample for sample in SAMPLES if SAMPLES_DICT[sample]["r2"]]




def get_all_possible_pathogens():
    """Get all possible pathogens from the spreadsheet for all samples"""
    pairs = []
    pathogen_names = sorted(spreadsheet_df["Krakenuniq name"].dropna().astype(str).unique())
    for sample in sorted(SAMPLES):
        for pathogen in pathogen_names:
            pairs.append((sample, pathogen))
    return pairs







def get_reference_path(wc):
    row = spreadsheet_df[spreadsheet_df["Krakenuniq name"] == wc.krakenuniq_name]
    if row.empty:
        raise ValueError(f"No reference path found for {wc.krakenuniq_name}")
    return row.iloc[0]["bwa index"]


# (Legacy get_bwa_targets removed — pathogen mapping uses checkpoint targets + pathogen_mapping/)



def check_bwa_index(fasta_path):
    """Check if BWA index exists for a FASTA file"""
    if not fasta_path or not os.path.exists(fasta_path):
        return False
    
    # BWA index files
    index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    base_path = fasta_path.replace('.fa', '').replace('.fasta', '').replace('.fna', '')
    
    for ext in index_extensions:
        if not os.path.exists(base_path + ext):
            return False
    return True

def check_bowtie2_index(fasta_path):
    """Check if Bowtie2 index exists for a FASTA file"""
    if not fasta_path or not os.path.exists(fasta_path):
        return False
    
    # Bowtie2 index files
    index_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    base_path = fasta_path.replace('.fa', '').replace('.fasta', '').replace('.fna', '')
    
    for ext in index_extensions:
        if not os.path.exists(base_path + ext):
            return False
    return True

def get_pathogen_references():
    """Get all pathogen reference paths from spreadsheet"""
    spreadsheet_df = pd.read_csv(config['pathogen_spreadsheet'])
    references = []
    for _, row in spreadsheet_df.iterrows():
        if pd.notna(row['bwa index']) and row['bwa index'].strip():
            references.append(row['bwa index'].strip())
    return list(set(references))  # Remove duplicates

def get_host_references():
    """Get all host reference paths from config"""
    references = []
    host_map = config.get('bwa_indices') or {}
    if not isinstance(host_map, dict):
        host_map = {}
    for host, path in host_map.items():
        if path and path.strip():
            references.append(path.strip())
    return references

def get_mtdna_references():
    """Get all mtDNA reference paths from config"""
    references = []
    mtdna_map = config.get('mtDNA_indices') or {}
    if not isinstance(mtdna_map, dict):
        mtdna_map = {}
    for host, path in mtdna_map.items():
        if path and path.strip():
            references.append(path.strip())
    return references

def safe_name(name):
    """Filesystem-safe token for pathogen names in output paths (stable across runs)."""
    s = str(name).strip()
    for ch in (" ", "/", "\\", ":"):
        s = s.replace(ch, "_")
    while "__" in s:
        s = s.replace("__", "_")
    return s


def _parse_pathogen_target_line(line: str):
    """Return (bio_sample, pathogen_safe) from a pathogen_targets.txt BAM path line."""
    line = str(line).strip().replace("\\", "/")
    if not line or not line.endswith(".dedup.bam"):
        return None
    parts = line.split("/")
    # v3: results/pathogen/{bio}/pathogen_mapping/{bio}_{pathogen_safe}.dedup.bam
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


def _parse_pathogen_targets_file(targets_file: str):
    """
    Parse results/workflow/pathogen_targets.txt into (sample, pathogen_safe) pairs.
    Each line is expected to be:
      results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.dedup.bam
    """
    pairs = []
    with open(targets_file, "r") as f:
        for line in f:
            parsed = _parse_pathogen_target_line(line)
            if parsed:
                pairs.append(parsed)

    # Deduplicate and sort for stable DAG + reproducible ordering.
    pairs = sorted(set(pairs), key=lambda x: (x[0], x[1]))
    return pairs


def get_sample_ref_pairs_safe():
    """Pairs for cleanup markers: prefer checkpoint targets file, never disk-scan at parse time."""
    targets_file = "results/workflow/pathogen_targets.txt"
    if os.path.isfile(targets_file):
        return _parse_pathogen_targets_file(targets_file)
    return []


def _reference_fasta_is_gzip(path):
    """True if FASTA is gzip-compressed (.gz suffix or gzip magic), not bgzip."""
    p = str(path).strip()
    if p.lower().endswith(".gz"):
        return True
    try:
        with open(p, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _materialized_ref_plain_path(ref_name_safe):
    clean = ref_name_safe.split("/")[0]
    return f"results/workflow/index_status/ref_plain/{clean}.fa"


def get_spreadsheet_ref_path_for_safe(ref_name_safe):
    """Path from Pathogen_spreadsheet 'bwa index' column (may be gzip)."""
    clean_ref_name = ref_name_safe.split("/")[0]
    pathogen_name = clean_ref_name.replace("_", " ")
    row = spreadsheet_df.loc[spreadsheet_df["Krakenuniq name"] == pathogen_name]
    if row.empty:
        raise ValueError(f"[ERROR] No bwa index found for pathogen name: '{pathogen_name}'")
    val = row.iloc[0]["bwa index"]
    if pd.isna(val) or str(val).strip() == "":
        raise ValueError(f"Empty bwa index in spreadsheet for: {pathogen_name}")
    return str(val).strip()


def get_reference_from_safe_name(ref_name_safe):
    """
    FASTA path for aligners / DamageProfiler / pysam.
    If the spreadsheet points to gzip-compressed FASTA, use materialized plain copy
    (written by rule fasta_index) so samtools faidx and tools accept it.
    """
    raw = get_spreadsheet_ref_path_for_safe(ref_name_safe)
    if _reference_fasta_is_gzip(raw):
        return _materialized_ref_plain_path(ref_name_safe)
    return raw

def get_zipped_pairs():
    """Get zipped_pairs dynamically"""
    pairs = get_sample_ref_pairs()
    return [(s, safe_name(r)) for s, r in pairs]


def get_pathogen_reads(wc):
    """
    Select input reads for pathogen BWA mapping based on PATHOGEN_MODE.
    - default:      use Bowtie2-unmapped per-sample merged reads (faster)
    - super_careful: use BWA host-unmapped per-PCR reads (more stringent)
    """
    if PATHOGEN_MODE == "default":
        # Bowtie2-unmapped reads merged per biological sample.
        # wc.sample may be a PCR ID or a biological sample ID.
        # If it's already a bio sample, use it directly; otherwise map PCR -> bio sample.
        sid = wc.sample
        if sid in BIO_SAMPLES:
            bio = sid
        else:
            # It's a PCR ID, map to bio sample
            bio = PCR_INFO.get(sid, {}).get("sample", sid)
        return f"results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz"
    else:
        # BWA host-unmapped reads per PCR (current super-careful behaviour)
        return f"results/libraries/{wc.sample}/host_mapping/{wc.sample}_host_unaligned.fastq.gz"

def get_pathogen_reports(checkpoint_output):
    """Get pathogen report PDFs after checkpoint completion"""
    try:
        # Read from checkpoint output to get list of mapped pathogens
        if os.path.exists(checkpoint_output.targets_file):
            reports = []
            with open(checkpoint_output.targets_file) as f:
                for line in f:
                    line = line.strip()
                    if line.endswith(".dedup.bam"):
                        parsed = _parse_pathogen_target_line(line)
                        if parsed:
                            sample, pathogen_safe = parsed
                            reports.append(
                                f"results/pathogen/{sample}/summary/{sample}_{pathogen_safe}_pathogen_report.pdf"
                            )
            return reports
        return []
    except:
        return []

def get_cleanup_markers():
    """Get cleanup marker files if cleanup is enabled"""
    if not CLEANUP_ENABLED:
        return []
    markers = []
    # Prinseq cleanup
    markers.extend([f"results/libraries/{pcr}/prinseq/.cleanup_done" for pcr in SAMPLES])
    # Host cleanup (only when host/mtDNA analysis is enabled)
    if HOST_MTDNA_ANALYSIS_ENABLED:
        markers.extend([f"results/libraries/{pcr}/host_mapping/.sai_cleanup_done" for pcr in SAMPLES])
        markers.extend([f"results/libraries/{pcr}/host_mapping/.intermediate_cleanup_done" for pcr in SAMPLES])
        # mtDNA cleanup
        markers.extend([f"results/libraries/{pcr}/mtdna_mapping/.sai_cleanup_done" for pcr in SAMPLES])
        markers.extend([f"results/libraries/{pcr}/mtdna_mapping/.intermediate_cleanup_done" for pcr in SAMPLES])
    # Pathogen cleanup
    for sample, ref_name_safe in get_sample_ref_pairs_safe():
        markers.append(f"results/pathogen/{sample}/pathogen_mapping/.{ref_name_safe}_sai_cleanup_done")
        markers.append(f"results/pathogen/{sample}/pathogen_mapping/.{ref_name_safe}_intermediate_cleanup_done")
    return markers


###--------------------------wrappers-----------------------------------------------------


###-----------------------------------------------rules-------------------------------------------------------------------------------

rule all:
    input:
        PIGSTI_VALIDATION_OK,
        expand("results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz", sample=SAMPLES),
        expand("results/libraries/{sample}/prinseq/{sample}-passed.fq.gz", sample=SAMPLES),
        expand("results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt", sample=BIO_SAMPLES),
        expand("results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.html", sample=SAMPLES),
        expand("results/libraries/{sample}/fastq_screen/{sample}_best_species.txt", sample=SAMPLES),
        # Disabled for now: by-type fastq_screen symlinks can race with cleanup on some filesystems
        # expand(f"{BY_TYPE_ROOT}/fastq_screen/{{sample}}/{{sample}}.collapsed_screen.html", sample=SAMPLES),

        # Host/mtDNA analysis outputs (skip in pathogen_screening_only mode)
        *(
            [
                # Host alignments (final outputs only; intermediates are cleaned up if cleanup_intermediates is enabled)
                *expand("results/libraries/{sample}/host_mapping/{sample}.dedup.bam", sample=SAMPLES),
                *expand("results/libraries/{sample}/host_mapping/{sample}.dedup.bam.bai", sample=SAMPLES),
                *expand("results/libraries/{sample}/host_mapping/{sample}.dedup_q30_softclipped.cram", sample=SAMPLES),
                *expand("results/libraries/{sample}/host_mapping/{sample}.dedup.metrics.txt", sample=SAMPLES),
                *expand("results/libraries/{sample}/host_mapping/{sample}.q30_metrics.txt", sample=SAMPLES),
                *expand(f"{BY_TYPE_ROOT}/host_mapping/{{sample}}/{{sample}}.dedup.bam", sample=SAMPLES),

                # Per-sample merged host BAMs
                *expand("results/host/host_mapping/{sample}.dedup.merged.bam", sample=BIO_SAMPLES),

                # Host QC and damage
                *expand("results/libraries/{sample}/qualimap/genome_results.txt", sample=SAMPLES),
                *expand("results/libraries/{sample}/qualimap_mtdna/genome_results.txt", sample=SAMPLES),
                *expand("results/libraries/{sample}/damageprofiler_host", sample=SAMPLES),

                # Genetic sexing (Cow, Goat, Sheep, Dog only; others skipped in rule)
                *(
                    expand(
                        [
                            "results/samples/{sample}/sexing/{sample}_sexing.tsv",
                            "results/samples/{sample}/sexing/{sample}_sexing.pdf",
                        ],
                        sample=BIO_SAMPLES,
                    )
                    if ENABLE_SEXING
                    else []
                ),

                # Sample-level Qualimap on merged BAMs
                *expand("results/samples/{sample}/qualimap/genome_results.txt", sample=BIO_SAMPLES),
                *expand("results/samples/{sample}/qualimap_mtdna/genome_results.txt", sample=BIO_SAMPLES),

                # mtDNA alignments (final outputs only)
                *expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam", sample=SAMPLES),
                *expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam.bai", sample=SAMPLES),
                *expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup.metrics.txt", sample=SAMPLES),
                *expand("results/libraries/{sample}/mtdna_mapping/{sample}.q30_metrics.txt", sample=SAMPLES),
                *expand(f"{BY_TYPE_ROOT}/mtdna_mapping/{{sample}}/{{sample}}.dedup.bam", sample=SAMPLES),

                # Per-sample merged mtDNA BAMs
                *expand("results/host/mtdna_mapping/{sample}.dedup.merged.bam", sample=BIO_SAMPLES),
                *expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup_q30_softclipped.cram", sample=SAMPLES),

                # mtDNA QC and damage
                *expand("results/libraries/{sample}/damageprofiler_mtdna", sample=SAMPLES),
            ]
            if HOST_MTDNA_ANALYSIS_ENABLED
            else []
        ),

        # E-Score outputs
        expand("results/pathogen/{sample}/evalue/genus/{sample}_genus.csv", sample=BIO_SAMPLES),
        expand("results/pathogen/{sample}/evalue/species/{sample}_species.csv", sample=BIO_SAMPLES),
        expand("results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES),
        expand(f"{BY_TYPE_ROOT}/evalue/pathogen/{{sample}}/{{sample}}_pathogen.csv", sample=BIO_SAMPLES),

        # Comparison outputs (only if HOPS is enabled)
        *(expand("results/pathogen/{sample}/comparison/{sample}_comparison.tsv", sample=BIO_SAMPLES) if ENABLE_HOPS else []),
        *(expand("results/pathogen/{sample}/comparison/{sample}_heatmap.html", sample=BIO_SAMPLES) if ENABLE_HOPS else []),

        # decOM (optional)
        *(
            [
                "results/metagenomics/decOM/p_sink.txt",
                *expand("results/metagenomics/decOM/p_keys/{sample}.fof", sample=BIO_SAMPLES),
                DECOM_OUTPUT_DIR,
            ]
            if ENABLE_DECOM
            else []
        ),

        # Shared/global files
        expand("results/metagenomics/krakenuniq/{sample}/{sample}_output.txt", sample=BIO_SAMPLES),
        *(["results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"] if ENABLE_HOPS else []),
        "lists/krakenuniq_pathogen_list.txt",
        "lists/hops_pathogen_list.txt",
        "config/config_hops_custom.txt",
        "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv",
        "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_normalized.csv",
        "results/metagenomics/kraken_abundance/heatmap_absolute.pdf",
        "results/metagenomics/kraken_abundance/heatmap_normalized.pdf",

        # Do NOT list results/workflow/pathogen_targets.txt here: it is the checkpoint
        # generate_pathogen_targets output. Snakemake#823 — having a checkpoint
        # file as a direct target in rule all can leave bowtie2_pathogen_align etc.
        # stuck ("no jobs ready"). Downstream rules already pull it in via inputs
        # or checkpoints.generate_pathogen_targets.get().

        # Final summary reports
        *(["results/final/host_mtdna_summary_all_samples.xlsx"] if HOST_MTDNA_ANALYSIS_ENABLED else []),
        "results/final/pathogen_summary_all_samples.xlsx",
        "results/final/pathogen_detection_scores_heatmap.png",
        "results/final/pathogen_detection_scores_heatmap.pdf",
        "results/final/comprehensive_summary_all_samples.xlsx",
        *(["results/final/multi_qc_dashboard.html"] if (ENABLE_MULTI_QC_DASHBOARD and HOST_MTDNA_ANALYSIS_ENABLED) else []),
        *(["results/final/pathogen_summary_all_samples_pcr.xlsx"] if PATHOGEN_MODE == "super_careful" else []),
        # Per-sample reports
        *[
            f"results/pathogen/{sample}/summary/{sample}_sample_report.pdf"
            for sample in BIO_SAMPLES
            if HOST_MTDNA_ANALYSIS_ENABLED
        ],
        # Global barrier: all pathogen align + QC (qualimap, damage, etc.) must finish
        # before summarize_pathogen_data (via pathogen_mapping_complete.txt). Restores correct
        # coverage/ANI/entropy in summaries; expand_downstream_targets matches checkpoint.
        "results/workflow/pathogen_mapping_complete.txt",
        # Pathogen-specific PDF reports (generated from checkpoint targets).
        "results/workflow/pathogen_reports_complete.txt",
        # Reproducibility manifest (file hashes + config snapshot).
        "results/run_manifest.json",
        *([
            "results/metagenomics/pathogen_detection/detection_scores_heatmap.pdf",
            "results/metagenomics/pathogen_detection/detection_scores_matrix.csv",
            "results/metagenomics/pathogen_detection/detailed_scores.csv",
        ] if ENABLE_HOPS else []),

        # Pipeline execution report and workflow diagram (always; independent of HOPS)
        "results/final/pipeline_execution_report.html",
        "results/final/pipeline_timing_data.csv",
        "results/final/pipeline_workflow_diagram.html",
        "results/final/pipeline_workflow_diagram.png",
        "results/final/pipeline_workflow_diagram.svg",

        # Index status
        "results/workflow/index_status/pathogen_indices_built.txt",
        "results/workflow/index_status/host_indices_built.txt",
        *(["results/workflow/index_status/mtdna_indices_built.txt"] if HOST_MTDNA_ANALYSIS_ENABLED else []),

        # Cleanup markers (only if cleanup_intermediates is enabled in config)
        *get_cleanup_markers(),

        # Browse catalog (grouped view; does not move canonical outputs)
        *(
            ["results/final/output_catalog.tsv"]
            if BUILD_RESULTS_CATALOG
            else []
        ),



# ---------------------------------------------------------------------------
# Cohort gate before pathogen mapping checkpoint
# Snakemake only runs this when EVERY listed input is from a *completed* upstream job.
# That blocks --rerun-incomplete from firing the checkpoint while some evalue jobs
# are still incomplete but leave older CSVs on disk (the old os.path.exists bug).
# ---------------------------------------------------------------------------
rule metagenomics_screening_cohort_ready:
    input:
        evalue_pathogen=expand(
            "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES
        ),
        evalue_genus=expand(
            "results/pathogen/{sample}/evalue/genus/{sample}_genus.csv", sample=BIO_SAMPLES
        ),
        evalue_species=expand(
            "results/pathogen/{sample}/evalue/species/{sample}_species.csv", sample=BIO_SAMPLES
        ),
        kraken_report=expand(
            "results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt", sample=BIO_SAMPLES
        ),
        *(
            ["results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"]
            if ENABLE_HOPS
            else []
        ),
    output:
        flag="results/workflow/metagenomics_screening_cohort_ready.json",
    run:
        seen = set()
        for p in input.evalue_pathogen:
            sample = _sample_from_evalue_pathogen_csv(p)
            seen.add(sample)
            if os.path.getsize(p) < 2:
                raise ValueError(f"E-value pathogen CSV empty or truncated: {p}")
            df = pd.read_csv(p)
            if "taxonomy" not in df.columns:
                raise ValueError(f"E-value pathogen CSV missing 'taxonomy' column: {p}")
        missing = sorted(set(BIO_SAMPLES) - seen)
        if missing:
            raise ValueError(
                "metagenomics_screening_cohort_ready: not all biological samples have "
                f"completed E-score outputs. Missing: {missing}"
            )
        os.makedirs(os.path.dirname(output.flag), exist_ok=True)
        with open(output.flag, "w", encoding="utf-8") as fh:
            json.dump(
                {
                    "status": "ready",
                    "cohort_samples": sorted(BIO_SAMPLES),
                    "n_samples": len(BIO_SAMPLES),
                    "hops_enabled": bool(ENABLE_HOPS),
                },
                fh,
                indent=2,
            )


# Checkpoint: pathogen BAM targets from cohort E-score (+ HOPS) inputs only
checkpoint generate_pathogen_targets:
    input:
        cohort_ready="results/workflow/metagenomics_screening_cohort_ready.json",
        evalue_files=expand(
            "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES
        ),
        *(["results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"] if ENABLE_HOPS else []),
    output:
        targets_file="results/workflow/pathogen_targets.txt",
        manifest="results/workflow/pathogen_targets.manifest.json",
    conda:
        "workflow/envs/python.yaml",
    params:
        bio_samples=" ".join(BIO_SAMPLES),
        hops_arg=(
            '--hops "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv"'
            if ENABLE_HOPS
            else ""
        ),
    shell:
        r"""
        set -euo pipefail
        python scripts/build_pathogen_target_pairs.py \
            --config config/config.yaml \
            --spreadsheet config/Pathogen_spreadsheet.csv \
            --bio-samples {params.bio_samples} \
            --evalue {input.evalue_files} \
            {params.hops_arg} \
            --strict-cohort \
            --output-targets {output.targets_file} \
            --output-manifest {output.manifest}
        """


# -------------------- checkpoint-driven pathogen helpers --------------------
# These helpers avoid relying on `os.path.exists(...)` during DAG parse time.
# Instead, they read the checkpoint-produced `results/workflow/pathogen_targets.txt`.
# (_parse_pathogen_target_line / _parse_pathogen_targets_file are defined earlier,
# before get_sample_ref_pairs_safe, so parse-time cleanup markers work on reruns.)

# Cache checkpoint parsing by targets file mtime (safe across checkpoint re-runs).
_pathogen_checkpoint_pairs_cache = {"mtime": None, "pairs": None}


def pathogen_pairs_from_checkpoint(wildcards=None):
    chk = checkpoints.generate_pathogen_targets.get()
    targets_file = chk.output.targets_file
    mtime = os.path.getmtime(targets_file)
    if (
        _pathogen_checkpoint_pairs_cache["mtime"] != mtime
        or _pathogen_checkpoint_pairs_cache["pairs"] is None
    ):
        _pathogen_checkpoint_pairs_cache["pairs"] = _parse_pathogen_targets_file(
            targets_file
        )
        _pathogen_checkpoint_pairs_cache["mtime"] = mtime
    return _pathogen_checkpoint_pairs_cache["pairs"]


def pathogen_safe_names_for_sample(sample: str, wildcards=None):
    pairs = pathogen_pairs_from_checkpoint()
    return sorted({ps for s, ps in pairs if s == sample})


def expand_downstream_targets(wildcards):
    """
    Files that must exist before pathogen_mapping_complete.txt is written.
    Pairs are taken from the checkpoint targets file (same set as mapping jobs).
    """
    chk = checkpoints.generate_pathogen_targets.get()
    pairs = _parse_pathogen_targets_file(chk.output.targets_file)
    targets = []
    for sample, pathogen_safe in pairs:
        targets.extend(
            [
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.dedup.bam",
                f"results/pathogen/{sample}/pathogen_mapping/qualimap_{pathogen_safe}",
                f"results/pathogen/{sample}/pathogen_mapping/damageprofiler_{pathogen_safe}",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.ani.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.depth.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.breadth.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.expected_breadth.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.breadth_ratio.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.entropy_plot.png",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.mean_entropy.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.entropy_100bp.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.entropy_1000bp.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.edit_distance_logr2.txt",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.edit_distance_logr2.png",
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.ani_distribution.txt",
                *(
                    [
                        f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.edit_distance_logr2_damaged.txt",
                        f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.edit_distance_logr2_default.txt",
                    ]
                    if EDIT_DISTANCE_DAMAGE_SPLIT
                    else []
                ),
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{pathogen_safe}.genus_ranking.txt",
            ]
        )
    return targets


def pathogen_report_pdfs_from_checkpoint(wildcards=None):
    """
    Return all expected pathogen report PDFs derived from checkpoint targets.
    This enables Snakemake to parallelize report generation across jobs.
    """
    pairs = pathogen_pairs_from_checkpoint()
    return [
        f"results/pathogen/{sample}/summary/{sample}_{pathogen_safe}_pathogen_report.pdf"
        for sample, pathogen_safe in pairs
    ]


rule pathogen_mapping_targets:
    input:
        manifest="results/workflow/pathogen_targets.manifest.json",
        checkpoint_output="results/workflow/pathogen_targets.txt",
        downstream_targets=expand_downstream_targets,
    output:
        touch("results/workflow/pathogen_mapping_complete.txt"),
    run:
        with open(input.manifest, encoding="utf-8") as fh:
            man = json.load(fh)
        expected = sorted(BIO_SAMPLES)
        got = sorted(man.get("cohort_samples") or [])
        if got != expected:
            raise ValueError(
                "pathogen_targets.manifest.json cohort_samples mismatch: "
                f"expected {expected}, got {got}. Re-run generate_pathogen_targets "
                "after all E-score jobs finish."
            )
        if man.get("n_mapping_pairs", 0) == 0:
            print(
                "[pathogen_mapping_targets] Warning: zero pathogen mapping pairs in manifest; "
                "summaries may be empty.",
                file=sys.stderr,
            )
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", encoding="utf-8") as fh:
            fh.write(
                "Pathogen alignment and downstream analysis completed for cohort "
                f"({len(expected)} samples, {man.get('n_mapping_pairs', 0)} pairs)\n"
            )


#--------------------list and configs files -------------------------------------------

rule generate_pathogen_lists:
    input:
        spreadsheet = "config/Pathogen_spreadsheet.csv"
    output:
        kraken_list = "lists/krakenuniq_pathogen_list.txt",
        hops_list = "lists/hops_pathogen_list.txt"
    conda:
        "workflow/envs/python.yaml"  # or your python env
    shell:
        """
        python scripts/generate_pathogen_lists.py {input.spreadsheet}
        """

rule create_hops_config:
    input:
        original_config = "config/config_hops2.0.txt",
        hops_list = "lists/hops_pathogen_list.txt",
        cfg = "config/config.yaml"
    output:
        new_config = "config/config_hops_custom.txt"
    conda:
        "workflow/envs/python.yaml"
    shell:
        """
        python scripts/create_hops_config.py {input.original_config} {output.new_config} {input.hops_list} {input.cfg}
        """


#-------------------bwa pathogen mapping----------------------------------------------

if PATHOGEN_ALIGNER == "bwa":

    rule bwa_aln:
        input:
            reads = get_pathogen_reads,
            pathogen ="results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
            # Wait for checkpoint to complete (which waits for HOPS when enabled)
            checkpoint_targets = "results/workflow/pathogen_targets.txt",
            fai_ready = "results/workflow/index_status/faidx_{ref_name_safe}.done",
            bwa_ready = "results/workflow/index_status/bwa_{ref_name_safe}.done",
        output:
            sai = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.sai"
        conda:
            "workflow/envs/bwa.yaml"
        threads: 6
        params:
            reference = lambda wc: get_reference_from_safe_name(wc.ref_name_safe)
        shell:
            "bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {params.reference} {input.reads} > {output.sai}"

    rule bwa_samse:
        input:
            sai = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.sai",
            reads = get_pathogen_reads,
            bwa_ready = "results/workflow/index_status/bwa_{ref_name_safe}.done",
        output:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4.bam"
        conda:
            "workflow/envs/bwa.yaml"
        params:
            reference =  lambda wc: get_reference_from_safe_name(wc.ref_name_safe),
            read_group = get_read_group
        shell:
            """
            bwa samse -r '{params.read_group}' {params.reference} {input.sai} {input.reads} | \
            samtools view -F 4 -Sb - > {output.bam}
            """

    rule bwa_sort:
        input: "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4.bam"
        output: "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4.sorted.bam"
        conda: "workflow/envs/bwa.yaml"
        shell: "samtools sort {input} -o {output}"

    rule bwa_q30:
        input: "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4.sorted.bam"
        output: "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4_q30.bam"
        conda: "workflow/envs/bwa.yaml"
        shell: "samtools view -q 30 -o {output} {input}"

    rule bwa_q30_sort:
        input:
            "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4_q30.bam"
        output:
            "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4_q30_sort.bam"
        conda:
            "workflow/envs/samtools.yaml"
        shell:
            "samtools sort {input} -o {output}"

else:

    def get_bowtie2_pathogen_prefix(wc):
        """
        Derive Bowtie2 index prefix from the pathogen reference FASTA path.
        We assume the spreadsheet 'bwa index' column points to the FASTA file,
        and that Bowtie2 indices will be (or have been) built with the same basename.
        """
        ref = get_reference_from_safe_name(wc.ref_name_safe)
        return ref.replace(".fa", "").replace(".fasta", "").replace(".fna", "")

    rule bowtie2_pathogen_align:
        input:
            reads = get_pathogen_reads,
            # Wait for checkpoint to complete (which waits for HOPS when enabled)
            checkpoint_targets = "results/workflow/pathogen_targets.txt",
            fai_ready = "results/workflow/index_status/faidx_{ref_name_safe}.done",
        output:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4_q30_sort.bam"
        threads: 6
        conda:
            "workflow/envs/bowtie2.yaml"
        params:
            prefix = get_bowtie2_pathogen_prefix,
            ref = lambda wc: get_reference_from_safe_name(wc.ref_name_safe)
        shell:
            """
            mkdir -p results/samples/{wildcards.sample}/pathogen_mapping

            # ---------------- Concurrency-safe Bowtie2 index build ----------------
            idx_prefix="{params.prefix}"
            ref_fa="{params.ref}"
            lock_file="${{idx_prefix}}.bt2_build_lock"

            index_complete() {{
                # Consider the index complete if all .bt2 files are present
                for suf in 1 2 3 4 rev.1 rev.2; do
                    if [ ! -f "${{idx_prefix}}.${{suf}}.bt2" ]; then
                        return 1
                    fi
                done
                return 0
            }}

            if ! index_complete; then
                echo "[bowtie2_pathogen_align] Bowtie2 index not complete for {wildcards.ref_name_safe} at $idx_prefix"
                # Try to acquire build lock
                if ( set -o noclobber; > "$lock_file" ) 2>/dev/null; then
                    trap 'rm -f "$lock_file"' EXIT
                    echo "[bowtie2_pathogen_align] Building Bowtie2 index from $ref_fa ..."
                    bowtie2-build "$ref_fa" "$idx_prefix"
                    rm -f "$lock_file"
                    trap - EXIT
                else
                    echo "[bowtie2_pathogen_align] Another job is building the index, waiting..."
                    # Wait until index is complete
                    while ! index_complete; do
                        sleep 10
                    done
                    echo "[bowtie2_pathogen_align] Index for {wildcards.ref_name_safe} is now complete."
                fi
            fi

            # ---------------- Pathogen alignment ----------------
            bowtie2 --end-to-end --sensitive -x {params.prefix} -U {input.reads} -p {threads} | \
            samtools view -b -q 30 - | \
            samtools sort -@ {threads} -o {output.bam}
            """

rule mark_duplicates:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}_F4_q30_sort.bam"
    output:
        dedup = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
        metrics = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.metrics.txt"
    threads: MARKDUP_THREADS
    conda: "workflow/envs/picard.yaml"
    params:
        dedup_tool=DEDUP_TOOL
    shell:
        """
        if [ "{params.dedup_tool}" = "samtools" ]; then
            # Use samtools markdup with removal (requires coordinate-sorted input)
            tmp_sort="results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{wildcards.ref_name_safe}.markdup_tmp.bam"
            samtools sort -@ {threads} -o "$tmp_sort" {input.bam}
            samtools index "$tmp_sort"
            samtools markdup -r -@ {threads} "$tmp_sort" {output.dedup}
            # samtools markdup does not produce Picard-style metrics; write a simple stub
            echo -e "tool\\treads\\tnotes\nsamtools_markdup\\tNA\\tNo Picard-style metrics available" > {output.metrics}
            rm -f "$tmp_sort" "$tmp_sort.bai"
        else
            {PICARD_CMD} MarkDuplicates \
                I={input.bam} \
                O={output.dedup} \
                M={output.metrics} \
                REMOVE_DUPLICATES=true \
                ASSUME_SORTED=true \
                VALIDATION_STRINGENCY=SILENT
        fi
        """

rule index_bam:
    input:
        "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
    output:
        "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input}"


rule qualimap_bamqc_bwa:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
    output:
        report = directory("results/pathogen/{sample}/pathogen_mapping/qualimap_{ref_name_safe}")
    log:
        "logs/qualimap/{sample}_{ref_name_safe}.log"
    conda: "workflow/envs/qualimap.yaml"
    shell:
        """
        mkdir -p {output.report}
        mapped_count=$(samtools view -c -F 4 {input.bam} 2>/dev/null || echo "0")
        if [ "$mapped_count" -eq 0 ]; then
            echo "No mapped reads in {input.bam}. Skipping Qualimap." > {log} 2>&1
            echo "No mapped reads. Qualimap skipped." > {output.report}/NO_MAPPED_READS.txt
            touch {output.report}/genome_results.txt
        else
            qualimap bamqc \
                -bam {input.bam} \
                -outdir {output.report} \
                -outformat html > {log} 2>&1 || true
        fi
        """
rule damageprofiler:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
        ref = lambda wc: get_reference_from_safe_name(wc.ref_name_safe),
        fai_ready = "results/workflow/index_status/faidx_{ref_name_safe}.done"
    output:
        directory("results/pathogen/{sample}/pathogen_mapping/damageprofiler_{ref_name_safe}")
    conda: "workflow/envs/damageprofiler.yaml"
    shell:
        """
        mkdir -p {output}
        mapped_count=$(samtools view -c -F 4 {input.bam})
        if [ "$mapped_count" -eq 0 ]; then
            echo "No mapped reads in {input.bam}. Skipping DamageProfiler." >&2
            : > {output}/misincorporation.txt
            echo "No mapped reads. DamageProfiler skipped." > {output}/NO_MAPPED_READS.txt
        else
            damageprofiler -i {input.bam} -o {output} -r {input.ref}
        fi
        """

rule fasta_index:
    input:
        ref = lambda wc: get_spreadsheet_ref_path_for_safe(wc.ref_name_safe)
    output:
        done = "results/workflow/index_status/faidx_{ref_name_safe}.done"
    params:
        ref_clean = lambda wc: wc.ref_name_safe.split("/")[0]
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        src="{input.ref}"
        plain="results/workflow/index_status/ref_plain/{params.ref_clean}.fa"
        mkdir -p results/workflow/index_status/ref_plain

        is_gzip_file() {{
            [ -r "$1" ] && [ "$(od -An -N2 -tx1 "$1" 2>/dev/null | tr -d ' \n')" = "1f8b" ]
        }}

        if is_gzip_file "$src" || [[ "$src" == *.gz ]]; then
            idx_target="$plain"
            fai="$plain.fai"
        else
            idx_target="$src"
            fai="$src.fai"
        fi

        mkdir -p results/index_status

        if [ -s "$fai" ]; then
            echo "[fasta_index] Found existing FAI: $fai"
            touch {output.done}
            exit 0
        fi

        lock_file="${{fai}}.lock"

        index_complete() {{
            [ -s "$fai" ]
        }}

        if ! index_complete; then
            if ( set -o noclobber; echo $$ > "$lock_file" ) 2>/dev/null; then
                trap 'rm -f "$lock_file"' EXIT
                if [[ "$idx_target" == "$plain" ]]; then
                    echo "[fasta_index] gzip-compressed reference; decompressing to $plain"
                    gunzip -c "$src" > "$plain"
                fi
                echo "[fasta_index] Building FAI: $fai"
                samtools faidx "$idx_target"
                rm -f "$lock_file"
                trap - EXIT
            else
                echo "[fasta_index] Waiting for another job to finish FAI build: $fai"
                while ! index_complete; do
                    sleep 10
                done
            fi
        fi

        touch {output.done}
        """

rule bwa_index_pathogen:
    """Build BWA aln index next to the pathogen FASTA (faidx is not enough)."""
    input:
        fai_ready = "results/workflow/index_status/faidx_{ref_name_safe}.done"
    output:
        done = "results/workflow/index_status/bwa_{ref_name_safe}.done"
    params:
        reference = lambda wc: get_reference_from_safe_name(wc.ref_name_safe)
    log:
        "logs/index_building/bwa_{ref_name_safe}.log"
    conda:
        "workflow/envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail
        ref="{params.reference}"
        mkdir -p results/workflow/index_status logs/index_building
        echo "[bwa_index_pathogen] {wildcards.ref_name_safe}  $ref" > {log}

        if [ ! -r "$ref" ] || [ ! -s "$ref" ]; then
            echo "ERROR: pathogen FASTA missing or empty: $ref" | tee -a {log} >&2
            exit 1
        fi

        index_complete() {{
            for suf in .amb .ann .bwt .pac .sa; do
                if [ ! -s "${{ref}}${{suf}}" ]; then
                    return 1
                fi
            done
            return 0
        }}

        if index_complete; then
            echo "[bwa_index_pathogen] Existing BWA index for $ref" >> {log}
            touch {output.done}
            exit 0
        fi

        lock_file="${{ref}}.bwa_index.lock"
        if ( set -o noclobber; echo $$ > "$lock_file" ) 2>/dev/null; then
            trap 'rm -f "$lock_file"' EXIT
            echo "[bwa_index_pathogen] Running: bwa index $ref" >> {log}
            bwa index "$ref" >> {log} 2>&1
            rm -f "$lock_file"
            trap - EXIT
        else
            echo "[bwa_index_pathogen] Waiting for another job to finish: $ref" >> {log}
            while ! index_complete; do
                sleep 10
            done
        fi

        if ! index_complete; then
            echo "ERROR: BWA index still missing after build: $ref (.amb/.ann/.bwt/.pac/.sa)" | tee -a {log} >&2
            exit 1
        fi
        touch {output.done}
        """

# AdnaPlotter rule removed - too time and memory consuming
rule Compute_ANI:
    input:
        bam="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
    output:
        ani="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.ani.txt"
    log:
        "logs/ANI/{sample}_{ref_name_safe}.log"
    threads: 1
    conda:
        "workflow/envs/samtools.yaml"  # ensure samtools is included
    message:
        "Compute_ANI: Calculating ANI for {input.bam}"
    shell:
        """
        samtools stats {input.bam} 2>> {log} | \
        awk '/^SN/ && /mismatches:/ {{mis=$3}} /^SN/ && /bases mapped:/ {{map=$4}} \
        END {{if (map > 0) printf("ANI ≈ %.2f%%\\n", (1 - mis/map)*100); else print "No mapped bases!"}}' \
        > {output.ani}
        """

rule MappingStats:
    input:
        bam="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
        bai="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam.bai"
    output:
        depth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.depth.txt",
        breadth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.breadth.txt",
        expected_breadth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.expected_breadth.txt",
        ratio="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.breadth_ratio.txt"
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_breadth_stats.py"


rule EntropyProfile:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
        bai = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam.bai"
    output:
        plot = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.entropy_plot.png",
        mean = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.mean_entropy.txt",
        entropy_100bp = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.entropy_100bp.txt",
        entropy_1000bp = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.entropy_1000bp.txt",
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_entropy_profile.py"


rule EditDistanceR2:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
        bai = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam.bai"
    output:
        r2 = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2.txt",
        plot = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2.png",
        ani_distribution = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.ani_distribution.txt",
        edit_distance_distribution = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2_distribution.txt",
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_edit_distance_r2.py"

if EDIT_DISTANCE_DAMAGE_SPLIT:
    rule filter_pathogen_bam_damage_split:
        # Create damage vs no-damage BAM subsets for downstream edit-distance metrics.
        input:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam",
            ref_fasta = lambda wc: get_reference_from_safe_name(wc.ref_name_safe),
            # Ensure FASTA index exists for pysam.FastaFile access
            fai_ready = "results/workflow/index_status/faidx_{ref_name_safe}.done"
        output:
            damaged_bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.damage_subset.bam",
            default_bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.default_subset.bam"
        params:
            window = EDIT_DISTANCE_DAMAGE_WINDOW
        conda:
            "workflow/envs/qc.yaml"
        threads: 1
        script:
            "scripts/filter_bam_by_end_damage.py"

    rule EditDistanceR2_damaged:
        input:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.damage_subset.bam"
        output:
            r2 = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2_damaged.txt",
            plot = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2_damaged.png"
        conda:
            "workflow/envs/qc.yaml"
        params:
            metric_start_ed=1
        script:
            "scripts/calculate_edit_distance_r2.py"

    rule EditDistanceR2_default:
        input:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.default_subset.bam"
        output:
            r2 = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2_default.txt",
            plot = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.edit_distance_logr2_default.png"
        conda:
            "workflow/envs/qc.yaml"
        script:
            "scripts/calculate_edit_distance_r2.py"


rule GenusRanking:
    input:
        kraken_report = "results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt",
        spreadsheet = "config/Pathogen_spreadsheet.csv"
    output:
        ranking = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.genus_ranking.txt"
    params:
        pathogen_name = lambda wildcards: wildcards.ref_name_safe.replace("_", " ")
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_genus_ranking.py"



##plus---------------




# -------------------- Pathogen Tracker Preprocessing -----------------------

def get_adapter_removal_inputs_pe(wc):
    """Return inputs for paired-end adapter removal (only if r2 exists).

    The wildcard 'sample' may be either a PCR ID or a biological sample ID.
    If it's a biological sample, we map it to the first PCR for that sample.
    """
    sid = wc.sample
    if sid not in SAMPLES_DICT and sid in SAMPLE_TO_PCRS:
        # Map biological sample -> first PCR ID
        sid = SAMPLE_TO_PCRS[sid][0]
    if sid not in SAMPLES_DICT:
        raise Exception(f"Sample '{wc.sample}' not found in samples.tsv. Available PCR IDs: {list(SAMPLES_DICT.keys())}")

    if not SAMPLES_DICT[sid]["r2"]:
        # If no r2, this rule should not be applicable
        raise Exception(f"Sample/PCR {sid} has no r2 reads, paired-end rule not applicable")
    r1 = _resolve_existing_fastq_path(SAMPLES_DICT[sid]["r1"][0], sid, "R1")
    r2 = _resolve_existing_fastq_path(SAMPLES_DICT[sid]["r2"][0], sid, "R2")
    return {"r1": r1, "r2": r2}

def get_adapter_removal_inputs_se(wc):
    """Return inputs for single-end adapter removal (only if r2 is absent).

    The wildcard 'sample' may be either a PCR ID or a biological sample ID.
    If it's a biological sample, we map it to the first PCR for that sample.
    """
    sid = wc.sample
    if sid not in SAMPLES_DICT and sid in SAMPLE_TO_PCRS:
        sid = SAMPLE_TO_PCRS[sid][0]
    if sid not in SAMPLES_DICT:
        raise Exception(f"Sample '{wc.sample}' not found in samples.tsv. Available PCR IDs: {list(SAMPLES_DICT.keys())}")

    if SAMPLES_DICT[sid]["r2"]:
        # If r2 exists, this rule should not be applicable
        raise Exception(f"Sample/PCR {sid} has r2 reads, single-end rule not applicable")
    r1 = _resolve_existing_fastq_path(SAMPLES_DICT[sid]["r1"][0], sid, "R1")
    return {"r1": r1}

rule calculate_raw_reads:
    """
    Per-PCR raw read count, stored in a tiny text file.

    Stores raw read counts in a small QC file for summaries (raw FASTQs are never deleted by the pipeline).
    """
    input:
        r1=lambda wc: PCR_INFO[wc.sample]["r1"]
    output:
        "results/libraries/{sample}/qc/{sample}.raw_reads.txt"
    conda:
        "workflow/envs/qc.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/libraries/{wildcards.sample}/qc

        RAW="{input.r1}"
        if [ ! -f "$RAW" ]; then
            echo "0" > {output}
            exit 0
        fi

        if command -v pigz >/dev/null 2>&1; then
            DECOMPRESS="pigz -dc"
        else
            DECOMPRESS="zcat"
        fi

        # Count lines / 4 to get number of reads.
        # If decompression fails (e.g. corrupted FASTQ), fall back to 0 reads
        # instead of killing the whole pipeline.
        if ! $DECOMPRESS "$RAW" | wc -l | awk '{{print int($1/4)}}' > {output}; then
            echo "0" > {output}
        fi
        """


rule calculate_collapsed_reads:
    """
    Per-PCR collapsed read count from the AdapterRemoval output, stored in a tiny text file.

    This mirrors calculate_raw_reads but operates on the collapsed reads, so that we can
    safely delete very large .collapsed.gz files once downstream steps are done.
    """
    input:
        collapsed="results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz"
    output:
        "results/libraries/{sample}/qc/{sample}.collapsed_reads.txt"
    conda:
        "workflow/envs/qc.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/libraries/{wildcards.sample}/qc

        COLLAPSED="{input.collapsed}"
        if [ ! -f "$COLLAPSED" ]; then
            echo "0" > {output}
            exit 0
        fi

        if command -v pigz >/dev/null 2>&1; then
            DECOMPRESS="pigz -dc"
        else
            DECOMPRESS="zcat"
        fi

        # Count lines / 4 to get number of reads.
        # If decompression fails (e.g. corrupted FASTQ), fall back to 0 reads.
        if ! $DECOMPRESS "$COLLAPSED" | wc -l | awk '{{print int($1/4)}}' > {output}; then
            echo "0" > {output}
        fi
        """

rule adapter_removal_pe:
    input:
        r1 = lambda wc: get_adapter_removal_inputs_pe(wc)["r1"],
        r2 = lambda wc: get_adapter_removal_inputs_pe(wc)["r2"]
    output:
        collapsed = maybe_temp("results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz")
    log:
        "logs/adapter_removal/{sample}_pe.log"
    conda:
        "workflow/envs/adapterremoval.yaml"
    params:
        strict_inputs=lambda wc: "true" if STRICT_INPUTS else "false",
        adapter1=AR_ADAPTER1,
        adapter2=AR_ADAPTER2,
    threads: 6
    shell:
        """
        set -eo pipefail
        mkdir -p results/libraries/{wildcards.sample}/adapter_removal
        mkdir -p results/failures
        
        echo "Starting paired-end adapter removal for {wildcards.sample} at $(date)" > {log}
        echo "Input R1: {input.r1}" >> {log}
        echo "Input R2: {input.r2}" >> {log}

        if [ ! -r "{input.r1}" ] || [ ! -s "{input.r1}" ] || [ ! -r "{input.r2}" ] || [ ! -s "{input.r2}" ]; then
            echo "ERROR: Missing or empty input FASTQ(s) for {wildcards.sample}" >> {log}
            echo "R1 readable: $( [ -r "{input.r1}" ] && echo yes || echo no )" >> {log}
            echo "R1 exists+nonempty: $( [ -s "{input.r1}" ] && echo yes || echo no )" >> {log}
            echo "R2 readable: $( [ -r "{input.r2}" ] && echo yes || echo no )" >> {log}
            echo "R2 exists+nonempty: $( [ -s "{input.r2}" ] && echo yes || echo no )" >> {log}
            if [ "{params.strict_inputs}" = "true" ]; then
                exit 1
            fi
            echo -n | gzip > {output.collapsed}
            exit 0
        fi

        # Build adapter arguments (omit if empty/unset in config)
        ADAPTER_ARGS=""
        if [ -n "{params.adapter1}" ]; then
            ADAPTER_ARGS="$ADAPTER_ARGS --adapter1 {params.adapter1}"
        fi
        if [ -n "{params.adapter2}" ]; then
            ADAPTER_ARGS="$ADAPTER_ARGS --adapter2 {params.adapter2}"
        fi
        
        if AdapterRemoval --file1 "{input.r1}" --file2 "{input.r2}" \
        --basename "results/libraries/{wildcards.sample}/adapter_removal/{wildcards.sample}" \
        --threads {threads} --collapse --minadapteroverlap 1 \
        $ADAPTER_ARGS \
        --minlength 30 --gzip --trimns --trimqualities >> {log} 2>&1; then
        echo "Paired-end adapter removal completed for {wildcards.sample} at $(date)" >> {log}
        else
            EXIT_CODE=$?
            echo "ERROR: Adapter removal failed for {wildcards.sample} with exit code $EXIT_CODE" >> {log}
            if [ ! -f results/workflow/failures/Failing_samples.tsv ]; then
                echo -e "Sample\tReason\tTimestamp" > results/workflow/failures/Failing_samples.tsv
            fi
            echo -e "{wildcards.sample}\tAdapterRemoval_PE_failed\t$(date -Iseconds)" >> results/workflow/failures/Failing_samples.tsv
            if [ "{params.strict_inputs}" = "true" ]; then
                exit $EXIT_CODE
            fi
            echo "Creating empty output file to allow pipeline to continue..." >> {log}
            echo -n | gzip > {output.collapsed}
            exit 0
        fi
        """

rule adapter_removal_se:
    input:
        r1 = lambda wc: get_adapter_removal_inputs_se(wc)["r1"]
    output:
        collapsed = maybe_temp("results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz")
    log:
        "logs/adapter_removal/{sample}_se.log"
    conda:
        "workflow/envs/adapterremoval.yaml"
    params:
        strict_inputs=lambda wc: "true" if STRICT_INPUTS else "false",
        cutadapt_adapter=CUTADAPT_ADAPTER,
    threads: 6
    shell:
        """
        set -eo pipefail
        mkdir -p results/libraries/{wildcards.sample}/adapter_removal
        mkdir -p results/failures
        
        echo "Starting single-end adapter removal for {wildcards.sample} at $(date)" > {log}
        echo "Input R1: {input.r1}" >> {log}

        if [ ! -r "{input.r1}" ] || [ ! -s "{input.r1}" ]; then
            echo "ERROR: Missing or empty input FASTQ for {wildcards.sample}" >> {log}
            echo "R1 readable: $( [ -r "{input.r1}" ] && echo yes || echo no )" >> {log}
            echo "R1 exists+nonempty: $( [ -s "{input.r1}" ] && echo yes || echo no )" >> {log}
            if [ "{params.strict_inputs}" = "true" ]; then
                exit 1
            fi
            echo -n | gzip > {output.collapsed}
            exit 0
        fi

        # Build adapter argument (omit if empty/unset in config — cutadapt uses built-in defaults)
        ADAPTER_ARG=""
        if [ -n "{params.cutadapt_adapter}" ]; then
            ADAPTER_ARG="-a {params.cutadapt_adapter}"
        fi
        
        if cutadapt $ADAPTER_ARG \
                 -O 1 -m 30 \
                 "{input.r1}" -o "{output.collapsed}" -j {threads} >> {log} 2>&1; then
        echo "Single-end adapter removal completed for {wildcards.sample} at $(date)" >> {log}
        else
            EXIT_CODE=$?
            echo "ERROR: Adapter removal failed for {wildcards.sample} with exit code $EXIT_CODE" >> {log}
            if [ ! -f results/workflow/failures/Failing_samples.tsv ]; then
                echo -e "Sample\tReason\tTimestamp" > results/workflow/failures/Failing_samples.tsv
            fi
            echo -e "{wildcards.sample}\tAdapterRemoval_SE_failed\t$(date -Iseconds)" >> results/workflow/failures/Failing_samples.tsv
            if [ "{params.strict_inputs}" = "true" ]; then
                exit $EXIT_CODE
            fi
            echo "Creating empty output file to allow pipeline to continue..." >> {log}
            echo -n | gzip > {output.collapsed}
            exit 0
        fi
        """


rule generate_fastq_screen_conf:
    """
    Generate a FastQ Screen configuration file from the configured indices.
    This keeps the pipeline portable across systems without requiring a pre-shipped conf file.
    """
    output:
        conf=FASTQ_SCREEN_CONF
    # conda is not allowed with `run:` (Snakemake); this block only needs stdlib.
    run:
        import os

        # Prefer Bowtie2 index prefixes for FastQ Screen DATABASE lines.
        # If only BWA FASTA paths are available, derive a plausible prefix by stripping extensions.
        indices = {}
        if isinstance(CFG.get("bowtie2_indices"), dict) and CFG["bowtie2_indices"]:
            indices = dict(CFG["bowtie2_indices"])
        elif isinstance(CFG.get("bwa_indices"), dict) and CFG["bwa_indices"]:
            indices = {
                k: str(v).replace(".fa", "").replace(".fasta", "").replace(".fna", "")
                for k, v in CFG["bwa_indices"].items()
            }

        os.makedirs(os.path.dirname(output.conf), exist_ok=True)
        with open(output.conf, "w", encoding="utf-8") as f:
            f.write("# Auto-generated by PIGSTI (Snakemake)\n")
            f.write("# FastQ Screen configuration\n\n")
            f.write("# Use conda-provided aligner executables\n")
            f.write("BOWTIE\tbowtie\n")
            f.write("BOWTIE2\tbowtie2\n")
            f.write("#BWA\tbwa\n\n")
            f.write("THREADS\t\t8\n\n")
            f.write("##############\n")
            f.write("## DATABASES #\n")
            f.write("##############\n")

            for name in sorted(indices.keys()):
                path = str(indices[name]).strip()
                if not path:
                    continue
                # FastQ Screen supports an optional 3rd column (aligner); we pin to BOWTIE2.
                f.write(f"DATABASE\t{name}\t{path}\tBOWTIE2\n")


# FastQ Screen on AdapterRemoval collapsed reads (before PRINSEQ); host mapping uses the same collapsed FASTQ.
rule fastq_screen:
    input:
        collapsed="results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz",
        conf=FASTQ_SCREEN_CONF
    output:
        html="results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.html",
        txt="results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.txt",
        run_mode="results/libraries/{sample}/fastq_screen/{sample}_fastq_screen_run_mode.txt",
    params:
        full_rescreen=1 if FASTQ_SCREEN_FULL_RESCREEN else 0,
        min_one_hit=FASTQ_SCREEN_MIN_ONE_HIT,
    threads: 4
    conda:
        "workflow/envs/fastq_screen.yaml"
    shell:
        r"""
        set -eo pipefail
        mkdir -p results/libraries/{wildcards.sample}/fastq_screen

        COLLAPSED="{input.collapsed}"
        CONF="{input.conf}"
        OUTDIR="results/libraries/{wildcards.sample}/fastq_screen"
        OUTTXT="{output.txt}"
        OUTHTML="{output.html}"
        RUNMODE="{output.run_mode}"
        FASTQ_SCREEN_LOG="$OUTDIR/fastq_screen_error.log"

        write_stub_outputs() {{
            cat << EOF > "$OUTTXT"
# FastQ Screen report for {wildcards.sample}
# $1
EOF
            cat << EOF > "$OUTHTML"
<html><body><h2>FastQ Screen report for {wildcards.sample}</h2><p>$1</p></body></html>
EOF
        }}

        run_fastq_screen() {{
            local pass_label="$1"
            local use_full="$2"
            rm -f "$SUBSET_TEMP"
            # FastQ Screen skips if *_screen.png already exists; remove prior report files only.
            FQS_STEM=$(basename "$COLLAPSED" .gz)
            rm -f "$OUTTXT" "$OUTHTML"
            rm -f "$OUTDIR/${{FQS_STEM}}_screen.png" "$OUTDIR/${{FQS_STEM}}_screen.html" "$OUTDIR/${{FQS_STEM}}_screen.txt" 2>/dev/null || true
            rm -f "$OUTDIR"/*.png 2>/dev/null || true
            local cmd=("{FASTQ_SCREEN_EXE}" -conf "$CONF" "$COLLAPSED" --outdir "$OUTDIR")
            if [ "$use_full" = "1" ]; then
                cmd+=(--subset 0)
            fi
            if ! "${{cmd[@]}}" > "$FASTQ_SCREEN_LOG" 2>&1; then
                echo "[fastq_screen] ERROR: fastq_screen failed ($pass_label) for {wildcards.sample}" >&2
                cat "$FASTQ_SCREEN_LOG" >&2
            fi
            if [ ! -s "$OUTTXT" ] || [ ! -s "$OUTHTML" ]; then
                echo "[fastq_screen] WARNING: no output after $pass_label for {wildcards.sample}" >&2
                if [ -f "$FASTQ_SCREEN_LOG" ]; then
                    tail -20 "$FASTQ_SCREEN_LOG" >&2
                fi
                # Do not overwrite a valid subset report when full-dataset rescreen fails.
                if [ "$pass_label" = "subset_default" ] || [ ! -s "$OUTTXT" ]; then
                    write_stub_outputs "fastq_screen failed ($pass_label); see fastq_screen_error.log"
                fi
                return 1
            fi
            return 0
        }}

        # If collapsed FASTQ is missing or empty, create stub outputs and continue.
        if [ ! -s "$COLLAPSED" ]; then
            echo "[fastq_screen] No reads in $COLLAPSED, creating empty FastQ Screen report for {wildcards.sample}." >&2
            write_stub_outputs "No reads available (collapsed FASTQ empty or missing)."
            echo "pass=none reason=empty_input input=adapter_removal_collapsed" > "$RUNMODE"
            exit 0
        fi

        FILE_SIZE=$(stat -c%s "$COLLAPSED" 2>/dev/null || stat -f%z "$COLLAPSED" 2>/dev/null || echo "unknown")
        READ_COUNT=$(zcat "$COLLAPSED" 2>/dev/null | wc -l | awk '{{print int($1/4)}}' || echo "unknown")

        # FastQ Screen leaves a temp subset file next to the input; remove so each run is fresh.
        SUBSET_TEMP="$(dirname "$COLLAPSED")/$(basename "$COLLAPSED")_temp_subset.fastq"

        # Pass 1: default subset (~100k reads, FastQ Screen built-in default).
        run_fastq_screen "subset_default" 0 || true

        BEST_ONE_HIT=0
        if [ -s "$OUTTXT" ]; then
            BEST_ONE_HIT=$(python scripts/parse_fastq_screen.py \
                --screen-txt "$OUTTXT" \
                --sample {wildcards.sample} \
                --print-best-one-hit 2>/dev/null || echo 0)
        fi

        PASS_USED="subset_default"
        REASON="best_one_hit=${{BEST_ONE_HIT}}_ge_min_{params.min_one_hit}"

        # Pass 2: full dataset when enabled and best #One_hit_one_genome is below threshold.
        if [ {params.full_rescreen} -eq 1 ] && [ -s "$OUTTXT" ]; then
            if ! python scripts/parse_fastq_screen.py \
                --screen-txt "$OUTTXT" \
                --sample {wildcards.sample} \
                --min-one-hit {params.min_one_hit} \
                --check-rescreen >/dev/null 2>&1; then
                echo "[fastq_screen] Best species one-hit=${{BEST_ONE_HIT}} < {params.min_one_hit}; re-running with --subset 0 (full dataset) for {wildcards.sample}" >&2
                if run_fastq_screen "full_dataset" 1; then
                    PASS_USED="full_dataset"
                    REASON="rescreen_best_one_hit_${{BEST_ONE_HIT}}_lt_{params.min_one_hit}"
                    BEST_ONE_HIT=$(python scripts/parse_fastq_screen.py \
                        --screen-txt "$OUTTXT" \
                        --sample {wildcards.sample} \
                        --print-best-one-hit 2>/dev/null || echo 0)
                else
                    REASON="full_dataset_failed_after_subset_best_one_hit_${{BEST_ONE_HIT}}"
                fi
            fi
        fi

        echo "pass=$PASS_USED best_one_hit=$BEST_ONE_HIT min_one_hit={params.min_one_hit} full_rescreen={params.full_rescreen} input=adapter_removal_collapsed $REASON" > "$RUNMODE"
        echo "file_size=$FILE_SIZE read_count=$READ_COUNT" >> "$RUNMODE"
        """


# PRINSEQ++ dedup/complexity filter for the pathogen / metagenomics branch (after FastQ Screen).
rule prinseq:
    input:
        collapsed="results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz"
    output:
        passed=maybe_temp("results/libraries/{sample}/prinseq/{sample}-passed.fq.gz")
    threads: 6
    conda:
        "workflow/envs/prinseq.yaml"
    shell:
        """
        mkdir -p results/libraries/{wildcards.sample}/prinseq
        prinseq++ -fastq {input.collapsed} -derep 14 \
        -out_good results/libraries/{wildcards.sample}/prinseq/{wildcards.sample}-passed.fq \
        -out_bad results/libraries/{wildcards.sample}/prinseq/{wildcards.sample}-bad.fq \
        -VERBOSE 2 -threads {threads}
        pigz results/libraries/{wildcards.sample}/prinseq/{wildcards.sample}-passed.fq
        rm -f results/libraries/{wildcards.sample}/prinseq/{wildcards.sample}-bad.fq
        """


rule fastq_screen_by_type:
    input:
        html = "results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.html",
        txt  = "results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.txt"
    output:
        html = f"{BY_TYPE_ROOT}/fastq_screen/{{sample}}/{{sample}}.collapsed_screen.html",
        txt  = f"{BY_TYPE_ROOT}/fastq_screen/{{sample}}/{{sample}}.collapsed_screen.txt"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/fastq_screen/{wildcards.sample}
        # Remove existing symlinks if they exist (stale or broken)
        rm -f {output.html} {output.txt} 2>/dev/null || true
        # Create new symlinks
        ln -sf $(realpath {input.html}) {output.html}
        ln -sf $(realpath {input.txt}) {output.txt}
        """


rule prinseq_by_type:
    input:
        fastq = "results/libraries/{sample}/prinseq/{sample}-passed.fq.gz"
    output:
        fastq = f"{BY_TYPE_ROOT}/prinseq/{{sample}}/{{sample}}-passed.fq.gz"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/prinseq/{wildcards.sample}
        rm -f {output.fastq} 2>/dev/null || true
        ln -sf $(realpath {input.fastq}) {output.fastq}
        """

rule bowtie2_unaligned:
    input:
        passed="results/libraries/{sample}/prinseq/{sample}-passed.fq.gz"
    output:
        bam="results/libraries/{sample}/bowtie2/{sample}_unaligned.bam"
    threads: 6
    conda:
        "workflow/envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x {config[host_index]} -U {input.passed} -p {threads} | \
        samtools view -Sb - -f4 > {output.bam}
        """


def sample_pcrs(wc):
    """Helper for merge rule: PCRs for a given biological sample wildcard."""
    return pcrs_for_sample(wc.sample)


rule merge_unaligned_fastq_per_sample:
    input:
        lambda wc: expand(
            "results/libraries/{pcr}/unaligned_fastq/{pcr}_unaligned.fastq.gz",
            pcr=sample_pcrs(wc)
        )
    output:
        # Keep merged unaligned FASTQ permanently so it can be reused by KrakenUniq,
        # decOM, and standalone tools.
        merged="results/pools/unaligned_fastq/{sample}_unaligned.fastq.gz"
    threads: 2
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        """
        mkdir -p results/pools/unaligned_fastq
        cat {input} > {output.merged}
        """


rule sample_unaligned_fastq_by_type:
    input:
        fastq = "results/pools/unaligned_fastq/{sample}_unaligned.fastq.gz"
    output:
        fastq = f"{BY_TYPE_ROOT}/pools_unaligned_fastq/{{sample}}/{{sample}}_unaligned.fastq.gz"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/pools_unaligned_fastq/{wildcards.sample}
        rm -f {output.fastq} 2>/dev/null || true
        ln -sf $(realpath {input.fastq}) {output.fastq}
        """

rule bam_to_fastq:
    input:
        bam="results/libraries/{sample}/bowtie2/{sample}_unaligned.bam"
    output:
        fastq="results/libraries/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz"
    threads: 2
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools bam2fq {input.bam} -@{threads} | pigz - > {output.fastq}"


rule unaligned_fastq_by_type:
    input:
        fastq = "results/libraries/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz"
    output:
        fastq = f"{BY_TYPE_ROOT}/unaligned_fastq/{{sample}}/{{sample}}_unaligned.fastq.gz"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/unaligned_fastq/{wildcards.sample}
        rm -f {output.fastq} 2>/dev/null || true
        ln -sf $(realpath {input.fastq}) {output.fastq}
        """

rule krakenuniq:
    input:
        fastq=kraken_input_fastq,
    output:
        report="results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt",
        output="results/metagenomics/krakenuniq/{sample}/{sample}_output.txt"
    log:
        "logs/krakenuniq/{sample}.log"
    threads: 8
    conda:
        "workflow/envs/krakenuniq.yaml"
    shell:
        """
        echo "Starting KrakenUniq classification for {wildcards.sample} at $(date)" > {log}
        krakenuniq --version >> {log} 2>&1 || true
        echo "Input FASTQ: {input.fastq}" >> {log}
        echo "Database: {config[kraken_db]}" >> {log}
        
        krakenuniq --db {config[kraken_db]} --fastq-input {input.fastq} \
        --threads {threads} --output {output.output} --report-file {output.report} \
        --gzip-compressed --only-classified-out >> {log} 2>&1
        
        echo "KrakenUniq classification completed for {wildcards.sample} at $(date)" >> {log}
        """


rule krakenuniq_by_type:
    input:
        report = "results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt",
        output_file = "results/metagenomics/krakenuniq/{sample}/{sample}_output.txt"
    output:
        report = f"{BY_TYPE_ROOT}/krakenuniq/{{sample}}/{{sample}}_kraken-report.txt",
        output_file = f"{BY_TYPE_ROOT}/krakenuniq/{{sample}}/{{sample}}_output.txt"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/krakenuniq/{wildcards.sample}
        rm -f {output.report} {output.output_file} 2>/dev/null || true
        ln -sf $(realpath {input.report}) {output.report}
        ln -sf $(realpath {input.output_file}) {output.output_file}
        """

rule evalue:
    input:
        report="results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt",
        config="config/config.yaml"
    output:
        genus="results/pathogen/{sample}/evalue/genus/{sample}_genus.csv",
        species="results/pathogen/{sample}/evalue/species/{sample}_species.csv",
        pathogen="results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv"
    conda:
        "workflow/envs/evalue.yaml"
    params:
        script="scripts/dExp_Escore.py"
    shell:
        """
        python {params.script} {input.report} {output.genus} {output.species} {output.pathogen} {input.config}
        """


rule evalue_by_type:
    input:
        genus   = "results/pathogen/{sample}/evalue/genus/{sample}_genus.csv",
        species = "results/pathogen/{sample}/evalue/species/{sample}_species.csv",
        pathogen= "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv"
    output:
        genus   = f"{BY_TYPE_ROOT}/evalue/genus/{{sample}}/{{sample}}_genus.csv",
        species = f"{BY_TYPE_ROOT}/evalue/species/{{sample}}/{{sample}}_species.csv",
        pathogen= f"{BY_TYPE_ROOT}/evalue/pathogen/{{sample}}/{{sample}}_pathogen.csv"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/evalue/genus/{wildcards.sample}
        mkdir -p {BY_TYPE_ROOT}/evalue/species/{wildcards.sample}
        mkdir -p {BY_TYPE_ROOT}/evalue/pathogen/{wildcards.sample}
        rm -f {output.genus} {output.species} {output.pathogen} 2>/dev/null || true
        ln -sf $(realpath {input.genus}) {output.genus}
        ln -sf $(realpath {input.species}) {output.species}
        ln -sf $(realpath {input.pathogen}) {output.pathogen}
        """


rule parse_fastq_screen:
    input:
        "results/libraries/{sample}/fastq_screen/{sample}.collapsed_screen.txt"
    output:
        "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    params:
        exclude_human=True
    script:
        "scripts/parse_fastq_screen.py"
    

if ENABLE_HOPS:

    _hops_evalue_inputs = dict(
        genus=expand("results/pathogen/{sample}/evalue/genus/{sample}_genus.csv", sample=BIO_SAMPLES),
        species=expand("results/pathogen/{sample}/evalue/species/{sample}_species.csv", sample=BIO_SAMPLES),
        pathogen=expand("results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES),
        config="config/config_hops_custom.txt",
    )

    if HOPS_PARALLEL:
        if HOPS_MALT_MMAP:
            rule hops_malt_parallel_mmap:
                input:
                    fastqs=[kraken_input_fastq_path(bio) for bio in BIO_SAMPLES],
                    config="config/config_hops_custom.txt",
                    script="scripts/run_parallel_malt.py",
                    evalue_gate=expand(
                        "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
                        sample=BIO_SAMPLES,
                    ),
                output:
                    done=expand(f"{HOPS_MALT_ROOT}/{{bio}}/.malt_done", bio=BIO_SAMPLES),
                params:
                    malt_root=HOPS_MALT_ROOT,
                    bio_samples=" ".join(BIO_SAMPLES),
                    parallel_jobs=HOPS_PARALLEL_JOBS,
                    threads_per_job=HOPS_THREADS_PER_JOB,
                    heap_gb=HOPS_HEAP_GB,
                threads:
                    HOPS_PARALLEL_JOBS * HOPS_THREADS_PER_JOB,
                conda:
                    "workflow/envs/hops.yaml",
                log:
                    "logs/hops_malt/parallel.log",
                shell:
                    r"""
                    set -eo pipefail
                    python {input.script} \
                      --bio-samples {params.bio_samples} \
                      --fastq {input.fastqs} \
                      --malt-root {params.malt_root} \
                      --config {input.config} \
                      --parallel {params.parallel_jobs} \
                      --threads {params.threads_per_job} \
                      --heap-gb {params.heap_gb} \
                      --memory-mode map \
                      --resume \
                      --log {log}
                    """

        else:
            rule hops_malt_per_sample:
                input:
                    fq=lambda wildcards: kraken_input_fastq_path(wildcards.bio),
                    config="config/config_hops_custom.txt",
                    # Scheduling gate: wait for cohort Kraken/E-value, not used by hops itself.
                    evalue_gate=expand(
                        "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
                        sample=BIO_SAMPLES,
                    ),
                output:
                    done=f"{HOPS_MALT_ROOT}/{{bio}}/.malt_done",
                params:
                    malt_root=HOPS_MALT_ROOT,
                    heap_gb=HOPS_HEAP_GB,
                    rma=f"{HOPS_MALT_ROOT}/{{bio}}/malt/{{bio}}_unaligned.rma6",
                resources:
                    hops_jobs=1,
                threads:
                    HOPS_THREADS_PER_JOB,
                conda:
                    "workflow/envs/hops.yaml",
                log:
                    "logs/hops_malt/{bio}.log",
                shell:
                    r"""
                    set -eo pipefail
                    out="{params.malt_root}/{wildcards.bio}"
                    rm -rf "$out"
                    hops -Xmx{params.heap_gb}g \
                      -input {input.fq} \
                      -output "$out" \
                      -m malt \
                      -c {input.config} \
                      > {log} 2>&1
                    if [ ! -s "{params.rma}" ]; then
                        echo "[hops_malt] expected RMA missing: {params.rma}" >> {log}
                        echo "[hops_malt] files under $out:" >> {log}
                        find "$out" -type f 2>/dev/null >> {log} || true
                        tail -n 40 {log} >&2 || true
                        exit 1
                    fi
                    touch {output.done}
                    """

        rule hops_stage_rma:
            input:
                malt_done=[hops_malt_done_path(b) for b in BIO_SAMPLES],
                script="scripts/stage_hops_rma.py",
            output:
                [hops_staged_rma_path(b) for b in BIO_SAMPLES],
            params:
                bio_samples=" ".join(BIO_SAMPLES),
                malt_root=HOPS_MALT_ROOT,
                rma_root=HOPS_RMA_ROOT,
            conda:
                "workflow/envs/python.yaml",
            log:
                "logs/hops_stage_rma.log",
            shell:
                r"""
                python {input.script} \
                  --bio-samples {params.bio_samples} \
                  --malt-root {params.malt_root} \
                  --rma-root {params.rma_root} \
                  --output {output} \
                  > {log} 2>&1
                """

        rule hops_maltex_post:
            input:
                rma=[hops_staged_rma_path(b) for b in BIO_SAMPLES],
                config="config/config_hops_custom.txt",
            output:
                heatmap=HOPS_HEATMAP,
            params:
                heap_gb=HOPS_HEAP_GB,
                rma_root=HOPS_RMA_ROOT,
            threads:
                HOPS_MALTEX_THREADS,
            conda:
                "workflow/envs/hops.yaml",
            log:
                "logs/hops_maltex_post.log",
            shell:
                r"""
                set -eo pipefail
                hops -Xmx{params.heap_gb}g \
                  -input {params.rma_root}/ \
                  -output results/metagenomics/hops \
                  -m me_po \
                  -c {input.config} \
                  > {log} 2>&1
                test -s {output.heatmap}
                """

    else:
        rule hops:
            input:
                **_hops_evalue_inputs,
                fq=[kraken_input_fastq_path(bio) for bio in BIO_SAMPLES],
            output:
                heatmap=HOPS_HEATMAP,
            params:
                heap_gb=HOPS_HEAP_GB,
            threads:
                HOPS_THREADS_PER_JOB,
            conda:
                "workflow/envs/hops.yaml",
            log:
                "logs/hops_cohort.log",
            shell:
                r"""
                set -eo pipefail
                hops -Xmx{params.heap_gb}g \
                  -input {input.fq} \
                  -output results/metagenomics/hops \
                  -m full \
                  -c {input.config} \
                  > {log} 2>&1
                test -s {output.heatmap}
                """

    rule compare_pathogens:
        input:
            evalue="results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
            hops="results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv",
            spreadsheet="config/Pathogen_spreadsheet.csv"
        output:
            tsv="results/pathogen/{sample}/comparison/{sample}_comparison.tsv",
            html="results/pathogen/{sample}/comparison/{sample}_heatmap.html"
        params:
            sample="{sample}"
        conda:
            "workflow/envs/compare_pathogens.yaml"
        shell:
            "python scripts/compare_pathogens.py {input.evalue} {input.hops} {params.sample} {input.spreadsheet}"


##------------------------------------------Microbiome--------------------#######################

rule gather_sinks:
    input:
        fastqs=[kraken_input_fastq_path(bio) for bio in BIO_SAMPLES],
    output:
        sink_txt="results/metagenomics/decOM/p_sink.txt",
        fof_files=expand("results/metagenomics/decOM/p_keys/{sample}.fof", sample=BIO_SAMPLES)
    run:
        import os
        os.makedirs("results/metagenomics/decOM/p_keys", exist_ok=True)
        missing = []
        with open(output.sink_txt, "w") as f_sink:
            for sample in BIO_SAMPLES:
                f_sink.write(sample + "\n")
                fq_path = os.path.abspath(kraken_input_fastq_path(sample))
                if not os.path.isfile(fq_path):
                    missing.append(fq_path)
                with open(f"results/metagenomics/decOM/p_keys/{sample}.fof", "w") as f_fof:
                    f_fof.write(f"{sample}: {fq_path}\n")
        if missing:
            raise ValueError(
                "decOM FASTQ inputs missing (run merge_unaligned_fastq_per_sample first): "
                + ", ".join(missing)
            )


rule decom_run:
    input:
        p_sink="results/metagenomics/decOM/p_sink.txt",
        fof_files=expand("results/metagenomics/decOM/p_keys/{sample}.fof", sample=BIO_SAMPLES),
        wrapper="scripts/run_decom.py",
        # Wait for HOPS before decOM so parallel MALT jobs do not compete for RAM.
        **({"hops_done": HOPS_HEATMAP} if ENABLE_HOPS else {}),
    output:
        directory(DECOM_OUTPUT_DIR),
    params:
        p_sources=config["decOM_sources"],
        memory=DECOM_MEMORY,
        threads=DECOM_THREADS,
        fail_flag="--fail-on-error" if DECOM_FAIL_ON_ERROR else "",
    resources:
        decom=1,
    threads:
        DECOM_THREADS,
    conda:
        "workflow/envs/decom.yaml"
    log:
        "logs/decom_run.log"
    shell:
        r"""
        set -eo pipefail
        mkdir -p results/metagenomics/decOM/p_keys
        python {input.wrapper} --p-sink {input.p_sink} --p-sources {params.p_sources} --p-keys-dir results/metagenomics/decOM/p_keys --memory {params.memory} --threads {params.threads} --output {output} --log {log} {params.fail_flag}
        """
## Krona plots removed (user no longer needs krona HTML outputs).
rule krakenuniq_abundance_matrix:
    input:
        reports = expand("results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt", sample=BIO_SAMPLES),
        script1 = "scripts/krakenuniq_abundance_matrix.R",
        script2 = "scripts/plot_krakenuniq_abundance_matrix.R"
    output:
        matrix = "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv",
        matrix_norm = "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_normalized.csv",
        abs_plot = "results/metagenomics/kraken_abundance/heatmap_absolute.pdf",
        norm_plot = "results/metagenomics/kraken_abundance/heatmap_normalized.pdf"
    conda:
        "workflow/envs/r.yaml"
    shell:
        """
        set -euo pipefail
        outdir="results/metagenomics/kraken_abundance"
        mkdir -p "$outdir"
        Rscript {input.script1} results/metagenomics/krakenuniq "$outdir" 1000 25
        Rscript {input.script2} "$outdir" {output.abs_plot} {output.norm_plot}
        """

# -------------------- Host bwa aln Mapping -----------------------------
rule bwa_aln_host:
    input:
        reads = "results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz",
        species = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4.bam",
        unmapped_bam = "results/libraries/{sample}/host_mapping/{sample}_host_unaligned.bam",
        unmapped_fastq = "results/libraries/{sample}/host_mapping/{sample}_host_unaligned.fastq.gz"
    threads: 6
    conda:
        "workflow/envs/bwa.yaml"
    params:
        lb   = lambda wc: PCR_INFO[wc.sample]["RGLB"],
        rgid = lambda wc: PCR_INFO[wc.sample]["RGLB"]
    script:
        "scripts/bwa_aln_host.py"


rule sort_bam_initial_host_mapping:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4.bam"
    output:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"

rule filter_q30_host_mapping:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4.sorted.bam"
    output:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools view -q 30 -@ {threads} -b {input.bam} -o {output.bam}"

rule index_q30_bam_host:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam"
    output:
        bai = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule sort_q30_bam_host_mapping:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.bam"
    output:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"

rule mark_duplicates_host:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam"
    output:
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam",
        metrics = "results/libraries/{sample}/host_mapping/{sample}.dedup.metrics.txt"
    threads: MARKDUP_THREADS
    conda: "workflow/envs/picard.yaml"
    params:
        dedup_tool=DEDUP_TOOL
    shell:
        """
        if [ "{params.dedup_tool}" = "samtools" ]; then
            tmp_sort="results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}.markdup_tmp.bam"
            samtools sort -@ {threads} -o "$tmp_sort" {input.bam}
            samtools index "$tmp_sort"
            samtools markdup -r -@ {threads} "$tmp_sort" {output.bam}
            echo -e "tool\\treads\\tnotes\nsamtools_markdup\\tNA\\tNo Picard-style metrics available" > {output.metrics}
            rm -f "$tmp_sort" "$tmp_sort.bai"
        else
            {PICARD_CMD} MarkDuplicates \
                I={input.bam} \
                O={output.bam} \
                M={output.metrics} \
                REMOVE_DUPLICATES=true \
                ASSUME_SORTED=true \
                VALIDATION_STRINGENCY=SILENT
        fi
        """
rule index_dedup_bam_host:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam"
    output:
        bai = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule save_host_q30_metrics:
    input:
        q30_bam = "results/libraries/{sample}/host_mapping/{sample}_F4_q30.sorted.bam",
        dedup_bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam"
    output:
        metrics_file = "results/libraries/{sample}/host_mapping/{sample}.q30_metrics.txt"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        """
        # Count Q30 reads and dedup reads, then calculate duplication rate
        q30_count=$(samtools view -c -F 4 {input.q30_bam} 2>/dev/null || echo "0")
        dedup_count=$(samtools view -c -F 4 {input.dedup_bam} 2>/dev/null || echo "0")
        
        # Calculate duplication rate
        if [ "$q30_count" -gt 0 ]; then
            dup_rate=$(echo "scale=4; 1 - ($dedup_count / $q30_count)" | bc 2>/dev/null || echo "NA")
        else
            dup_rate="NA"
        fi
        
        # Save to file
        echo "q30_reads=$q30_count" > {output.metrics_file}
        echo "dedup_reads=$dedup_count" >> {output.metrics_file}
        echo "duplication_rate=$dup_rate" >> {output.metrics_file}
        """


rule damage_profiler_host:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam",
        bai = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam.bai",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt",
        config = "config/config.yaml"
    output:
        dir = directory("results/libraries/{sample}/damageprofiler_host/")
    log:
        "logs/damageprofiler/{sample}.log"
    conda:
        "workflow/envs/damageprofiler.yaml"
    params:
        # Use host reference map that matches selected host aligner.
        index_key=("bowtie2_indices" if HOST_ALIGNER == "bowtie2" else "bwa_indices")
    script:
        "scripts/run_damageprofiler.py"

def _mtdna_reference_map_from_cfg(cfg):
    """
    Return species -> mtDNA reference map.
    Prefer mtDNA_indices (FASTA paths). If absent, fallback to
    bowtie2_mtDNA_indices (prefixes/paths resolved later by pick_fasta).
    """
    mt_map = cfg.get("mtDNA_indices")
    if isinstance(mt_map, dict) and mt_map:
        return dict(mt_map)
    bt_mt_map = cfg.get("bowtie2_mtDNA_indices")
    if isinstance(bt_mt_map, dict) and bt_mt_map:
        return dict(bt_mt_map)
    return {}

MTDNA_REF_MAP = _mtdna_reference_map_from_cfg(CFG)

def _host_reference_map_from_cfg(cfg):
    """
    Return species -> FASTA path map used by soft-clip/CRAM and qualimap host rules.
    Prefer bwa_indices (FASTA paths). If absent, fallback to bowtie2_indices by
    trying common FASTA suffixes next to the Bowtie2 prefix.
    """
    bwa_map = cfg.get("bwa_indices")
    if isinstance(bwa_map, dict) and bwa_map:
        return dict(bwa_map)

    bt_map = cfg.get("bowtie2_indices")
    out = {}
    if isinstance(bt_map, dict):
        for species, prefix in bt_map.items():
            p = str(prefix).strip()
            if not p:
                continue
            candidates = [p + ext for ext in (".fa", ".fasta", ".fna", ".fas")]
            chosen = next((c for c in candidates if os.path.exists(c)), "")
            if chosen:
                out[species] = chosen
    return out

HOST_REF_MAP = _host_reference_map_from_cfg(CFG)

def get_host_ref(wc):
    """Look up host reference based on species file."""
    species_file = f"results/{wc.sample}/fastq_screen/{wc.sample}_best_species.txt"
    # Defer reading until the file exists to avoid DAG build-time crashes
    try:
        if not os.path.exists(species_file):
            return ""
        with open(species_file) as f:
            species = f.read().strip()
    except FileNotFoundError:
        return ""
    return HOST_REF_MAP.get(species, "")
rule softclip_bam_host:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam",
        # DamageProfiler + Qualimap finish before soft-clipping (same dedup BAM; sequential QC order).
        damage = "results/libraries/{sample}/damageprofiler_host/",
        qualimap = "results/libraries/{sample}/qualimap/genome_results.txt",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        cram = "results/libraries/{sample}/host_mapping/{sample}.dedup_q30_softclipped.cram"
    threads: 4
    conda: "workflow/envs/soft_clip.yaml"  # only needs python=2.7 + samtools
    params:
        ref_map = json.dumps(HOST_REF_MAP)
    shell:
        r"""
        set -eo pipefail
        # Resolve reference from best-species file via bwa_indices mapping
        ref=$(python - << 'PY'
import json
m = json.loads(r'''{params.ref_map}''')
species = open("{input.species_file}").read().strip()
print(m.get(species, ""))
PY
        )
        # If no reference was found (e.g. species == "No species found"), create an empty-but-valid CRAM and exit cleanly.
        if [ -z "$ref" ]; then
            echo "[softclip_bam_host] No host reference for species in {input.species_file}; creating stub CRAM." >&2
            samtools view -H {input.bam} | samtools view -O CRAM -o {output.cram}
            exit 0
        fi

        samtools view -h {input.bam} \
        | python2 scripts/softclip_mod.py - {threads} \
        | samtools view -@ {threads} -T "$ref" -O CRAM -o {output.cram}
        """


rule qualimap_bamqc_host_mapping:
    input:
        # Mapping QC on dedup BAM (before soft-clipping).
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam",
        damage = "results/libraries/{sample}/damageprofiler_host/",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        txt = "results/libraries/{sample}/qualimap/genome_results.txt"
    log:
        "logs/qualimap/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"  # needs qualimap + samtools + java
    params:
        ref_map = json.dumps(HOST_REF_MAP)
    shell:
        r"""
        set -eo pipefail
        qdir="results/libraries/{wildcards.sample}/qualimap"
        mkdir -p "$(dirname "$qdir")"

        species_ref=$(python - << 'PY'
import json
m = json.loads(r'''{params.ref_map}''')
species = open("{input.species_file}", "r", encoding="utf-8", errors="ignore").read().strip()
ref = m.get(species, "")
print(species)
print(ref)
PY
        )

        species=$(echo "$species_ref" | sed -n '1p')
        ref=$(echo "$species_ref" | sed -n '2p')

        # If no reference (e.g. species == "No species found"), create stub genome_results.txt and exit cleanly.
        if [ -z "$ref" ]; then
            echo "[qualimap_bamqc_host_mapping] No host reference for species '$species' in {input.species_file}; creating stub genome_results.txt." >> {log}
            mkdir -p "$qdir"
            echo "No host reference found for species $species (file: {input.species_file}). Skipping qualimap_host." > {output.txt}
            exit 0
        fi

        # Check if BAM has mapped reads (qualimap fails on empty BAMs)
        mapped_count=$(samtools view -c -F 4 {input.bam})

        if [ "$mapped_count" -eq 0 ]; then
            echo "No mapped reads in {input.bam} for species $species. Skipping qualimap_host." >> {log}
            mkdir -p "$qdir"
            echo "No mapped reads found in {input.bam} for species $species" > {output.txt}
        else
            # Ensure FASTA index exists for -T usage
            if [ ! -s "$ref.fai" ]; then samtools faidx "$ref"; fi

            rm -rf "$qdir"
            mkdir -p "$qdir"
            samtools view -hb -T "$ref" {input.bam} \
            | qualimap --java-mem-size=9G bamqc -bam /dev/stdin -outdir "$qdir" > {log} 2>&1

            # Defensive check
            if [ ! -f {output.txt} ]; then
                echo "Qualimap did not produce {output.txt}" >&2
                exit 1
            fi
        fi
        """

####-----------------------------------------------------########mtDNA mapping#####----------------------------------------------------------------------------
rule bwa_aln_mtdna:
    input:
        reads = "results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz",
        species = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4.bam"
    threads: 6
    conda:
        "workflow/envs/bwa.yaml"
    params:
        lb   = lambda wc: PCR_INFO[wc.sample]["RGLB"],
        rgid = lambda wc: PCR_INFO[wc.sample]["RGLB"]
    script:
        "scripts/bwa_aln_mtdna.py"


rule sort_bam_initial_mtdna_mapping:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4.bam"
    output:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"


rule filter_q30_mtdna_mapping:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4.sorted.bam"
    output:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools view -q 30 -@ {threads} -b {input.bam} -o {output.bam}"

rule index_q30_bam_mtdna:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam"
    output:
        bai = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule sort_q30_bam_mtdna_mapping:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.bam"
    output:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"


rule mark_duplicates_mtdna:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam"
    output:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam",
        metrics = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.metrics.txt"
    threads: MARKDUP_THREADS
    conda: "workflow/envs/picard.yaml"
    params:
        dedup_tool=DEDUP_TOOL
    shell:
        """
        if [ "{params.dedup_tool}" = "samtools" ]; then
            tmp_sort="results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}.markdup_tmp.bam"
            samtools sort -@ {threads} -o "$tmp_sort" {input.bam}
            samtools index "$tmp_sort"
            samtools markdup -r -@ {threads} "$tmp_sort" {output.bam}
            echo -e "tool\\treads\\tnotes\nsamtools_markdup\\tNA\\tNo Picard-style metrics available" > {output.metrics}
            rm -f "$tmp_sort" "$tmp_sort.bai"
        else
            {PICARD_CMD} MarkDuplicates \
                I={input.bam} \
                O={output.bam} \
                M={output.metrics} \
                REMOVE_DUPLICATES=true \
                ASSUME_SORTED=true \
                VALIDATION_STRINGENCY=SILENT
        fi
        """
rule index_dedup_bam_mtdna:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam"
    output:
        bai = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"

rule save_mtdna_q30_metrics:
    input:
        q30_bam = "results/libraries/{sample}/mtdna_mapping/{sample}_F4_q30.sorted.bam",
        dedup_bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam"
    output:
        metrics_file = "results/libraries/{sample}/mtdna_mapping/{sample}.q30_metrics.txt"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        """
        # Count Q30 reads and dedup reads, then calculate duplication rate
        q30_count=$(samtools view -c -F 4 {input.q30_bam} 2>/dev/null || echo "0")
        dedup_count=$(samtools view -c -F 4 {input.dedup_bam} 2>/dev/null || echo "0")
        
        # Calculate duplication rate
        if [ "$q30_count" -gt 0 ]; then
            dup_rate=$(echo "scale=4; 1 - ($dedup_count / $q30_count)" | bc 2>/dev/null || echo "NA")
        else
            dup_rate="NA"
        fi
        
        # Save to file
        echo "q30_reads=$q30_count" > {output.metrics_file}
        echo "dedup_reads=$dedup_count" >> {output.metrics_file}
        echo "duplication_rate=$dup_rate" >> {output.metrics_file}
        """

rule damage_profiler_mtdna:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam",
        bai = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam.bai",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt",
        config = "config/config.yaml"  # Correct path here
    output:
        dir = directory("results/libraries/{sample}/damageprofiler_mtdna/")
    log:
        "logs/damageprofiler_mtdna/{sample}.log"
    conda:
        "workflow/envs/damageprofiler.yaml"
    params:
        # Match mtDNA reference map to selected aligner.
        index_key=("bowtie2_mtDNA_indices" if HOST_ALIGNER == "bowtie2" else "mtDNA_indices")
    script:
        "scripts/run_damageprofiler.py"
import yaml

# Load config once globally
with open("config/config.yaml") as f:
    CFG = yaml.safe_load(f)

def get_mtDNA_ref(wc):
    """Look up mtDNA reference based on species file."""
    species_file = f"results/{wc.sample}/fastq_screen/{wc.sample}_best_species.txt"
    # Defer reading until the file exists to avoid DAG build-time crashes
    try:
        if not os.path.exists(species_file):
            return ""
        with open(species_file) as f:
            species = f.read().strip()
    except FileNotFoundError:
        return ""
    return MTDNA_REF_MAP.get(species, "")
rule softclip_bam_mtdna:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam",
        # DamageProfiler + Qualimap finish before soft-clipping (same dedup BAM; sequential QC order).
        damage = "results/libraries/{sample}/damageprofiler_mtdna/",
        qualimap = "results/libraries/{sample}/qualimap_mtdna/genome_results.txt",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        cram = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup_q30_softclipped.cram"
    threads: 4
    conda: "workflow/envs/soft_clip.yaml"  # only needs python=2.7 + samtools
    params:
        ref_map = json.dumps(MTDNA_REF_MAP)
    shell:
        r"""
        set -eo pipefail
        mkdir -p "$(dirname '{output.cram}')"
        # Resolve FASTA for samtools -T (mtDNA_indices may list Bowtie2 prefix; use sibling .fa if needed)
        ref=$(python - << 'PY'
import codecs
import json
import os

def pick_fasta(p):
    if not p:
        return p
    if os.path.isfile(p) and os.path.getsize(p) > 0:
        return p
    base = p.rstrip("/")
    for suf in (".fa", ".fasta", ".fna"):
        c = base + suf
        if os.path.isfile(c) and os.path.getsize(c) > 0:
            return c
    return p

def ref_lookup(m, species):
    if species in m:
        return m[species] or ""
    sl = species.strip().lower()
    for k, v in m.items():
        if str(k).strip().lower() == sl:
            return v or ""
    return ""

m = json.loads(r'''{params.ref_map}''')
species = codecs.open("{input.species_file}", "r", encoding="utf-8", errors="ignore").read().strip()
print(pick_fasta(ref_lookup(m, species)))
PY
        )
        # If no reference was found (e.g. species == "No species found"), create an empty-but-valid CRAM and exit cleanly.
        if [ -z "$ref" ]; then
            echo "[softclip_bam_mtdna] No mtDNA reference for species in {input.species_file}; creating stub CRAM." >&2
            samtools view -H {input.bam} | samtools view -O CRAM -o {output.cram}
            exit 0
        fi

        samtools view -h {input.bam} \
        | python2 scripts/softclip_mod.py - {threads} \
        | samtools view -@ {threads} -T "$ref" -O CRAM -o {output.cram}
        """


# -------------------- Per-sample merged BAMs (host, mtDNA, pathogen) -----------------------------

def host_dedup_bams_for_sample(wc):
    """All host dedup BAMs for PCRs belonging to a biological sample."""
    return expand(
        "results/libraries/{pcr}/host_mapping/{pcr}.dedup.bam",
        pcr=pcrs_for_sample(wc.sample),
    )


def host_dedup_metrics_for_sample(wc):
    """All host dedup metrics files for PCRs belonging to a biological sample."""
    return expand(
        "results/libraries/{pcr}/host_mapping/{pcr}.dedup.metrics.txt",
        pcr=pcrs_for_sample(wc.sample),
    )


def host_species_files_for_sample(wc):
    """All species files for PCRs belonging to a biological sample (for species mismatch check)."""
    return expand(
        "results/libraries/{pcr}/fastq_screen/{pcr}_best_species.txt",
        pcr=pcrs_for_sample(wc.sample),
    )


rule merge_host_dedup_per_sample:
    input:
        bams=host_dedup_bams_for_sample,
        metrics_files=host_dedup_metrics_for_sample,
        species_files=host_species_files_for_sample
    output:
        merged="results/host/host_mapping/{sample}.dedup.merged.bam",
        metrics="results/host/host_mapping/{sample}.dedup.merged.metrics.txt",
        warning="results/host/host_mapping/{sample}.species_mismatch_warning.txt"
    conda:
        "workflow/envs/picard.yaml"
    threads: 6
    params:
        merge_tool=MERGE_TOOL
    shell:
        """
        mkdir -p results/host/host_mapping
        
        # Check for species mismatch across PCRs
        species_list=""
        pcr_list=""
        mismatch_detected=false
        
        for species_file in {input.species_files}; do
            if [ -f "$species_file" ]; then
                species=$(cat "$species_file" | head -1 | tr -d '[:space:]')
                pcr=$(basename $(dirname $(dirname "$species_file")))
                if [ -n "$species" ] && [ "$species" != "NA" ] && [ "$species" != "No species found" ]; then
                    if [ -z "$species_list" ]; then
                        species_list="$species"
                        pcr_list="$pcr"
                    elif [ "$species_list" != "$species" ]; then
                        mismatch_detected=true
                        echo "WARNING: Species mismatch detected for sample {wildcards.sample}" > {output.warning}
                        echo "PCRs assigned to different species:" >> {output.warning}
                        echo "  $pcr_list: $species_list" >> {output.warning}
                        echo "  $pcr: $species" >> {output.warning}
                        # Check all remaining PCRs to list all mismatches
                        for remaining_file in {input.species_files}; do
                            if [ "$remaining_file" != "$species_file" ]; then
                                remaining_species=$(cat "$remaining_file" 2>/dev/null | head -1 | tr -d '[:space:]' || echo "NA")
                                remaining_pcr=$(basename $(dirname $(dirname "$remaining_file")))
                                if [ -n "$remaining_species" ] && [ "$remaining_species" != "NA" ] && [ "$remaining_species" != "No species found" ]; then
                                    echo "  $remaining_pcr: $remaining_species" >> {output.warning}
                                fi
                            fi
                        done
                        echo "" >> {output.warning}
                        echo "Host/mtDNA BAM merging skipped for this sample due to species mismatch." >> {output.warning}
                        break
                    fi
                fi
            fi
        done
        
        if [ "$mismatch_detected" = "true" ]; then
            # Species mismatch: create stub outputs and skip merging
            echo "[WARNING] Species mismatch detected for {wildcards.sample}. Skipping host BAM merge." >&2
            cat {output.warning} >&2
            : > {output.merged}
            echo "Species mismatch - merge skipped" > {output.metrics}
        elif [ $(echo {input.bams} | wc -w) -eq 1 ]; then
            # Only one PCR BAM: just link it as the merged BAM and copy its metrics file
            ln -sf $(realpath {input.bams}) {output.merged}
            # Copy the metrics file from the single PCR
            first_metrics=$(echo {input.metrics_files} | awk '{{print $1}}')
            if [ -f "$first_metrics" ] && [ -s "$first_metrics" ]; then
                cp "$first_metrics" {output.metrics}
            else
                : > {output.metrics}
            fi
            : > {output.warning}  # No warning for single PCR
        else
            # Multiple PCRs with consistent species: proceed with merge
            if [ "{params.merge_tool}" = "samtools" ]; then
                samtools merge -@ {threads} {output.merged} {input.bams}
            else
                picard MergeSamFiles \
                    $(for bam in {input.bams}; do echo -n " INPUT=$bam"; done) \
                    OUTPUT={output.merged} \
                    ASSUME_SORTED=true \
                    VALIDATION_STRINGENCY=SILENT
            fi
            # Copy metrics from the first PCR (or aggregate if needed in the future)
            first_metrics=$(echo {input.metrics_files} | awk '{{print $1}}')
            if [ -f "$first_metrics" ] && [ -s "$first_metrics" ]; then
                cp "$first_metrics" {output.metrics}
            else
                : > {output.metrics}
            fi
            : > {output.warning}  # No warning
        fi
        """


rule index_host_merged_bam:
    input:
        bam="results/host/host_mapping/{sample}.dedup.merged.bam",
    output:
        bai="results/host/host_mapping/{sample}.dedup.merged.bam.bai",
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        """
        if [ ! -s {input.bam} ]; then
            : > {output.bai}
            exit 0
        fi
        samtools index {input.bam}
        """


rule sexing_residual:
    """Chromosome-residual sexing on merged host BAM per biological sample (Cow, Goat, Sheep, Dog)."""
    input:
        bam="results/host/host_mapping/{sample}.dedup.merged.bam",
        bai="results/host/host_mapping/{sample}.dedup.merged.bam.bai",
        species_files=host_species_files_for_sample,
        species_mismatch_warning="results/host/host_mapping/{sample}.species_mismatch_warning.txt",
        r_script=SEXING_R_SCRIPT,
    output:
        idx="results/samples/{sample}/sexing/{sample}_sexing.idx",
        pdf="results/samples/{sample}/sexing/{sample}_sexing.pdf",
        tsv="results/samples/{sample}/sexing/{sample}_sexing.tsv",
    params:
        enabled=1 if ENABLE_SEXING else 0,
        host_ref_map=json.dumps(HOST_REF_MAP),
    conda:
        "workflow/envs/sexing.yaml"
    log:
        "logs/sexing/{sample}.log"
    script:
        "scripts/run_pigsti_sexing.py"


def mtdna_dedup_bams_for_sample(wc):
    """All mtDNA dedup BAMs for PCRs belonging to a biological sample."""
    return expand(
        "results/libraries/{pcr}/mtdna_mapping/{pcr}.dedup.bam",
        pcr=pcrs_for_sample(wc.sample),
    )


def mtdna_dedup_metrics_for_sample(wc):
    """All mtDNA dedup metrics files for PCRs belonging to a biological sample."""
    return expand(
        "results/libraries/{pcr}/mtdna_mapping/{pcr}.dedup.metrics.txt",
        pcr=pcrs_for_sample(wc.sample),
    )


def mtdna_species_files_for_sample(wc):
    """All species files for PCRs belonging to a biological sample (for species mismatch check)."""
    return expand(
        "results/libraries/{pcr}/fastq_screen/{pcr}_best_species.txt",
        pcr=pcrs_for_sample(wc.sample),
    )


rule merge_mtdna_dedup_per_sample:
    input:
        bams=mtdna_dedup_bams_for_sample,
        metrics_files=mtdna_dedup_metrics_for_sample,
        species_files=mtdna_species_files_for_sample
    output:
        merged="results/host/mtdna_mapping/{sample}.dedup.merged.bam",
        metrics="results/host/mtdna_mapping/{sample}.dedup.merged.metrics.txt",
        warning="results/host/mtdna_mapping/{sample}.species_mismatch_warning.txt"
    conda:
        "workflow/envs/picard.yaml"
    threads: 6
    params:
        merge_tool=MERGE_TOOL
    shell:
        """
        mkdir -p results/host/mtdna_mapping
        
        # Check for species mismatch across PCRs (same logic as host merge)
        species_list=""
        pcr_list=""
        mismatch_detected=false
        
        for species_file in {input.species_files}; do
            if [ -f "$species_file" ]; then
                species=$(cat "$species_file" | head -1 | tr -d '[:space:]')
                pcr=$(basename $(dirname $(dirname "$species_file")))
                if [ -n "$species" ] && [ "$species" != "NA" ] && [ "$species" != "No species found" ]; then
                    if [ -z "$species_list" ]; then
                        species_list="$species"
                        pcr_list="$pcr"
                    elif [ "$species_list" != "$species" ]; then
                        mismatch_detected=true
                        echo "WARNING: Species mismatch detected for sample {wildcards.sample}" > {output.warning}
                        echo "PCRs assigned to different species:" >> {output.warning}
                        echo "  $pcr_list: $species_list" >> {output.warning}
                        echo "  $pcr: $species" >> {output.warning}
                        # Check all remaining PCRs to list all mismatches
                        for remaining_file in {input.species_files}; do
                            if [ "$remaining_file" != "$species_file" ]; then
                                remaining_species=$(cat "$remaining_file" 2>/dev/null | head -1 | tr -d '[:space:]' || echo "NA")
                                remaining_pcr=$(basename $(dirname $(dirname "$remaining_file")))
                                if [ -n "$remaining_species" ] && [ "$remaining_species" != "NA" ] && [ "$remaining_species" != "No species found" ]; then
                                    echo "  $remaining_pcr: $remaining_species" >> {output.warning}
                                fi
                            fi
                        done
                        echo "" >> {output.warning}
                        echo "Host/mtDNA BAM merging skipped for this sample due to species mismatch." >> {output.warning}
                        break
                    fi
                fi
            fi
        done
        
        if [ "$mismatch_detected" = "true" ]; then
            # Species mismatch: create stub outputs and skip merging
            echo "[WARNING] Species mismatch detected for {wildcards.sample}. Skipping mtDNA BAM merge." >&2
            cat {output.warning} >&2
            : > {output.merged}
            echo "Species mismatch - merge skipped" > {output.metrics}
        elif [ $(echo {input.bams} | wc -w) -eq 1 ]; then
            # Only one PCR BAM: just link it as the merged BAM and copy its metrics file
            ln -sf $(realpath {input.bams}) {output.merged}
            # Copy the metrics file from the single PCR
            first_metrics=$(echo {input.metrics_files} | awk '{{print $1}}')
            if [ -f "$first_metrics" ] && [ -s "$first_metrics" ]; then
                cp "$first_metrics" {output.metrics}
            else
                : > {output.metrics}
            fi
            : > {output.warning}  # No warning for single PCR
        else
            # Multiple PCRs with consistent species: proceed with merge
            if [ "{params.merge_tool}" = "samtools" ]; then
                samtools merge -@ {threads} {output.merged} {input.bams}
            else
                picard MergeSamFiles \
                    $(for bam in {input.bams}; do echo -n " INPUT=$bam"; done) \
                    OUTPUT={output.merged} \
                    ASSUME_SORTED=true \
                    VALIDATION_STRINGENCY=SILENT
            fi
            # Copy metrics from the first PCR (or aggregate if needed in the future)
            first_metrics=$(echo {input.metrics_files} | awk '{{print $1}}')
            if [ -f "$first_metrics" ] && [ -s "$first_metrics" ]; then
                cp "$first_metrics" {output.metrics}
            else
                : > {output.metrics}
            fi
            : > {output.warning}  # No warning
        fi
        """


def pathogen_dedup_bams_for_sample(wc):
    """
    All pathogen dedup BAMs for PCRs belonging to a biological sample and a given pathogen
    (ref_name_safe wildcard).
    """
    return expand(
        "results/libraries/{pcr}/pathogen_mapping/{pcr}_{ref}.dedup.bam",
        pcr=pcrs_for_sample(wc.sample),
        ref=wc.ref_name_safe,
    )


rule merge_pathogen_dedup_per_sample:
    input:
        pathogen_dedup_bams_for_sample
    output:
        merged="results/pools/pathogen_mapping/{sample}_{ref_name_safe}.dedup.merged.bam",
        metrics="results/pools/pathogen_mapping/{sample}_{ref_name_safe}.dedup.merged.metrics.txt"
    conda:
        "workflow/envs/picard.yaml"
    threads: 6
    params:
        merge_tool=MERGE_TOOL
    shell:
        """
        mkdir -p results/pools/pathogen_mapping
        if [ $(echo {input} | wc -w) -eq 1 ]; then
            # Only one PCR BAM: just link it as the merged BAM and create an empty metrics file
            ln -sf $(realpath {input}) {output.merged}
            : > {output.metrics}
        else
            if [ "{params.merge_tool}" = "samtools" ]; then
                samtools merge -@ {threads} {output.merged} {input}
            else
                picard MergeSamFiles \
                    $(for bam in {input}; do echo -n " INPUT=$bam"; done) \
                    OUTPUT={output.merged} \
                    ASSUME_SORTED=true \
                    VALIDATION_STRINGENCY=SILENT
            fi
            : > {output.metrics}
        fi
        """
rule host_bam_by_type:
    input:
        bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam"
    output:
        bam = f"{BY_TYPE_ROOT}/host_mapping/{{sample}}/{{sample}}.dedup.bam"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/host_mapping/{wildcards.sample}
        rm -f {output.bam} 2>/dev/null || true
        ln -sf $(realpath {input.bam}) {output.bam}
        """


rule mtdna_bam_by_type:
    input:
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam"
    output:
        bam = f"{BY_TYPE_ROOT}/mtdna_mapping/{{sample}}/{{sample}}.dedup.bam"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/mtdna_mapping/{wildcards.sample}
        rm -f {output.bam} 2>/dev/null || true
        ln -sf $(realpath {input.bam}) {output.bam}
        """


rule pathogen_bam_by_type:
    input:
        bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
    output:
        bam = f"{BY_TYPE_ROOT}/pathogen_mapping/{{sample}}/{{sample}}_{{ref_name_safe}}.dedup.bam"
    shell:
        """
        mkdir -p {BY_TYPE_ROOT}/pathogen_mapping/{wildcards.sample}
        rm -f {output.bam} 2>/dev/null || true
        ln -sf $(realpath {input.bam}) {output.bam}
        """


rule qualimap_bamqc_mtdna_mapping:
    input:
        # Mapping QC on dedup BAM (before soft-clipping).
        bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam",
        damage = "results/libraries/{sample}/damageprofiler_mtdna/",
        species_file = "results/libraries/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        txt = "results/libraries/{sample}/qualimap_mtdna/genome_results.txt"
    log:
        "logs/qualimap_mtdna/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"  # needs qualimap + samtools + java
    params:
        ref_map = json.dumps(MTDNA_REF_MAP)
    shell:
        r"""
        set -eo pipefail
        qdir="results/libraries/{wildcards.sample}/qualimap_mtdna"
        mkdir -p "$(dirname "$qdir")"

        species_ref=$(python - << 'PY'
import json
import os

def pick_fasta(p):
    if not p:
        return p
    if os.path.isfile(p) and os.path.getsize(p) > 0:
        return p
    base = p.rstrip("/")
    for suf in (".fa", ".fasta", ".fna"):
        c = base + suf
        if os.path.isfile(c) and os.path.getsize(c) > 0:
            return c
    return p

def ref_lookup(m, species):
    if species in m:
        return m[species] or ""
    sl = species.strip().lower()
    for k, v in m.items():
        if str(k).strip().lower() == sl:
            return v or ""
    return ""

m = json.loads(r'''{params.ref_map}''')
species = open("{input.species_file}", "r", encoding="utf-8", errors="ignore").read().strip()
ref = pick_fasta(ref_lookup(m, species))
print(species)
print(ref)
PY
        )

        species=$(echo "$species_ref" | sed -n '1p')
        ref=$(echo "$species_ref" | sed -n '2p')

        # If no reference (e.g. species == "No species found"), create stub genome_results.txt and exit cleanly.
        if [ -z "$ref" ]; then
            echo "[qualimap_bamqc_mtdna_mapping] No mtDNA reference for species '$species' in {input.species_file}; creating stub genome_results.txt." >> {log}
            mkdir -p "$qdir"
            echo "No mtDNA reference found for species $species (file: {input.species_file}). Skipping qualimap_mtdna." > {output.txt}
            exit 0
        fi

        # Check if BAM has mapped reads (qualimap fails on empty BAMs)
        mapped_count=$(samtools view -c -F 4 {input.bam})

        if [ "$mapped_count" -eq 0 ]; then
            echo "No mapped reads in {input.bam} for species $species. Skipping qualimap_mtdna." >> {log}
            mkdir -p "$qdir"
            echo "No mapped reads found in {input.bam} for species $species" > {output.txt}
        else
            # Ensure FASTA index exists for -T usage
            if [ ! -s "$ref.fai" ]; then samtools faidx "$ref"; fi

            rm -rf "$qdir"
            mkdir -p "$qdir"
            samtools view -hb -T "$ref" {input.bam} \
            | qualimap --java-mem-size=9G bamqc -bam /dev/stdin -outdir "$qdir" > {log} 2>&1

            # Defensive check
            if [ ! -f {output.txt} ]; then
                echo "Qualimap did not produce {output.txt}" >&2
                exit 1
            fi
        fi
        """

# Qualimap on merged per-sample host BAMs (for sample-level summary)
rule qualimap_bamqc_host_mapping_merged:
    input:
        bam = "results/host/host_mapping/{sample}.dedup.merged.bam",
        warning = "results/host/host_mapping/{sample}.species_mismatch_warning.txt"
    output:
        txt = "results/samples/{sample}/qualimap/genome_results.txt"
    log:
        "logs/qualimap_sample/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"
    shell:
        """
        set -eo pipefail
        qdir="results/samples/{wildcards.sample}/qualimap"
        
        # Check for species mismatch warning - skip qualimap if present
        if [ -s {input.warning} ]; then
            echo "[WARNING] Species mismatch detected for {wildcards.sample}. Skipping Qualimap on merged host BAM." > {log}
            mkdir -p "$qdir"
            echo "Species mismatch - Qualimap skipped" > {output.txt}
            exit 0
        fi
        
        # Check if BAM is empty (0 bytes or very small).
        # Use stat -Lc%s to dereference symlinks and get the target size.
        if [ ! -s {input.bam} ] || [ $(stat -Lc%s {input.bam} 2>/dev/null || echo 0) -lt 1000 ]; then
            echo "[WARNING] Empty or very small BAM file: {input.bam}. Creating stub output." > {log}
            mkdir -p "$qdir"
            echo "No mapped reads found in {input.bam}" > {output.txt}
            exit 0
        fi
        
        rm -rf "$qdir"
        mkdir -p "$qdir"
        # Run qualimap but don't fail if it exits non-zero (e.g., empty BAM)
        if ! qualimap bamqc \
            -bam {input.bam} \
            -outdir "$qdir" \
            -outformat html > {log} 2>&1; then
            echo "[WARNING] Qualimap failed for {input.bam} (likely no mapped reads). Creating stub output." >> {log}
        fi
        # Qualimap writes genome_results.txt directly in the outdir (same as output path)
        if [ ! -f {output.txt} ]; then
            echo "No mapped reads found in {input.bam}" > {output.txt}
        fi
        """

# Qualimap on merged per-sample mtDNA BAMs (for sample-level summary)
rule qualimap_bamqc_mtdna_mapping_merged:
    input:
        bam = "results/host/mtdna_mapping/{sample}.dedup.merged.bam",
        warning = "results/host/mtdna_mapping/{sample}.species_mismatch_warning.txt"
    output:
        txt = "results/samples/{sample}/qualimap_mtdna/genome_results.txt"
    log:
        "logs/qualimap_sample_mtdna/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"
    shell:
        r"""
        set -eo pipefail
        qdir="results/samples/{wildcards.sample}/qualimap_mtdna"
        
        # Check for species mismatch warning - skip qualimap if present
        if [ -s {input.warning} ]; then
            echo "[WARNING] Species mismatch detected for {wildcards.sample}. Skipping Qualimap on merged mtDNA BAM." > {log}
            mkdir -p "$qdir"
            echo "Species mismatch - Qualimap skipped" > {output.txt}
            exit 0
        fi
        
        # Check if BAM is empty (0 bytes or very small).
        # Use stat -Lc%s to dereference symlinks and get the target size.
        if [ ! -s {input.bam} ] || [ $(stat -Lc%s {input.bam} 2>/dev/null || echo 0) -lt 1000 ]; then
            echo "[WARNING] Empty or very small BAM file: {input.bam}. Creating stub output." > {log}
            mkdir -p "$qdir"
            echo "No mapped reads found in {input.bam}" > {output.txt}
            exit 0
        fi
        
        rm -rf "$qdir"
        mkdir -p "$qdir"
        # Run qualimap, but tolerate failures (e.g. no mapped reads) and fall back to a stub QC file.
        if ! qualimap bamqc \
            -bam {input.bam} \
            -outdir "$qdir" \
            -outformat html > {log} 2>&1; then
            echo "[qualimap_bamqc_mtdna_mapping_merged] qualimap failed for {input.bam}, creating stub genome_results.txt" >> {log} 2>&1
        fi
        # Qualimap writes genome_results.txt directly in the outdir (same as output path).
        if [ ! -f {output.txt} ]; then
            echo "No mapped reads found in {input.bam}" > {output.txt}
        fi
        """


###----------------------------------------wrappers------------------------------------------------######################
rule merge_pathogen_summaries:
    input:
        summaries=expand("results/pathogen/{sample}/summary/{sample}_pathogen_summary.csv", sample=BIO_SAMPLES)
    output:
        excel="results/final/pathogen_summary_all_samples.xlsx"
    conda:
        "workflow/envs/summary.yaml"
    script:
        "scripts/merge_summaries.py"

# Create heatmap of detection scores for all pathogens across all samples
rule create_pathogen_heatmap:
    input:
        excel="results/final/pathogen_summary_all_samples.xlsx"
    output:
        png="results/final/pathogen_detection_scores_heatmap.png",
        pdf="results/final/pathogen_detection_scores_heatmap.pdf"
    conda:
        "workflow/envs/summary.yaml"
    script:
        "scripts/create_pathogen_heatmap.py"

# Create comprehensive summary Excel with all samples
if HOST_MTDNA_ANALYSIS_ENABLED:
    rule create_comprehensive_summary:
        input:
            host_mtdna="results/final/host_mtdna_summary_all_samples.xlsx",
            pathogen="results/final/pathogen_summary_all_samples.xlsx"
        output:
            excel="results/final/comprehensive_summary_all_samples.xlsx"
        conda:
            "workflow/envs/summary.yaml"
        script:
            "scripts/create_comprehensive_summary.py"
else:
    # Screening-only / pathogen-only mode: still produce a "comprehensive" spreadsheet
    # with basic sample/pathogen metrics, even if host mapping/QC is disabled.
    rule create_comprehensive_summary:
        input:
            pathogen="results/final/pathogen_summary_all_samples.xlsx"
        output:
            excel="results/final/comprehensive_summary_all_samples.xlsx"
        conda:
            "workflow/envs/summary.yaml"
        script:
            "scripts/create_comprehensive_summary_pathogen_only.py"

# Generate per-sample PDF reports
rule generate_sample_report:
    input:
        host_mtdna="results/final/host_mtdna_summary_all_samples.xlsx",
        pathogen="results/final/pathogen_summary_all_samples.xlsx"
    output:
        pdf="results/pathogen/{sample}/summary/{sample}_sample_report.pdf"
    params:
        animal_icons_dir="config/animal_icons"  # Directory containing animal SVG icons
    conda:
        "workflow/envs/report_env.yaml"
    shell:
        """
        Rscript scripts/generate_sample_report.R {wildcards.sample} {input.host_mtdna} {input.pathogen} {output.pdf} {params.animal_icons_dir}
        """

def refs_for_sample(sample):
    pairs = get_sample_ref_pairs()
    return [safe_name(ref) for s, ref in pairs if s == sample]

def summarize_inputs(wc):
    refs = refs_for_sample(wc.sample)
    # Instead of specific files inside these directories, just input the directories:
    damage_dirs = [f"results/pathogen/{wc.sample}/pathogen_mapping/damageprofiler_{ref}" for ref in refs]
    qualimap_dirs = [f"results/pathogen/{wc.sample}/pathogen_mapping/qualimap_{ref}" for ref in refs]

    return [
        f"results/{wc.sample}/evalue/pathogen/{wc.sample}_pathogen.csv",
        "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv",
        "config/Pathogen_spreadsheet.csv",
    ] + damage_dirs + qualimap_dirs

if ENABLE_HOPS:
    rule summarize_pathogen_data:
        input:
            pathogen_complete="results/workflow/pathogen_mapping_complete.txt",
            evalue = "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
            hops = "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv",
            spreadsheet = "config/Pathogen_spreadsheet.csv",
            qualimaps = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/qualimap_{ps}"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            damage = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/damageprofiler_{ps}"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            ani = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.ani.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            breadth = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.breadth.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.mean_entropy.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy_100 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.entropy_100bp.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy_1000 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.entropy_1000bp.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            edit_r2 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.edit_distance_logr2.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            genus_ranking = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.genus_ranking.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
        output:
            "results/pathogen/{sample}/summary/{sample}_pathogen_summary.csv"
        params:
            cleanup_flag = "--cleanup" if CLEANUP_ENABLED else ""
        conda:
            "workflow/envs/summary.yaml"
        shell:
            """
            python scripts/summarize_pathogen_data.py \
                --sample {wildcards.sample} \
                --evalue {input.evalue} \
                --hops {input.hops} \
                --spreadsheet {input.spreadsheet} \
                --config config/config.yaml \
                --bam_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --qualimap_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --damage_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --output {output} \
                {params.cleanup_flag}
            """
else:
    # HOPS disabled: run pathogen summary using Escore + BWA metrics only
    rule summarize_pathogen_data:
        input:
            pathogen_complete="results/workflow/pathogen_mapping_complete.txt",
            evalue = "results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
            spreadsheet = "config/Pathogen_spreadsheet.csv",
            qualimaps = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/qualimap_{ps}"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            damage = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/damageprofiler_{ps}"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            ani = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.ani.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            breadth = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.breadth.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.mean_entropy.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy_100 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.entropy_100bp.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            entropy_1000 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.entropy_1000bp.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            edit_r2 = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.edit_distance_logr2.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
            genus_ranking = lambda wildcards: [
                f"results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{ps}.genus_ranking.txt"
                for ps in pathogen_safe_names_for_sample(wildcards.sample)
            ],
        output:
            "results/pathogen/{sample}/summary/{sample}_pathogen_summary.csv"
        params:
            cleanup_flag = "--cleanup" if CLEANUP_ENABLED else ""
        conda:
            "workflow/envs/summary.yaml"
        shell:
            """
            python scripts/summarize_pathogen_data.py \
                --sample {wildcards.sample} \
                --evalue {input.evalue} \
                --hops '' \
                --spreadsheet {input.spreadsheet} \
                --config config/config.yaml \
                --bam_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --qualimap_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --damage_dir results/pathogen/{wildcards.sample}/pathogen_mapping \
                --output {output} \
                {params.cleanup_flag}
            """

# PCR-level pathogen summaries (only in super_careful mode)
if PATHOGEN_MODE == "super_careful":
    rule summarize_pathogen_data_pcr:
        input:
            evalue=lambda wc: (
                f"results/pathogen/{PCR_INFO[wc.pcr]['sample']}/evalue/pathogen/"
                f"{PCR_INFO[wc.pcr]['sample']}_pathogen.csv"
            ),
            spreadsheet = "config/Pathogen_spreadsheet.csv",
            qualimaps = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/qualimap_{ps}"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            damage = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/damageprofiler_{ps}"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            ani = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/{wc.pcr}_{ps}.ani.txt"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            breadth = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/{wc.pcr}_{ps}.breadth.txt"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            entropy = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/{wc.pcr}_{ps}.mean_entropy.txt"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            entropy_100 = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/{wc.pcr}_{ps}.entropy_100bp.txt"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
            entropy_1000 = lambda wc: [
                f"results/libraries/{wc.pcr}/pathogen_mapping/{wc.pcr}_{ps}.entropy_1000bp.txt"
                for s, ps in pathogen_pairs_from_checkpoint()
                if s == PCR_INFO[wc.pcr]["sample"]
            ],
        output:
            "results/libraries/{pcr}/summary/{pcr}_pathogen_summary.csv"
        params:
            # Map PCR ID back to bio sample
            sample = lambda wc: PCR_INFO[wc.pcr]["sample"]
        conda:
            "workflow/envs/summary.yaml"
        shell:
            """
            python scripts/summarize_pathogen_data_pcr.py \
                --sample {params.sample} \
                --pcr {wildcards.pcr} \
                --evalue {input.evalue} \
                --spreadsheet {input.spreadsheet} \
                --bam_dir results/libraries/{wildcards.pcr}/pathogen_mapping \
                --qualimap_dir results/libraries/{wildcards.pcr}/pathogen_mapping \
                --damage_dir results/libraries/{wildcards.pcr}/pathogen_mapping \
                --output {output}
            """

    # Merge PCR-level pathogen CSVs into an Excel (PCR sheet)
    rule merge_pathogen_summaries_pcr:
        input:
            summaries = expand("results/libraries/{pcr}/summary/{pcr}_pathogen_summary.csv", pcr=PCRS)
        output:
            excel = "results/final/pathogen_summary_all_samples_pcr.xlsx"
        conda:
            "workflow/envs/summary.yaml"
        script:
            "scripts/merge_summaries.py"


rule summarize_host_mtdna:
    input:
        host_bams = expand("results/libraries/{sample}/host_mapping/{sample}.dedup.bam", sample=SAMPLES),
        mtdna_bams = expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam", sample=SAMPLES),
        qualimap_host = expand("results/libraries/{sample}/qualimap/genome_results.txt", sample=SAMPLES),
        qualimap_mtdna = expand("results/libraries/{sample}/qualimap_mtdna/genome_results.txt", sample=SAMPLES),
        qualimap_host_merged = expand("results/samples/{sample}/qualimap/genome_results.txt", sample=BIO_SAMPLES),
        qualimap_mtdna_merged = expand("results/samples/{sample}/qualimap_mtdna/genome_results.txt", sample=BIO_SAMPLES),
        samples_tsv = SAMPLES_TSV,
        species = expand("results/libraries/{sample}/fastq_screen/{sample}_best_species.txt", sample=SAMPLES),
        collapsed_qc = expand("results/libraries/{sample}/qc/{sample}.collapsed_reads.txt", sample=SAMPLES),
        raw_qc = expand("results/libraries/{sample}/qc/{sample}.raw_reads.txt", sample=SAMPLES),
        sexing_tsv=expand(
            "results/samples/{sample}/sexing/{sample}_sexing.tsv", sample=BIO_SAMPLES
        )
        if ENABLE_SEXING
        else [],
    output:
        "results/final/host_mtdna_summary_all_samples.xlsx"
    params:
        samples = ",".join(SAMPLES)
    threads: SUMMARY_THREADS
    conda:
        "workflow/envs/summary.yaml"
    script:
        "scripts/summarize_host_mtdna.py"


rule pathogen_report:
    input:
        mapping_complete="results/workflow/pathogen_mapping_complete.txt",
        evalue="results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv",
        spreadsheet="config/Pathogen_spreadsheet.csv",
        summary="results/pathogen/{sample}/summary/{sample}_pathogen_summary.csv",
        # Pathogen-specific outputs
        ani="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.ani.txt",
        entropy_plot="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.entropy_plot.png",
        entropy_mean="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.mean_entropy.txt",
        damageprofiler_dir="results/pathogen/{sample}/pathogen_mapping/damageprofiler_{ref_name_safe}",
        # AdnaPlotter removed - too time and memory consuming
        qualimap_dir="results/pathogen/{sample}/pathogen_mapping/qualimap_{ref_name_safe}",
        breadth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.breadth.txt",
        expected_breadth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.expected_breadth.txt",
        breadth_ratio="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.breadth_ratio.txt",
        depth="results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.depth.txt"
    output:
        pdf="results/pathogen/{sample}/summary/{sample}_{ref_name_safe}_pathogen_report.pdf"
    params:
        pathogen=lambda wc: wc.ref_name_safe.replace("_", " ")
    conda:
        "workflow/envs/report_env.yaml"
    shell:
        """
        Rscript scripts/generate_pathogen_report.R {wildcards.sample} '{params.pathogen}' {input.spreadsheet} {output.pdf}
        """

# Generate all pathogen-specific PDF reports (checkpoint-driven).
# This avoids relying on `get_sample_ref_pairs()` at DAG parse time.
rule generate_all_pathogen_reports:
    input:
        mapping_complete="results/workflow/pathogen_mapping_complete.txt",
        pathogen_summaries=expand("results/pathogen/{sample}/summary/{sample}_pathogen_summary.csv", sample=BIO_SAMPLES),
        pdfs=pathogen_report_pdfs_from_checkpoint,
    output:
        touch("results/workflow/pathogen_reports_complete.txt")
    conda:
        "workflow/envs/report_env.yaml"
    shell:
        r"""
        # Reports are generated by the per-target `pathogen_report` rule (one job per PDF).
        # This rule only acts as a barrier/marker once all PDFs exist.
        touch {output}
        """


rule write_run_manifest:
    """
    Write a lightweight reproducibility manifest (file hashes + key config).
    This is intended for paper supplementary information and for long-term reuse.
    """
    input:
        # Ensure the workflow ran far enough that the checkpoint produced targets.
        "results/workflow/pathogen_reports_complete.txt",
        "config/config.yaml",
        "config/samples.tsv",
        "config/Pathogen_spreadsheet.csv",
        # Pathogen detection criteria thresholds (embedded in config/config.yaml; legacy file optional).
        # Keep the legacy file as optional input to avoid breaking existing setups.
        *(["config/pathogen_detection_criteria.yaml"] if os.path.exists("config/pathogen_detection_criteria.yaml") else []),
    output:
        "results/run_manifest.json"
    conda:
        "workflow/envs/python.yaml"
    shell:
        r"""
        set -euo pipefail
        python scripts/write_run_manifest.py --output "{output}" --snakefile "Snakefile" --config "config/config.yaml"
        """


if BUILD_RESULTS_CATALOG:

    rule build_results_catalog:
        """
        Write results/final/output_catalog.tsv and optionally populate results/browse/
        with symlinks/pointers so all Kraken, E-score, etc. files are easy to find together.
        Canonical paths: results/libraries/, results/samples/, results/host/, results/pools/unaligned_fastq/, results/pathogen/, results/metagenomics/, results/final/, results/workflow/.
        """
        input:
            pathogen_complete="results/workflow/pathogen_mapping_complete.txt",
            manifest="results/workflow/pathogen_targets.manifest.json",
        output:
            catalog="results/final/output_catalog.tsv",
        params:
            symlink_flag="--symlinks" if BUILD_RESULTS_CATALOG_SYMLINKS else "--no-symlinks",
        conda:
            "workflow/envs/python.yaml",
        shell:
            """
            python scripts/build_results_catalog.py --config config/config.yaml {params.symlink_flag}
            """


rule create_multi_qc_dashboard:
    """
    Produce a single static HTML dashboard aggregating host/mtDNA QC
    and pathogen summary information per biological sample.
    """
    input:
        host_mtdna="results/final/host_mtdna_summary_all_samples.xlsx",
        pathogen="results/final/pathogen_summary_all_samples.xlsx"
    output:
        html="results/final/multi_qc_dashboard.html"
    conda:
        "workflow/envs/summary.yaml"
    params:
        top_pathogens=MULTI_QC_TOP_PATHOGENS,
        results_dir="results/pathogen"
    script:
        "scripts/create_multi_qc_dashboard.py"

# Pathogen Detection Scoring System
if ENABLE_HOPS:
    rule calculate_pathogen_detection_scores:
        input:
            # E-Score results (already filtered by user thresholds)
            escore_files = expand("results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES),
            # Hops results
            hops_results = "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv",
            # BWA / damageprofiler results for the checkpoint-selected pathogens
            bwa_results = lambda wildcards: [
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{ps}.ani.txt"
                for sample, ps in pathogen_pairs_from_checkpoint()
            ],
            damage_results = lambda wildcards: [
                f"results/pathogen/{sample}/pathogen_mapping/damageprofiler_{ps}"
                for sample, ps in pathogen_pairs_from_checkpoint()
            ],
            breadth_results = lambda wildcards: [
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{ps}.breadth_ratio.txt"
                for sample, ps in pathogen_pairs_from_checkpoint()
            ],
            entropy_results = lambda wildcards: [
                f"results/pathogen/{sample}/pathogen_mapping/{sample}_{ps}.mean_entropy.txt"
                for sample, ps in pathogen_pairs_from_checkpoint()
            ],
            # Comparison results (for k-mer ranking)
            comparison_results = expand("results/pathogen/{sample}/comparison/{sample}_comparison.tsv", sample=BIO_SAMPLES)
        output:
            scores_matrix = "results/metagenomics/pathogen_detection/detection_scores_matrix.csv",
            scores_heatmap = "results/metagenomics/pathogen_detection/detection_scores_heatmap.pdf",
            detailed_scores = "results/metagenomics/pathogen_detection/detailed_scores.csv"
        log:
            "logs/pathogen_detection/calculation.log"
        conda:
            "workflow/envs/python.yaml"
        script:
            "scripts/calculate_detection_scores.py"

# Pipeline execution report and workflow diagram (runs with or without HOPS)
rule generate_pipeline_report:
    input:
        adapter_removal=expand("results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz", sample=SAMPLES),
        kraken_reports=expand("results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt", sample=BIO_SAMPLES),
        evalue_files=expand("results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv", sample=BIO_SAMPLES),
        pathogen_complete="results/workflow/pathogen_mapping_complete.txt",
        *(
            expand("results/libraries/{sample}/host_mapping/{sample}.dedup.bam", sample=SAMPLES)
            + expand("results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam", sample=SAMPLES)
            if HOST_MTDNA_ANALYSIS_ENABLED
            else []
        ),
        *(
            ["results/metagenomics/pathogen_detection/detection_scores_matrix.csv"]
            if ENABLE_HOPS
            else []
        ),
    output:
        execution_report="results/final/pipeline_execution_report.html",
        timing_data="results/final/pipeline_timing_data.csv",
        workflow_html="results/final/pipeline_workflow_diagram.html",
        workflow_png="results/final/pipeline_workflow_diagram.png",
        workflow_svg="results/final/pipeline_workflow_diagram.svg",
    log:
        "logs/pipeline_report_generation.log",
    conda:
        "workflow/envs/pipeline_report.yaml",
    script:
        "scripts/generate_pipeline_report.py"

if ENABLE_HOPS:
    # Interactive Multi-Method Dashboard
    rule create_interactive_dashboard:
        input:
            scores_matrix = "results/metagenomics/pathogen_detection/detection_scores_matrix.csv",
            detailed_scores = "results/metagenomics/pathogen_detection/detailed_scores.csv",
            comparison_files = expand("results/pathogen/{sample}/comparison/{sample}_comparison.tsv", sample=BIO_SAMPLES),
            hops_results = "results/metagenomics/hops/maltExtract/heatmap_overview_Wevid.tsv",
            abundance_matrix = "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv"
        output:
            dashboard = "results/pathogen_detection_dashboard.html"
        log:
            "logs/interactive_dashboard.log"
        conda:
            "workflow/envs/pipeline_report.yaml"
        shell:
            """
            echo "Creating interactive pathogen detection dashboard at $(date)" > {log}
            python scripts/create_interactive_dashboard.py >> {log} 2>&1
            echo "Interactive dashboard creation completed at $(date)" >> {log}
            """

    # Publication-Ready Output Suite
    rule generate_publication_figures:
        input:
            scores_matrix = "results/metagenomics/pathogen_detection/detection_scores_matrix.csv",
            detailed_scores = "results/metagenomics/pathogen_detection/detailed_scores.csv",
            comparison_files = expand("results/pathogen/{sample}/comparison/{sample}_comparison.tsv", sample=BIO_SAMPLES),
            abundance_matrix = "results/metagenomics/kraken_abundance/krakenuniq_abundance_matrix_absolute.csv"
        output:
            multi_method = "results/final/publication_figures/multi_method_comparison.png",
            method_performance = "results/final/publication_figures/method_performance_analysis.png",
            sample_comparison = "results/final/publication_figures/sample_comparison.png",
            summary_stats = "results/final/publication_figures/summary_statistics.png",
            figure_index = "results/final/publication_figures/README.md"
        log:
            "logs/publication_figures.log"
        conda:
            "workflow/envs/pipeline_report.yaml"
        shell:
            """
            echo "Generating publication-ready figures at $(date)" > {log}
            python scripts/generate_publication_figures.py >> {log} 2>&1
            echo "Publication figures generation completed at $(date)" >> {log}
            """

# Index Building Rules - Simplified approach
# Build BWA index for any reference
rule build_bwa_index:
    input:
        fasta = "{fasta_path}"
    output:
        touch("results/workflow/index_status/bwa_index_{fasta_name}.built")
    log:
        "logs/index_building/bwa_{fasta_name}.log"
    conda:
        "workflow/envs/bwa.yaml"
    shell:
        """
        echo "Building BWA index for {wildcards.fasta_path} at $(date)" > {log}
        bwa index {input.fasta} >> {log} 2>&1
        echo "BWA index building completed for {wildcards.fasta_path} at $(date)" >> {log}
        """

# Build Bowtie2 index for any reference
rule build_bowtie2_index:
    input:
        fasta = "{fasta_path}"
    output:
        touch("results/workflow/index_status/bowtie2_index_{fasta_name}.built")
    log:
        "logs/index_building/bowtie2_{fasta_name}.log"
    conda:
        "workflow/envs/bowtie2.yaml"
    shell:
        """
        echo "Building Bowtie2 index for {wildcards.fasta_path} at $(date)" > {log}
        base_name=$(echo {wildcards.fasta_path} | sed 's/\\.fa$//; s/\\.fasta$//; s/\\.fna$//')
        bowtie2-build {input.fasta} $base_name >> {log} 2>&1
        echo "Bowtie2 index building completed for {wildcards.fasta_path} at $(date)" >> {log}
        """

# Index validation and building wrapper
rule check_and_build_indices:
    output:
        pathogen_indices = "results/workflow/index_status/pathogen_indices_built.txt",
        host_indices = "results/workflow/index_status/host_indices_built.txt",
        mtdna_indices = "results/workflow/index_status/mtdna_indices_built.txt"
    log:
        "logs/index_building/index_check.log"
    script:
        "scripts/check_indices.py"

# ==================== Cleanup Rules for Intermediate Files ====================

if CLEANUP_ENABLED:
    # Cleanup prinseq uncompressed intermediate files
    rule cleanup_prinseq_intermediates:
        input:
            gz = "results/libraries/{sample}/prinseq/{sample}-passed.fq.gz"
        output:
            marker = "results/libraries/{sample}/prinseq/.cleanup_done"
        shell:
            """
            # Remove uncompressed prinseq .fq file if it exists
            rm -f results/libraries/{wildcards.sample}/prinseq/{wildcards.sample}-passed.fq
            touch {output.marker}
            """

    # Cleanup BWA .sai files (temporary alignment index files)
    rule cleanup_bwa_sai_host:
        input:
            bam = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam"
        output:
            marker = "results/libraries/{sample}/host_mapping/.sai_cleanup_done"
        shell:
            """
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_F4.sai
            touch {output.marker}
            """

    rule cleanup_bwa_sai_mtdna:
        input:
            bam = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam"
        output:
            marker = "results/libraries/{sample}/mtdna_mapping/.sai_cleanup_done"
        shell:
            """
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_F4.sai
            touch {output.marker}
            """

    rule cleanup_bwa_sai_pathogen:
        input:
            bam = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
        output:
            marker = "results/pathogen/{sample}/pathogen_mapping/.{ref_name_safe}_sai_cleanup_done"
        shell:
            """
            rm -f results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{wildcards.ref_name_safe}.sai
            touch {output.marker}
            """

    # Cleanup intermediate host BAMs (before dedup)
    rule cleanup_intermediate_host_bams:
        input:
            dedup = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam",
            dedup_bai = "results/libraries/{sample}/host_mapping/{sample}.dedup.bam.bai",
            q30_metrics = "results/libraries/{sample}/host_mapping/{sample}.q30_metrics.txt"
        output:
            marker = "results/libraries/{sample}/host_mapping/.intermediate_cleanup_done"
        shell:
            """
            # Remove intermediate BAMs that are no longer needed after dedup
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_F4.bam
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_F4_q30.bam
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_F4_q30.sorted.bam
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_F4_q30.bam.bai
            # Remove empty host-unmapped files if in default mode (they're just placeholders)
            if [ ! -s results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_host_unaligned.fastq.gz ]; then
                rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_host_unaligned.bam
                rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_host_unaligned.fastq.gz
            fi
            # Remove Bowtie2 log files
            rm -f results/libraries/{wildcards.sample}/host_mapping/{wildcards.sample}_bowtie2.log
            touch {output.marker}
            """

    # Cleanup intermediate mtDNA BAMs (before dedup)
    rule cleanup_intermediate_mtdna_bams:
        input:
            dedup = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam",
            dedup_bai = "results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam.bai",
            q30_metrics = "results/libraries/{sample}/mtdna_mapping/{sample}.q30_metrics.txt"
        output:
            marker = "results/libraries/{sample}/mtdna_mapping/.intermediate_cleanup_done"
        shell:
            """
            # Remove intermediate BAMs that are no longer needed after dedup
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_F4.bam
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_F4_q30.bam
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_F4_q30.sorted.bam
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_F4_q30.bam.bai
            # Remove Bowtie2 log files
            rm -f results/libraries/{wildcards.sample}/mtdna_mapping/{wildcards.sample}_bowtie2.log
            touch {output.marker}
            """

    # Cleanup intermediate pathogen BAMs (before dedup)
    rule cleanup_intermediate_pathogen_bams:
        input:
            dedup = "results/pathogen/{sample}/pathogen_mapping/{sample}_{ref_name_safe}.dedup.bam"
        output:
            marker = "results/pathogen/{sample}/pathogen_mapping/.{ref_name_safe}_intermediate_cleanup_done"
        shell:
            """
            # Remove intermediate pathogen BAMs that are no longer needed after dedup
            rm -f results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{wildcards.ref_name_safe}_F4.bam
            rm -f results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{wildcards.ref_name_safe}_F4_q30.bam
            rm -f results/pathogen/{wildcards.sample}/pathogen_mapping/{wildcards.sample}_{wildcards.ref_name_safe}_F4_q30.sorted.bam
            touch {output.marker}
            """

    # Add cleanup markers to rule all if cleanup is enabled
    # Note: These are added conditionally in the rule all section
