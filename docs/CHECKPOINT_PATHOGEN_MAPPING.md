# Pathogen mapping checkpoint — design and `--rerun-incomplete` fix

## Problem

`checkpoint generate_pathogen_targets` decided which BAMs to map by calling `get_sample_ref_pairs()`, which used **`os.path.exists()`** per sample. During `snakemake --rerun-incomplete`:

1. Some samples still had **incomplete** `escore` jobs, but **old CSVs** remained on disk.
2. The checkpoint could run when **only a subset** of `BIO_SAMPLES` had readable files, writing a partial `results/workflow/pathogen_targets.txt`.
3. Pathogen mapping then started for that subset while other samples were still running Kraken/E-score.

Snakemake’s declared checkpoint inputs (`expand(..._pathogen.csv)`) were correct in theory, but the **run block ignored those inputs** and re-scanned disk, weakening the guarantee.

## Solution (three layers)

### 1. `rule metagenomics_screening_cohort_ready`

**Inputs (all must be finished jobs):**

- Every `results/{bio}/evalue/pathogen|genus|species/*.csv`
- Every `results/{bio}/krakenuniq/*_kraken-report.txt`
- HOPS heatmap (if `enable_hops`)

**Output:** `results/workflow/metagenomics_screening_cohort_ready.json`

Validates: non-empty pathogen CSVs, `taxonomy` column, **sample set == `BIO_SAMPLES`**.

### 2. `checkpoint generate_pathogen_targets`

**Additional input:** `metagenomics_screening_cohort_ready.json`

**Run:** `scripts/build_pathogen_target_pairs.py` with:

- `--evalue` (alias `--escore`) = Snakemake `input.evalue_files` only (no disk glob)
- `--strict-cohort` → fails if any biological sample missing
- Writes `pathogen_targets.txt` + **`pathogen_targets.manifest.json`**

### 3. `rule pathogen_mapping_targets`

Requires manifest + all QC outputs from `expand_downstream_targets`. Run block checks `cohort_samples` in manifest matches `BIO_SAMPLES`.

## Operational notes

| Situation | Action |
|-----------|--------|
| Stuck after partial rerun | Remove `results/workflow/pathogen_targets.txt`, `pathogen_targets.manifest.json`, `pathogen_mapping_complete.txt`, then `snakemake --rerun-incomplete` |
| E-score changed for one sample | Cohort JSON + checkpoint invalidate → targets regenerated |
| Zero pathogens in cohort | Allowed; manifest `n_mapping_pairs: 0`; warning in barrier rule |

## Files

| File | Role |
|------|------|
| `results/workflow/metagenomics_screening_cohort_ready.json` | Cohort gate |
| `results/workflow/pathogen_targets.txt` | One BAM path per line (checkpoint output); v3 paths use `results/pathogen/{bio}/pathogen_mapping/...` |
| `results/workflow/pathogen_targets.manifest.json` | Cohort metadata + pair list |
| `results/workflow/pathogen_mapping_complete.txt` | Barrier before summaries |
| `scripts/build_pathogen_target_pairs.py` | Pair logic (shared with Snakemake) |

## PCR `super_careful` summaries

PCR-level summarize rules now use `pathogen_pairs_from_checkpoint()` (not `get_sample_ref_pairs()` disk scan) so they stay aligned with the checkpoint target list.
