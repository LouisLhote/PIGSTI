# PIGSTI pipeline audit (publication readiness)

Audit date: 2026-05-21. Scope: `Snakefile`, `config/`, `scripts/`, `workflow/envs/`.

**Consolidated overall audit (latest):** [`OVERALL_AUDIT_2026-05-21.md`](OVERALL_AUDIT_2026-05-21.md).  
**Outputs/summaries:** [`OUTPUTS_AND_SUMMARIES_AUDIT.md`](OUTPUTS_AND_SUMMARIES_AUDIT.md).

## Executive summary

PIGSTI is **feature-rich and scientifically thorough**, but **publication friction** comes from:

1. A **3700+ line monolithic Snakefile** with duplicated config loading and conditional rule variants.
2. **Documentation drift** (README features not wired in the workflow).
3. **Monitoring / reporting gated on HOPS** even when HOPS is off (common default in `config.yaml`).
4. **Legacy directory names** (`bwa_*`) that confuse reviewers when Bowtie2 is used.

None of these block science, but they should be addressed before a Zenodo/GitHub release.

---

## Workflow map

See [`docs/images/pipeline_workflow.png`](images/pipeline_workflow.png) (regenerate: `python scripts/render_workflow_diagram.py`).

Rough scale: **~70 named rules**, **15+ config toggles**, **2 wildcard namespaces** (PCR vs biological sample).

---

## Config toggles (what the diagram should show)

| Key | Default in code | Your `config.yaml` | Effect |
|-----|-----------------|-------------------|--------|
| `pathogen_mapping_mode` | `super_careful` | `default` | Pathogen reads: merged Bowtie2-unmapped vs BWA host-unaligned |
| `pathogen_aligner` | `bwa` | `bowtie2` | BWA aln/samse vs Bowtie2 end-to-end for pathogens |
| `host_aligner` | `bwa` | `bowtie2` | Host/mtDNA via `bwa_aln_*.py` scripts (internal branch) |
| `enable_hops` | `True` | `false` | HOPS, compare_pathogens, detection scores, pipeline HTML report |
| `enable_decom` | `True` | `true` | decOM (note: `decom_run` uses `\|\| true`) |
| `pathogen_screening_only` | `false` | `false` | Skips host/mtDNA arm |
| `cleanup_intermediates` | `true` | `true` | Deletes intermediate BAMs/SAI |
| `dedup_tool` / `merge_tool` | `picard` | `samtools` | Markdup / merge implementation |
| `edit_distance_damage_split` | `true` | `true` | Extra damage vs no-damage BAMs + R² |
| `enable_multi_qc_dashboard` | `false` | `true` | Static HTML QC dashboard |
| `strict_inputs` | `true` | (not set) | Fail-fast on missing inputs |

---

## Strengths (keep)

- **Checkpoint-driven pathogen mapping** avoids mapping every spreadsheet row.
- **Sorted pathogen target lists** and memoized checkpoint parsing (Snakemake#823 workaround).
- **Per-rule conda envs** — good reproducibility primitive.
- **`run_manifest.json`** with file hashes (now includes git + effective config).
- **Dual output layout** (`results/` + `results/browse/` symlinks).
- **Rich pathogen QC**: ANI, breadth, entropy, edit-distance decay, damage split, PDF reports.
- **ENA integration** via `samples.tsv` (`source=ENA` metadata); raw FASTQs are retained.

---

## High-priority improvements

### 1. Split the Snakefile into modules

**Problem:** Single 3.7k-line file is hard to review, test, and cite.

**Change:** Snakemake 8 style layout:

```
workflow/
  Snakefile          # rule all + imports
  rules/preprocess.smk
  rules/host_mtdna.smk
  rules/metagenomics.smk
  rules/pathogen.smk
  rules/reporting.smk
  rules/cleanup.smk
common.smk           # shared functions, CFG, samples
```

### 2. Fix documentation vs code drift

| Documented | In Snakefile? |
|------------|-------------|
| Diamond screening (`enable_diamond`) | **No rules** — README only |
| Krona HTML | **Removed** (comment in Snakefile) |
| `host_aligner: bowtie2` | **Yes** — inside `bwa_aln_host.py` / `bwa_aln_mtdna.py`, not separate rules |
| Pipeline monitoring report | **Only if `enable_hops`** |

**Change:** Either implement Diamond as optional rules in `rule all`, or remove from README until implemented.

### 3. Un-gate monitoring from HOPS

**Problem:** `generate_pipeline_report`, detection scores, and `pipeline_workflow_diagram.*` are in `rule all` only when `ENABLE_HOPS`. Most runs with `enable_hops: false` get no workflow diagram in `results/`.

**Change:** Move `generate_pipeline_report` and diagram outputs outside the `if ENABLE_HOPS` block; keep HOPS-specific inputs optional inside the rule.

### 4. Remove duplicate `config.yaml` load

**Problem:** Lines ~2246–2250 reload `CFG` mid-file after host rules start — risk of inconsistent globals.

**Change:** Delete second load; use top-level `CFG` only.

### 5. Align defaults between Snakefile and `config.yaml`

**Problem:** `ENABLE_HOPS = CFG.get("enable_hops", True)` but config has `enable_hops: false`. New users following code defaults ≠ config file.

**Change:** Document that **config.yaml is authoritative**; consider default `False` in code to match typical faster runs.

### 6. `decom_run` silent failure

```snakemake
decOM ... -o {output} || true
```

**Change:** Remove `|| true` or write a marker file + log warning; add `enable_decom: false` to `rule all` when disabled (already conditional).

### 7. Orphan / dead rules

- `select_best_pathogens` → `scripts/group_by_genus.py` and `sample_pathogen.csv` — **not in `rule all`**; likely dead.
- `get_bwa_targets()` → `bwa_final/q30/` paths — **legacy**, unused by current pathogen flow.

**Change:** Delete or wire into workflow; reduces reviewer confusion.

---

## Medium-priority (quality & naming)

### 8. Rename output folders (major version)

| Current | Clearer |
|---------|---------|
| `host_mapping/` | `host_align/` or `align_host/` |
| `pathogen_mapping/` | `pathogen_align/` |
| `Escore/` | `evalue/` (v2; see [`V2_MIGRATION.md`](V2_MIGRATION.md)) |

Provide symlink aliases for one release cycle.

### 9. `PICARD_CMD` vs hardcoded `picard`

Pathogen `mark_duplicates` uses `picard MarkDuplicates` but host rules use `{PICARD_CMD}`.

**Change:** Use `{PICARD_CMD}` everywhere.

### 10. Wildcard `EntropyProfile` inconsistency

Rule uses `{pathogen}` wildcard; other pathogen rules use `{ref_name_safe}`.

**Change:** Standardize on `ref_name_safe`.

### 11. Pin conda envs for publication

`PIGSTI_snakemake.yaml` pins Snakemake env; **per-rule envs** use floating `python=3.11` etc.

**Change:** For release tag, `conda env export` per env or use `conda-lock` / `--conda-create-envs-only` once and ship lockfiles.

### 12. Machine-specific paths in committed `config.yaml`

**Change:** Ship `config/config.example.yaml` with placeholders; gitignore real `config.yaml` or use `config_overlay` only on HPC.

### 13. `summarize_host_mtdna` threads: 999

Can oversubscribe shared clusters.

**Change:** `threads: min(32, workflow.cores)` or config key `summary_threads`.

---

## Low-priority (nice for v2)

- **Schema validation** for `samples.tsv` and `Pathogen_spreadsheet.csv` at DAG build (pandera/cerberus).
- **CI:** `snakemake -n` on tiny fixture dataset in GitHub Actions.
- **nf-core style profiles:** `profiles/local.yaml`, `profiles/cluster.yaml`.
- **CITATION.cff** + version tag aligned with manifest `manifest_version`.
- **Test dataset** (<500 MB) for reviewers to run end-to-end in <2 h.

---

## Reproducibility checklist (for the paper)

1. Git tag + `results/run_manifest.json`
2. Archived `config.yaml`, `samples.tsv`, `Pathogen_spreadsheet.csv`
3. Kraken + reference genome versions (not only paths)
4. Note `enable_hops`, `pathogen_mapping_mode`, `host_aligner`, `pathogen_aligner`
5. State that BAMs may differ slightly across CPUs but **tables/checkpoint targets** should match

---

## Suggested release order

1. Fix HOPS-gated monitoring + README drift (Diamond).
2. Modularize Snakefile (no science change).
3. `config.example.yaml` + test profile.
4. Tag v1.0.0 + Zenodo DOI.
5. Folder renames in v2.0.0.
