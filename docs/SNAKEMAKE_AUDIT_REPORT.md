# Snakemake pipeline audit — PIGSTI_publication

**Audited copy:** `PIGSTI_publication/` (fork of `PIGSTI/` for publication fixes)  
**Date:** 2026-05-21  
**Scope:** ~80 rules, 26 conda envs, ancient DNA pathogen + host/mtDNA workflow  

**Pipeline purpose:** Metagenomic pathogen screening and host/mtDNA characterization from degraded DNA (multi-PCR per biological sample).  
**Inputs:** Paired/single FASTQ per PCR (`config/samples.tsv`), pathogen reference spreadsheet.  
**Infrastructure:** Conda via `--use-conda`; designed for Linux HPC (paths in config).  

---

## Reproducibility checklist (filled)

### Random seeds & stochasticity
- [x] **Partial** — Production rules use deterministic aligners/filters; `evaluate_escore_functions.py` uses `np.random` (dev-only, not in `rule all`).
- [ ] **Gap** — No global `seed:` in config; Bowtie2/BWA thread order can cause tiny BAM differences across CPUs.
- [x] **Mitigation** — Sorted pathogen lists, `sort_keys` in manifest JSON, deterministic workflow diagram.

### Input dependencies
- [x] Most inputs declared in rules; checkpoint reads `pathogen_targets.txt`.
- [ ] **Gap** — `get_sample_ref_pairs()` reads HOPS TSV via `os.path.exists` at parse time (early DAG).
- [x] **Fixed in publication copy** — `validate_pigsti_setup` checks FASTQ paths before heavy steps.
- [x] External: ENA download scripts (optional, outside Snakemake DAG).

### Conda/container environments
- [x] ~83 `conda:` directives for ~80 rules; no bare shell-only rules found.
- [ ] **Gap** — Env YAMLs use `python=3.11` / floating pins, not full lockfiles.
- [ ] No Singularity/Docker profiles (conda-only).

### Workflow parameters
- [x] Central `config/config.yaml` + optional `config_overlay` / `PIGSTI_CONFIG_OVERLAY`.
- [ ] **Gap** — Committed `config.yaml` contains `/raid_md0/...` paths (machine-specific).
- [x] **Fixed in publication copy** — `config/config.example.yaml` + `.gitignore` for local config.

### Logging & provenance
- [x] Per-rule logs under `logs/`; `results/run_manifest.json` (SHA256, git, effective config).
- [x] `results/workflow/pathogen_targets.txt` checkpoint audit trail.
- [ ] Timestamps in HTML reports (non-scientific fields).

---

## Pipeline coherence checklist (filled)

### Rule dependencies
- [x] Checkpoint `generate_pathogen_targets` correctly gates pathogen mapping.
- [ ] **Risk** — Memoized `_pathogen_checkpoint_pairs_memo` must not stale across checkpoint re-runs in same Python process (rare Snakemake edge case).
- [x] No circular dependencies identified.
- [ ] Implicit: host reference maps built from files read at runtime in shell/python.

### Workflow logic
- [x] Clear PCR → sample merge → detection → mapping → summary flow.
- [ ] **Orphan removed in publication copy** — `select_best_pathogens` (was unreachable).
- [ ] **Legacy removed** — `get_bwa_targets()` → obsolete `bwa_final/` paths.
- [x] Conditional blocks documented: `ENABLE_HOPS`, `PATHOGEN_ALIGNER`, `HOST_MTDNA_ANALYSIS_ENABLED`, `PATHOGEN_MODE`, `CLEANUP_ENABLED`.

### Naming conventions
- [x] Rule names mostly descriptive (`adapter_removal_pe`, `merge_unaligned_fastq_per_sample`).
- [ ] **Confusion** — `{sample}` wildcard = PCR ID or bio sample depending on rule (see `OUTPUT_SCHEMA.md`).
- [ ] **Legacy dirs** — `host_mapping` when aligner is Bowtie2.

### Modularity
- [ ] **Gap** — Single 3.7k-line `Snakefile` (maintainability risk).
- [ ] Repeated `*_by_type` symlink rules (could be helper function or Snakemake 8 module).

---

## Output folder naming checklist (filled)

### Directory structure
- [x] `results/{pcr}/` vs `results/{bio_sample}/` vs `results/sample_*` documented.
- [x] `results/final/` for publication tables.
- [x] `results/browse/` mirror via symlinks.
- [ ] Depth is deep (4–5 levels) but consistent.

### Naming patterns
- [x] `{sample}_{artifact}.{ext}` pattern dominant.
- [x] `pathogen_safe` sanitization unified.
- [ ] No tool versions in filenames (acceptable if manifest records versions).

### Path conflicts
- [x] Checkpoint prevents duplicate pathogen BAM targets when sorted/deduped.
- [ ] **Risk** — `rm -rf {output}` in `decom_run` before run (destructive but intentional).
- [x] Cleanup rules only run after dedup BAMs exist.

### Documentation
- [x] `docs/OUTPUT_SCHEMA.md`, `README.md`, workflow PNG.
- [x] `run_manifest.json` as machine-readable index.

---

## Bugs & edge cases checklist (filled)

### Input handling
- [x] `strict_inputs` + `_resolve_existing_fastq_path` for FASTQs.
- [x] **Fixed** — Empty `samples.tsv` now fails at DAG build.
- [x] **Fixed** — `validate_pigsti_setup` for missing files.
- [ ] Empty BAMs: many rules stub outputs (QualiMap, DamageProfiler) — good for robustness, document in Methods.

### Wildcard handling
- [x] `global_wildcard_constraints` on `sample`, `ref_name_safe`.
- [ ] Pathogens with special chars rely on `safe_name()` — tested for space/slash/colon.
- [x] **Fixed** — `EntropyProfile` used `{pathogen}` vs `{ref_name_safe}` (DAG breakage risk).

### Error handling
- [x] Many shell blocks use `set -euo pipefail`.
- [ ] **Was critical** — `decom_run ... || true` hid failures → **fixed** in publication copy.
- [ ] Some rules use `|| true` on qualimap (intentional stub path).

### Data loss risks
- [x] `cleanup_intermediates` deletes intermediates only after final BAMs (marker-gated).
- [x] Raw input FASTQs are **not** deleted by the pipeline (ENA cleanup rule removed).
- [ ] Checkpoint re-run requires understanding Snakemake checkpoint semantics.

### Tool-specific
- [ ] README documents Diamond but **no Snakemake rules** (documentation bug).
- [x] HOPS optional; pipeline completes without it.

### Sample/batch
- [x] **Fixed** — zero-sample guard.
- [x] Single-sample: should work (untested in CI).
- [ ] HOPS rule passes all sample FASTQs in one job (memory bottleneck).

---

## Efficiency checklist (filled)

### Computation
- [x] Most align/classify rules set `threads`.
- [ ] **Was issue** — `summarize_host_mtdna` `threads: 999` → **fixed** `summary_threads` config.
- [ ] Pathogen mapping one job per (sample × pathogen) — correct but can explode job count.
- [x] Checkpoint avoids mapping all spreadsheet pathogens.

### I/O
- [ ] `cat` merge of all PCR FASTQs per sample (disk-heavy for large cohorts).
- [x] `cleanup_intermediates` reduces disk after completion.
- [ ] Symlinks in `results_by_type` — good for UX, can confuse mtime on NFS.

### Workflow efficiency
- [ ] **Bottleneck** — `rule hops` waits for all samples + runs monolithic HOPS (800G RAM in shell).
- [ ] `check_and_build_indices` serializes index checks at end of `rule all` (acceptable).

### Caching
- [x] Snakemake timestamp-based rerun works.
- [x] `--rerun-incomplete` documented.

---

## Output format

### Critical issues (reproducibility / data loss)

| # | Issue | Location | Why critical | Fix (publication copy) |
|---|--------|----------|--------------|------------------------|
| C1 | `decom_run \|\| true` | `Snakefile` `decom_run` | Silent failure → false-positive “completed” decOM | Removed `\|\| true`; `decom_fail_on_error` |
| C2 | `EntropyProfile` wrong wildcard | `Snakefile` ~1341 | Rule may never match checkpoint outputs | Renamed to `{ref_name_safe}` |
| C3 | Duplicate config load overwrote overlay | `Snakefile` ~2246 (original) | Wrong reference maps mid-run | Removed second load (in both trees) |
| C4 | Machine paths in git config | `config/config.yaml` | Others cannot reproduce | `config.example.yaml` + `.gitignore` |
| C5 | ~~ENA raw FASTQ deletion~~ | removed | Was irreversible if `source=ENA` mis-set | Use external disk management if needed |

### Warnings (coherence / efficiency / maintainability)

| # | Issue | Impact | Suggestion |
|---|--------|--------|------------|
| W1 | 3700-line monolithic Snakefile | Hard review/cite | Split into `workflow/rules/*.smk` |
| W2 | Diamond in README, not in workflow | User confusion | Implement or remove from docs |
| W3 | `bwa_*` dir names with Bowtie2 aligner | Reviewer confusion | Rename in v2 + symlinks |
| W4 | `enable_hops` default True in code, False in config | Ambiguous for new users | Default False in publication copy |
| W5 | HOPS monolithic rule | RAM/queue bottleneck | Per-sample HOPS or external workflow |
| W6 | Unpinned conda envs | Version drift | Ship lock files per release |
| W7 | `get_sample_ref_pairs` parse-time file reads | Fragile DAG | Checkpoint-only pathogen selection only |
| W8 | Monitoring was HOPS-gated (original) | Missing reports without HOPS | Always run `generate_pipeline_report` (fixed earlier) |

### Suggestions (optimization / best practices)

| # | Suggestion | Benefit |
|---|------------|---------|
| S1 | CI: `snakemake -n` on 2-sample fixture | Catch DAG breaks |
| S2 | `profiles/hpc.yaml` with `default-resources` | Cluster portability |
| S3 | Combine small QC rules per pathogen where I/O-bound | Fewer job submissions |
| S4 | Schema validate `samples.tsv` with Pandera | Early typo detection |
| S5 | Tag releases + Zenodo DOI linked in manifest | Publication standard |

### Positive findings

- Checkpoint-driven pathogen targeting (scientifically sound + efficient).
- Rich QC layer (ANI, breadth, entropy, edit-distance damage split, PDF reports).
- Conda isolation per step.
- `run_manifest.json` + sorted checkpoint outputs.
- PCR/sample merge model matches lab reality.
- Fail-open stubs for empty alignments avoid pipeline stalls.
- `config_overlay` for user-friendly deployments.

### Summary & recommendations

**Top 3 priority fixes**
1. **Publication config** — Use `PIGSTI_publication` with `config.example.yaml`, never commit private paths.  
2. **Modularize Snakefile** before Zenodo release (no science change).  
3. **Resolve README/tool drift** (Diamond, Krona) and add minimal CI dry-run.

**Overall assessment**

| Dimension | Grade | Notes |
|-----------|-------|-------|
| Reproducibility | **B+** | Strong manifest + conda; weak env pinning + BAM byte identity |
| Coherence | **B** | Logical flow; wildcard naming needs discipline |
| Output naming | **B+** | Documented; legacy `bwa_*` names |
| Bugs / edge cases | **B** | After publication-copy fixes; HOPS RAM still risky |
| Efficiency | **B-** | Good parallelism; HOPS + merge CAT are bottlenecks |

**Estimated effort**

| Task | Effort |
|------|--------|
| Use publication copy + example config | 1 hour |
| Full Snakefile modularization | 2–4 days |
| Conda lock + test dataset + CI | 2–3 days |
| Folder rename v2 | 1–2 days |

---

## Changes applied in `PIGSTI_publication/` only

- `validate_pigsti_setup` rule + script  
- `config/config.example.yaml`, `.gitignore`  
- Removed orphan `select_best_pathogens`, legacy `get_bwa_targets`  
- `EntropyProfile` wildcard fix  
- decOM strict failure (`decom_fail_on_error`)  
- `summary_threads`, empty-sample guards, `enable_hops` default `false`  
- `PUBLICATION_README.md`  

Original pipeline in `../PIGSTI/` is unchanged unless you merge these commits.
