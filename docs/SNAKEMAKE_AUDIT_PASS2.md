# Snakemake audit — Pass 2 (scripts + naming)

**Target:** `PIGSTI_publication/`  
**Date:** 2026-05-21  
**Focus:** Reproducibility, coherence, **output naming**, **scripts/** cross-check vs `Snakefile`

Companion docs: [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md), [`SCRIPTS_INVENTORY.md`](../scripts/SCRIPTS_INVENTORY.md), [Pass 1](SNAKEMAKE_AUDIT_REPORT.md) (superseded in places by this pass).

---

## Naming convention matrix (authoritative)

| Layer | Convention | Example |
|-------|------------|---------|
| Snakefile `safe_name()` | space `/` `\` `:` → `_`; collapse `__` | `Yersinia pestis` → `Yersinia_pestis` |
| `scripts/pigsti_naming.py` | **Same** (use everywhere) | `safe_pathogen_name()` |
| E-value CSV column | `taxonomy` (not filename) | matches `Krakenuniq name` in spreadsheet |
| Kraken report file | `{bio}_kraken-report.txt` | not `kraken-report.txt` alone |
| Pathogen BAM | `results/{bio}/pathogen_mapping/{bio}_{safe}.dedup.bam` | dir name legacy `bwa_*` |
| HOPS heatmap column | `{bio}_unaligned.rma6` | from merged FASTQ |
| HOPS row labels | spaces → underscores in comparisons | `hops_species_token()` |

### Pass-2 naming issues found & fixed

| Issue | Where | Fix in `PIGSTI_publication` |
|-------|-------|------------------------------|
| `EntropyProfile` used `{pathogen}` wildcard | Snakefile | → `{ref_name_safe}` |
| `.entropy_100bp.txt` / `.entropy_1000bp.txt` not Snakemake outputs | Snakefile + entropy script | Added as rule outputs + `expand_downstream_targets` |
| `.ani_distribution.txt` side effect only | `calculate_edit_distance_r2.py` | Declared as `EditDistanceR2` output |
| `safe_name` only replaced spaces in summaries | `summarize_pathogen_data.py`, PCR script | Import `pigsti_naming` |
| HOPS threshold `>=2` in compare script vs `>1` in Snakefile | `compare_pathogens.py` | Aligned to `>1` for comparison TSV |
| Scoring still uses `>=2` for “HOPS detected” point | `calculate_detection_scores.py` | **Intentional** (stricter score than mapping trigger) — document in Methods |

---

## Scripts audit summary

### Path alignment: scripts ↔ Snakefile

| Script expects | Snakefile produces | Status |
|----------------|-------------------|--------|
| `{sample}_kraken-report.txt` | ✓ `krakenuniq` rule | OK |
| `Escore/pathogen/{sample}_pathogen.csv` | ✓ `escore` rule | OK |
| `{sample}_{safe}.dedup.bam` | ✓ checkpoint mapping | OK (use `safe_pathogen_name`) |
| `{safe}.entropy_100bp.txt` | ✓ after pass-2 | **Fixed** |
| `{safe}.ani_distribution.txt` | ✓ after pass-2 | **Fixed** |
| `heatmap_overview_Wevid.tsv` | ✓ HOPS rule | OK (upstream name) |
| `kraken-report.txt` glob (R) | matches `*_kraken-report.txt` | OK |
| `group_by_genus.py` | N/A | **Removed** orphan rule (script never existed) |

### Scripts with hardcoded `results/` paths (acceptable if CWD = repo root)

Most summary/R scripts assume execution from pipeline root — matches Snakemake. **Risk** if users run scripts manually from another directory.

Notable hardcoders: `generate_pipeline_report.py`, `compare_pathogens.py`, `create_interactive_dashboard.py`, `generate_pathogen_report.R` (uses `results_dir` arg from Rscript).

### Scripts with fallback / legacy paths (robust but complex)

`summarize_pathogen_data.py` tries multiple BAM locations (merged, direct, legacy `_F4_q30_sort.bam`) — good for migration, harder to audit. Prefer single canonical path from `pigsti_paths.pathogen_bam()`.

### Non-deterministic / dev-only scripts

| Script | Note |
|--------|------|
| `evaluate_escore_functions.py` | `np.random` — not in DAG |
| `generate_pipeline_report.py` | `datetime.now()` in HTML — monitoring only |
| `create_multi_qc_dashboard.py` | timestamp in HTML |

---

## Checklist (pass 2 highlights)

### Reproducibility
- [x] Manifest + sorted checkpoint (pass 1)
- [ ] Tool seeds (BWA/Bowtie2) not pinned — **BAM bytes may differ across CPUs**
- [x] Canonical naming module `pigsti_naming.py`
- [ ] Conda envs still floating versions

### Pipeline coherence
- [x] 51 scripts catalogued — see [`SCRIPTS_INVENTORY.md`](../scripts/SCRIPTS_INVENTORY.md)
- [x] Orphan rules removed; Diamond still README-only
- [x] HOPS/compare/mapping threshold documented (`>1` vs scoring `>=2`)

### Output naming
- [x] `OUTPUT_SCHEMA.md` + entropy/ANI files now first-class outputs
- [ ] Column names in Excel mix styles (`Relative_entropy_100bp`, `Entropy100` in `merge_pathogen_plots.py`)

### Bugs & edge cases
- [x] Empty samples.tsv → DAG error
- [x] validate_pigsti_setup
- [ ] `summarize_host_mtdna` references `species_mismatch_warning.txt` — verify rule produces it
- [ ] `compare_pathogens` writes stub HTML on missing HOPS column (fail-open)

### Efficiency
- [ ] `summarize_pathogen_data.py` opens BAM twice per pathogen
- [ ] HOPS single job all samples — RAM wall

---

## Output format (required)

### Critical issues

**C1 — Inconsistent pathogen sanitization across scripts**  
- **Location:** `summarize_pathogen_data.py` (was `.replace(" ","_")` only), `Snakefile` `safe_name()`, R reports  
- **Why critical:** Pathogens with `/` or `:` in spreadsheet could map to different files in Snakefile vs summaries  
- **Fix:** `pigsti_naming.safe_pathogen_name()` shared module; updated key Python scripts  

**C2 — Untracked side-output files (entropy 100/1000bp, ANI distribution)**  
- **Location:** `calculate_entropy_profile.py`, `calculate_edit_distance_r2.py` vs `expand_downstream_targets`  
- **Why critical:** Snakemake could mark pathogen step complete before side files exist → summaries show `NA`  
- **Fix:** Added explicit outputs on rules + checkpoint barrier list  

**C3 — HOPS comparison threshold mismatch**  
- **Location:** `compare_pathogens.py` used `>=2`, Snakefile/summaries use `>1`  
- **Why critical:** Comparison TSV disagrees with mapping decisions  
- **Fix:** `compare_pathogens.py` aligned to `>1` (scoring script still uses `>=2` by design)  

### Warnings

**W1 — `{sample}` wildcard means PCR or bio-sample**  
- Scripts using `results/{sample}/` must match rule wildcards — document in every methods paragraph  

**W2 — Legacy directory `bwa_pathogen` / `host_mapping`**  
- All aligner scripts write under `bwa_*` paths  

**W3 — 15+ scripts not in DAG**  
- Diamond, `merge_pathogen_plots`, ENA tools — see inventory  

**W4 — R/Python column name drift**  
- E-score: `taxonomy`, `taxReads`, `# of reads`, `reads` — handled with fallbacks in R/Python  
- Risk when adding new columns  

**W5 — `calculate_detection_scores` glob discovery**  
- Uses `glob("results/*/evalue/...")` — can pick non-bio sample dirs if present  

**W6 — `krakenuniq_abundance_matrix.R` pattern**  
- `kraken-report.txt$` matches `{id}_kraken-report.txt` — OK  

### Suggestions

**S1** — Migrate all scripts to `from pigsti_paths import ...`  
**S2** — Add `scripts/test_naming.py` unit tests for `safe_pathogen_name` edge cases  
**S3** — Pin conda envs at release tag  
**S4** — Move standalone utilities to `scripts/cli/`  

### Positive findings

- Rich fallback logic in `summarize_pathogen_data.py` and `generate_pathogen_report.R`  
- `dExp_Escore.py` reads thresholds from config  
- `compare_pathogens.py` fail-open stubs prevent DAG halt on missing HOPS column  
- `validate_pigsti_setup.py` + manifest for publication  
- Spreadsheet-driven pathogen list with checkpoint union (Kraken + HOPS)  

### Summary & recommendations

**Top 3 priorities**
1. Use **`PIGSTI_publication`** for release (naming + tracked outputs + validation).  
2. Methods text: define **`>1` HOPS** for mapping vs **`>=2`** for scoring.  
3. Add **CI**: `python -m pytest scripts/test_naming.py` + `snakemake -n`.

**Overall:** Reproducibility **B+**; naming **B→A-** after pass-2 fixes; scripts **B** (many utilities, good fallbacks, some drift).  

**Effort:** Pass-2 code fixes ~2h; full script refactor to `pigsti_paths` ~1 day; modular Snakefile ~2–4 days.

---

## Changes applied in pass 2 (`PIGSTI_publication` only)

- `scripts/pigsti_naming.py`, `scripts/pigsti_paths.py`
- `scripts/SCRIPTS_INVENTORY.md`
- Updated: `summarize_pathogen_data.py`, `summarize_pathogen_data_pcr.py`, `calculate_detection_scores.py`, `create_multi_qc_dashboard.py`, `compare_pathogens.py`
- Snakefile: entropy + ANI distribution outputs in barrier list
- `calculate_entropy_profile.py`, `calculate_edit_distance_r2.py` wired to new outputs
