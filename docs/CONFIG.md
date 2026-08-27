# PIGSTI configuration guide

## Quick reference

| Mandatory | Optional |
|-----------|----------|
| `config/samples.tsv` (`sample`, `r1`, optional `pcr`/`r2`) | HOPS (`enable_hops`, malt index, resources) |
| `config/Pathogen_spreadsheet.csv` (`Krakenuniq name`, `Hops name`, `bwa index`) | decOM (`enable_decom`, `decOM_sources`) |
| `kraken_db`, `host_index` | `pathogen_screening_only`, per-pathogen spreadsheet thresholds |
| `bwa_indices` + `mtDNA_indices` if `host_aligner: bwa` | Most thresholds under `pathogen_detection_criteria` |
| `bowtie2_indices` + `bowtie2_mtDNA_indices` if `host_aligner: bowtie2` | `cleanup_intermediates`, thread counts |

Pathogen references for **mapping** come from the spreadsheet `bwa index` column, not from YAML.

Copy the template and edit paths for your machine:

```bash
cp config/config.example.yaml config/config.yaml
```

### Results directory (`results_root`)

By default, outputs go under **`results/`** in the pipeline folder. To write elsewhere (large disk, shared RAID, etc.):

```bash
# Snakemake CLI (recommended for one-off runs)
snakemake -j 8 --config results_root=/raid_md0/louis/PIGSTI_runs/run_2026_05

# Or in config/config.yaml
# results_root: /raid_md0/louis/PIGSTI_runs/run_2026_05

# Or environment variable
export PIGSTI_RESULTS_ROOT=/raid_md0/louis/PIGSTI_runs/run_2026_05
snakemake -j 8
```

At startup the Snakefile creates a symlink **`results` → your `results_root`** so all existing paths (`results/libraries/…`, `results/final/…`) still work without editing hundreds of rules. The resolved path is recorded in `.pigsti/results_root`.

**Verify after launch:** `ls -la results` should show `results -> /your/results_root`. If `results` is a normal directory (`drwx…`) and `PIGSTI_run_output` does not exist, the redirect did not apply — see troubleshooting below.

**Still in the pipeline directory** (not moved): `config/`, `scripts/`, `lists/`, `logs/`, `.snakemake/`. Launch Snakemake from the pipeline root (`PIGSTI_publication/`).

If `./results` already exists as a **non-empty real folder** and you set a different `results_root`, Snakemake will stop with a clear error — rename or archive `./results` first.

**Troubleshooting:** `.pigsti/results_root` can list your external path while `./results` is still a real folder (outputs under `PIGSTI_publication/results/`). That happens if a run started before `results/` had data, or if `PIGSTI_RESULTS_ROOT` was set in the shell but `./results` was never replaced. Fix:

```bash
mv results results_local_backup_$(date +%Y%m%d)
mkdir -p /raid_md0/.../PIGSTI_run_output
rsync -a results_local_backup_*/ /raid_md0/.../PIGSTI_run_output/   # optional
snakemake ... --config results_root=/raid_md0/.../PIGSTI_run_output
ls -la results    # must be: results -> /raid_md0/.../PIGSTI_run_output
```

Validation runs automatically when Snakemake loads the Snakefile (`scripts/validate_pigsti_setup.py`). You can also test manually:

```bash
python scripts/validate_pigsti_setup.py --config config/config.yaml --samples config/samples.tsv
```

---

## Mandatory inputs (three files + config keys)

### Files (must exist before a real run)

| File | Required columns / content |
|------|----------------------------|
| **`config/samples.tsv`** | Tab-separated. **`sample`** (bio ID), **`pcr`** (library ID; optional — defaults to `sample`), **`r1`** (and **`r2`** if paired-end). Paths must exist on disk. Optional: `RGLB`, `sequencing_run`, `source` (`LOCAL` / `ENA` metadata only). |
| **`config/Pathogen_spreadsheet.csv`** | **`Krakenuniq name`**, **`Hops name`**, **`bwa index`** (pathogen reference FASTA or index path). One row per pathogen you may map. |
| **`config/config.yaml`** | See mandatory keys below. |

Snakefile paths default to `config/samples.tsv` and `config/Pathogen_spreadsheet.csv` unless you override in YAML.

### Config keys enforced by validation (errors if missing/wrong)

| Key | What it is |
|-----|------------|
| **`kraken_db`** | Existing **KrakenUniq** database directory (folder built with `krakenuniq-build`). |
| **`host_index`** | Bowtie2 index **prefix** for chimera / unaligned extraction (`bowtie2_unaligned`). Files must exist, e.g. `.../chimera.1.bt2l` → set `host_index: ".../chimera"`. |
| **`bwa_indices`** | Map species → host nuclear FASTA (e.g. `Pig`, `Human`). Required when `host_aligner: bwa` (default). |
| **`mtDNA_indices`** | Map species → mtDNA FASTA. Required when `host_aligner: bwa`. |
| **`bowtie2_indices`** | Map species → Bowtie2 prefix. Required when `host_aligner: bowtie2` (or use FASTA in `bwa_indices` to build at runtime — validator warns). |
| **`bowtie2_mtDNA_indices`** | Same for mtDNA when `host_aligner: bowtie2`. |

Pathogen reference FASTAs for **mapping** come from the spreadsheet `bwa index` column, not from `config.yaml`.

---

## Strongly recommended (defaults exist)

| Key | Default | Notes |
|-----|---------|--------|
| `pathogen_mapping_mode` | `default` *(recommended)* | `default` = pooled unaligned FASTQs; `super_careful` = per-PCR host-unmapped reads (needs host mapping). |
| `pathogen_aligner` | `bwa` | `bwa` or `bowtie2` for pathogen reference mapping. |
| `host_aligner` | `bwa` | `bwa` or `bowtie2` for host/mtDNA. |
| `pathogen_screening_only` | `false` | `true` skips host/mtDNA and forces `pathogen_mapping_mode: default`. |
| `strict_inputs` | `true` | Fail fast on empty/missing critical inputs. |
| `fastq_screen_full_dataset_rescreen` | `true` | If best species `#One_hit_one_genome` on the default subset is below the threshold, re-run FastQ Screen with `--subset 0` (full collapsed FASTQ). |
| `fastq_screen_full_dataset_min_one_hit` | `50` | Minimum `#One_hit_one_genome` reads on the subset pass before triggering a full-dataset re-screen. |
| `enable_sexing` | `true` | Run chromosome-residual sexing after host mapping (requires `pathogen_screening_only: false`). |

Each library writes `results/libraries/{pcr}/fastq_screen/{pcr}_fastq_screen_run_mode.txt` (`pass=subset_default` or `pass=full_dataset`).

### Genetic sexing (Cow, Goat, Sheep, Dog)

After **merged host BAM** per biological sample (`merge_host_dedup_per_sample`), rule `sexing_residual` runs when `enable_sexing: true`. Host species is resolved from FastQ Screen across PCR libraries (skipped if species mismatch).

| FastQ Screen species | Typical autosome count (QC hint) |
|----------------------|----------------------------------|
| Cow, Goat | 29 |
| Sheep | 26 |
| Dog | 38 |

**Scripts:** `scripts/sexing_residual_method.R` (or `scripts/sexing/sexing_residual_method.R`) — declared as a Snakemake input; either layout is accepted.

**Reference handling:** Python reads `samtools idxstats` from the merged host BAM, finds numbered autosomes (`1`…`N` or `chr1`…`chrN`) and **X** by contig name (Y ignored). The R script does **not** hard-code row numbers — it uses all rows except the last as autosomes and the last row as X.

Outputs per **biological sample** (`{bio}`):

- `results/samples/{bio}/sexing/{bio}_sexing.idx`
- `results/samples/{bio}/sexing/{bio}_sexing.pdf`
- `results/samples/{bio}/sexing/{bio}_sexing.tsv`

Host reference chromosomes must be named so autosomes parse as `1`…`N` and X as `X` / `chrX` (standard Ensembl-style references).

**Final outputs**

| Artefact | Path |
|----------|------|
| Sexing plot | `results/samples/{bio}/sexing/{bio}_sexing.pdf` |
| Host/mtDNA table | `results/final/host_mtdna_summary_all_samples.xlsx` — `Sample_level` + `Sexing` sheet |
| Comprehensive table | `results/final/comprehensive_summary_all_samples.xlsx` — `sexing_call`, `sexing_plot`, … on `Sample_level` |

---

## Optional features (off unless you set paths)

| Key | Requires when `true` |
|-----|----------------------|
| `enable_hops: true` | `hops_malt_index`, `hops_resources`, HOPS config under `config/` |
| `enable_decom: true` | `decOM_sources` directory; optional `decom_memory` (default `64GB`), `decom_threads` (default `8`) |

### HOPS parallel MALT (`hops_parallel: true`)

By default HOPS runs as **one cohort job** (`-m full`) on all samples. With parallel mode:

1. **MALT** runs per biological sample, then **MaltExtract + postprocessing** run once (`-m me_po`).
2. With **`hops_malt_mmap: true`** (recommended on a single large node): one batch job warms the MALT index into the Linux page cache, then runs `malt-run --memoryMode map` for N samples concurrently. Processes share the cached index instead of each loading a private copy (`memoryMode=load` via `hops -m malt`).

| Key | Default | Meaning |
|-----|---------|---------|
| `hops_parallel` | `false` | Enable parallel MALT (per-sample HOPS or mmap batch) |
| `hops_malt_mmap` | `false` | Use page-cache warmup + `malt-run --memoryMode map` batch (`scripts/run_parallel_malt.py`) |
| `hops_parallel_jobs` | `2` | Max concurrent MALT jobs (`N_PARALLEL` in mmap mode) |
| `hops_threads_per_job` | `threadsMalt / hops_parallel_jobs` | Threads per MALT job |
| `hops_heap_gb` | `800` | Java `-Xmx` per malt-run / HOPS invocation (GB) |
| `hops_max_memory_malt_per_job` | `maxMemoryMalt / hops_parallel_jobs` | Written to `config_hops_custom.txt` when using `hops -m malt` |

**mmap mode (single node only):** `vmtouch` is included in `workflow/envs/hops.yaml` for fast index warm-up. Total cores ≈ `hops_parallel_jobs × hops_threads_per_job`. Logs: `logs/hops_malt/parallel.log` (batch) and `logs/hops_malt/{bio}.log` (per sample).

**Per-sample HOPS mode (`hops_malt_mmap: false`):** limit Snakemake concurrency with `--resources hops_jobs=N`.

**Important:** HOPS refuses to start MALT if `maxMemoryMalt` exceeds `hops_heap_gb` when using `hops -m malt`. `create_hops_config.py` auto-splits and caps memory for that path.

Example (mmap batch on one 42-core node):

```yaml
enable_hops: true
hops_parallel: true
hops_malt_mmap: true
hops_parallel_jobs: 4
hops_threads_per_job: 10
hops_heap_gb: 400
```

---

## Pathogen screening → mapping (`evalue` rule)

After KrakenUniq, the pipeline runs the **`evalue`** rule:

- **Tool:** `scripts/dExp_Escore.py`
- **Outputs:** `results/pathogen/{bio}/evalue/pathogen/{bio}_pathogen.csv` (plus genus/species tables)

The CSV still contains an **`Escore`** column from the script; the **directory and rule name are `evalue/`**.

**Who gets reference-mapped?** At the **`generate_pathogen_targets`** checkpoint, `scripts/build_pathogen_target_pairs.py` reads those CSVs (and optional HOPS):

| `pathogen_detection_criteria` | Behaviour |
|--------------------------------|-----------|
| `map_all_escore_pathogens: true` (default) | Every taxon in `{bio}_pathogen.csv` with a spreadsheet row → mapping target. |
| `map_all_escore_pathogens: false` | Only taxa passing `escore_threshold` / `reads_threshold` (or Guellil E-value if `use_evalue_for_detection: true`). |
| HOPS enabled | Union with HOPS-positive taxa for that sample. |

Mapping starts only after **`metagenomics_screening_cohort_ready.json`** confirms **all** biological samples finished Kraken + E-value.

---

## decOM troubleshooting

decOM uses **Dask** and allocates about **`decom_memory × decom_threads`** RAM (not `decom_memory` total). Legacy PIGSTI used **`64GB × 8`** (~512 GB budget) for full-cohort runs. Lower values (e.g. `8GB × 2`) cause Dask workers to die on large cohorts.

Defaults: `decom_memory: 64GB`, `decom_threads: 8`, **`decom_fail_on_error: false`**. When `enable_hops: true`, `decom_run` waits for the HOPS heatmap so parallel MALT jobs do not compete for RAM. Limit concurrent decOM with `--resources decom=1` (default).

Before decOM, `gather_sinks` checks that `results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz` exists.

| If decOM keeps failing | Action |
|------------------------|--------|
| Workers died / OOM | Restore `decom_memory: 64GB`, `decom_threads: 8`; run after HOPS finishes |
| HOPS + decOM on same node | Ensure HOPS completes first (automatic when `enable_hops: true`) |
| decOM not required | `enable_decom: false` |
| Continue without decOM | `decom_fail_on_error: false` (writes `decOM_out/FAILED.txt`, rest of pipeline continues) |

---

## Optional tuning (safe to leave as template)

- `cleanup_intermediates`, `dedup_tool`, `merge_tool`, `markdup_threads`, `summary_threads`
- `fastq_screen_exe`, `fastq_screen_conf`, `picard_cmd`
- `build_results_catalog`, `build_results_catalog_symlinks`
- `defaults.min_reads`, `defaults.min_escore` (legacy spreadsheet defaults)
- `pathogen_detection_criteria.*` thresholds

---

## Minimal config.yaml example

```yaml
kraken_db: "/data/krakenuniq/mydb"
host_index: "/data/refs/chimera/chimera"   # prefix only; .1.bt2l on disk

host_aligner: bwa
pathogen_aligner: bwa
pathogen_mapping_mode: default

bwa_indices:
  Pig: "/data/refs/sus_scrofa.fa"
  Human: "/data/refs/human.fa"
mtDNA_indices:
  Pig: "/data/refs/sus_scrofa_mito.fa"
  Human: "/data/refs/human_mito.fa"

enable_hops: false
enable_decom: false
pathogen_screening_only: false

pathogen_detection_criteria:
  map_all_escore_pathogens: true
  use_evalue_for_detection: true
  escore_threshold: 5
  reads_threshold: 50
```

---

## See also

- [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md) — result paths  
- [`METRICS_DEFINITIONS.md`](METRICS_DEFINITIONS.md) — authentication scoring  
