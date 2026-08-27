# PIGSTI results layout — migration notes

## v3 (current)

| Role | Path |
|------|------|
| Per-PCR QC | `results/libraries/{pcr}/` |
| Merged host/mtDNA sample QC | `results/samples/{bio}/qualimap/`, `qualimap_mtdna/` |
| Merged unaligned FASTQs | `results/pools/unaligned_fastq/` |
| Merged host / mtDNA BAMs | `results/host/host_mapping/`, `results/host/mtdna_mapping/` |
| Kraken abundance (legacy) | `results/KRAKENUNIQ_ABUNDANCE_MATRIX/` → `results/metagenomics/kraken_abundance/` |
| KrakenUniq | `results/metagenomics/krakenuniq/{bio}/` |
| Pathogen evalue + mapping + summaries | `results/pathogen/{bio}/` |
| Cohort tables | `results/final/` |
| HOPS / decOM / abundance | `results/metagenomics/` |
| Checkpoints | `results/workflow/` |

## v2 → v3

| v2 | v3 |
|----|-----|
| `results/samples/{bio}/krakenuniq/` | `results/metagenomics/krakenuniq/{bio}/` |
| `results/samples/{bio}/evalue/` | `results/pathogen/{bio}/evalue/` |
| `results/samples/{bio}/pathogen_mapping/` | `results/pathogen/{bio}/pathogen_mapping/` |
| `results/samples/{bio}/summary/` | `results/pathogen/{bio}/summary/` |
| `results/samples/{bio}/comparison/` | `results/pathogen/{bio}/comparison/` |
| `results/pools/qualimap/{bio}/` | `results/samples/{bio}/qualimap/` |
| `results/pools/qualimap_mtdna/{bio}/` | `results/samples/{bio}/qualimap_mtdna/` |
| `results/libraries/{pcr}/bwa_host/` | `results/libraries/{pcr}/host_mapping/` |
| `results/libraries/{pcr}/bwa_mtdna/` | `results/libraries/{pcr}/mtdna_mapping/` |
| `results/pools/bwa_host/` | `results/host/host_mapping/` |
| `results/pools/bwa_mtdna/` | `results/host/mtdna_mapping/` |
| `results/pools/host_mapping/` | `results/host/host_mapping/` |
| `results/pools/mtdna_mapping/` | `results/host/mtdna_mapping/` |
| `results/workflow/pathogen_bwa_complete.txt` | `results/workflow/pathogen_mapping_complete.txt` |

Optional: `scripts/rename_aligner_output_dirs.py --move-results` only if you must keep an old `results/` tree.

## v1 → v3 (shortcut)

| v1 | v3 |
|----|-----|
| `results/{bio}/Escore/` | `results/pathogen/{bio}/evalue/` |
| `results/{bio}/krakenuniq/` | `results/metagenomics/krakenuniq/{bio}/` |
| `results/final/` | `results/final/` (unchanged) |
| `results/sample_host_mapping/` | `results/pools/host_mapping/` |

After moving directories, delete stale `results/workflow/pathogen_targets*` and run `snakemake --rerun-incomplete`.
