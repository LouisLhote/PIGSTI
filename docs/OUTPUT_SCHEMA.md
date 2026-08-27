# PIGSTI output naming schema (v3)

See [`V2_MIGRATION.md`](V2_MIGRATION.md) for older layout mappings.

**Pipeline visuals:** [`PIPELINE_OVERVIEW.md`](PIPELINE_OVERVIEW.md), [`pipeline_overview.html`](pipeline_overview.html), optional `docs/images/snakemake_dag.svg`.

## Tree overview

```
results/
  libraries/{pcr}/              # adapter removal, fastq_screen, host_mapping/, mtdna_mapping/
  samples/{bio}/                  # merged host/mtDNA QC (qualimap on merged BAMs)
  pools/unaligned_fastq/          # merged unaligned FASTQs (Kraken/decOM input)
  host/                           # merged host and mtDNA BAMs per biological sample
    host_mapping/
    mtdna_mapping/
  metagenomics/
    krakenuniq/{bio}/             # KrakenUniq report + output per sample
    hops/                         # optional HOPS
    decOM/                        # optional decOM (p_sink, p_keys/; tool output in decOM/decOM_out/)
    kraken_abundance/             # cohort matrices
    pathogen_detection/           # optional detection-score matrix (HOPS mode)
  pathogen/{bio}/
    evalue/                       # genus, species, pathogen CSVs
    pathogen_mapping/             # alignments + QC metrics per pathogen
    summary/                      # per-pathogen PDFs, pathogen_summary CSV, sample_report PDF
    comparison/                   # HOPS comparison TSV/HTML (optional)
  final/                          # cohort Excel, heatmaps, manifest, catalog
  workflow/                       # checkpoints, barriers
  browse/                         # optional symlinks
```

## Wildcards

| Rule context | `{sample}` means | Base path |
|--------------|------------------|-----------|
| PCR-level | Library ID (`pcr`) | `results/libraries/{pcr}/` |
| Bio-level (pathogen, kraken) | Biological sample ID | `results/pathogen/{bio}/` or `results/metagenomics/...` |
| Merged host QC | Biological sample ID | `results/samples/{bio}/` |

## Key paths

| Artifact | Path |
|----------|------|
| Kraken report | `results/metagenomics/krakenuniq/{bio}/{bio}_kraken-report.txt` |
| E-value pathogen CSV | `results/pathogen/{bio}/evalue/pathogen/{bio}_pathogen.csv` |
| Pathogen BAM | `results/pathogen/{bio}/pathogen_mapping/{bio}_{safe}.dedup.bam` |
| Pathogen PDF | `results/pathogen/{bio}/summary/{bio}_{safe}_pathogen_report.pdf` |
| Merged host qualimap | `results/samples/{bio}/qualimap/genome_results.txt` |
| Merged mtDNA qualimap | `results/samples/{bio}/qualimap_mtdna/genome_results.txt` |
| Cohort pathogen Excel | `results/final/pathogen_summary_all_samples.xlsx` |
| Kraken input FASTQ | `results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz` |
| Host dedup BAM (PCR) | `results/libraries/{pcr}/host_mapping/{pcr}.dedup.bam` |
| mtDNA dedup BAM (PCR) | `results/libraries/{pcr}/mtdna_mapping/{pcr}.dedup.bam` |
| Merged host BAM | `results/host/host_mapping/{bio}.dedup.merged.bam` |
| Merged mtDNA BAM | `results/host/mtdna_mapping/{bio}.dedup.merged.bam` |
| Kraken abundance matrices | `results/metagenomics/kraken_abundance/` |
| Mapping barrier | `results/workflow/pathogen_mapping_complete.txt` |
| Sexing plot (bio sample) | `results/samples/{bio}/sexing/{bio}_sexing.pdf` |
| Sexing plot | `results/samples/{bio}/sexing/{bio}_sexing.pdf` |
| Sexing TSV | `results/samples/{bio}/sexing/{bio}_sexing.tsv` |

Sexing columns (`sexing_call`, `sexing_plot`, …) appear in `results/final/host_mtdna_summary_all_samples.xlsx` and `comprehensive_summary_all_samples.xlsx` when `enable_sexing: true` (Cow, Goat, Sheep, Dog).

Aligner (`host_aligner` / `pathogen_aligner`: `bwa` or `bowtie2`) is config-only; directory names are aligner-neutral.

## Helpers

- `scripts/pigsti_paths.py` — canonical path functions
- `scripts/pigsti_naming.safe_pathogen_name()` — filename-safe pathogen tokens

CSV column **`Escore`** from `dExp_Escore.py` is unchanged.
