# PIGSTI pipeline overview

Concise map of the workflow and configuration options. See **Figure 1** in the [README](../README.md) for the full schematic.

## Workflow (core steps)

| Step | What it does | Tool |
|------|----------------|------|
| 1. Trim | Adapter removal and quality filtering | **AdapterRemoval** (PE) · **cutadapt** (SE) |
| 2. Host ID | Best host species from collapsed reads | **FastQ Screen** |
| 3. Host map | Nuclear genome mapping | **BWA** or **Bowtie2** |
| 4. mtDNA map | Mitochondrial mapping | **BWA** or **Bowtie2** |
| 5. Merge | Combine PCR replicates per biological sample | **samtools** |
| 6. Host QC | Damage, mapping stats, optional sexing | **DamageProfiler** · **Qualimap** · sexing scripts |
| 7. Metagenomics prep | Complexity filter, unaligned reads | **PRINSEQ++** · **Bowtie2** |
| 8. Screen | Taxonomic classification | **KrakenUniq** |
| 9. Detect | Pathogen candidates | Guellil **E-value** (default) · optional **HOPS** |
| 10. Map | Reference mapping for candidates | **BWA** or **Bowtie2** |
| 11. Authenticate | Multi-criteria scoring | Edit-distance split · breadth · damage · genus rank |
| 12. Report | Cohort Excel/PDF and manifest | Pipeline summary scripts |

Dashed / optional branches: HOPS, decOM, genetic sexing, Multi-QC dashboard.

## Optional modules

| Config flag | Effect |
|-------------|--------|
| `enable_hops: true` | HOPS screening (MALT + MaltExtract); union with E-value hits |
| `hops_parallel: true` | Per-sample parallel MALT |
| `hops_malt_mmap: true` | Single-node mmap batch mode |
| `enable_decom: true` | decOM source tracking |
| `enable_sexing: true` | Residual sexing (Cow, Goat, Sheep, Dog) |
| `enable_multi_qc_dashboard: true` | Per-sample HTML QC dashboard |
| `pathogen_screening_only: true` | Skip host/mtDNA mapping and host QC |
| `cleanup_intermediates: true` | Remove large intermediates after final outputs |

## Key configuration

| Setting | Options / notes |
|---------|-----------------|
| `host_aligner` | `bwa` (default) · `bowtie2` |
| `pathogen_aligner` | `bwa` (default) · `bowtie2` |
| `pathogen_mapping_mode` | `default` · `super_careful` |
| `dedup_tool` / `merge_tool` | `samtools` · `picard` |
| `bwa_indices` / `mtDNA_indices` | Host reference FASTAs (species names match `samples.tsv`) |
| `kraken_db` | KrakenUniq database path |
| `host_index` | Bowtie2 chimera index prefix (unaligned read extraction) |
| `pathogen_spreadsheet` | CSV: KrakenUniq name, HOPS name, BWA index per pathogen |
| `results_root` | Optional output directory (CLI or env override) |

### Pathogen detection thresholds (`pathogen_detection_criteria`)

| Key | Default | Role |
|-----|---------|------|
| `use_evalue_for_detection` | `true` | Use Guellil E-value; `false` → E-score |
| `guellil_evalue_threshold` | `0.001` | E-value cutoff |
| `escore_threshold` | `5` | E-score cutoff (if E-score mode) |
| `reads_threshold` | `50` | Minimum supporting reads |
| `map_all_escore_pathogens` | `false` | Map all E-score hits vs filtered set |

Edit-distance damage split is **always on**.

## Main outputs

```text
results/
├── libraries/{pcr}/     per-library preprocessing and mapping
├── host/                merged host/mtDNA alignments
├── samples/{bio}/         per-sample QC
├── metagenomics/        KrakenUniq, HOPS, decOM
├── pathogen/{bio}/      detection, mapping, authentication
└── final/               cohort reports, heatmaps, run manifest
```

Full output schema: [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md) · config details: [`CONFIG.md`](CONFIG.md) · metrics: [`METRICS_DEFINITIONS.md`](METRICS_DEFINITIONS.md)
