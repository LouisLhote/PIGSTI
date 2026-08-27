# PIGSTI scripts inventory

Scripts under `scripts/` and whether they are invoked by `Snakefile` (`PIGSTI_publication`).

## Wired in Snakemake DAG

| Script | Role |
|--------|------|
| `validate_pigsti_setup.py` | Pre-flight config/samples/paths |
| `dExp_Escore.py` | KrakenUniq → E-score tables |
| `parse_fastq_screen.py` | Best host species from FastQ Screen |
| `bwa_aln_host.py` / `bwa_aln_mtdna.py` | Host/mtDNA align (BWA or Bowtie2 branch) |
| `run_damageprofiler.py` | Host/mtDNA/pathogen damage |
| `softclip_mod.py` | Host CRAM soft-clip |
| `calculate_breadth_stats.py` | Pathogen depth/breadth |
| `calculate_entropy_profile.py` | Pathogen entropy (100bp + 1000bp + plot) |
| `calculate_edit_distance_r2.py` | Edit-distance R² + ANI distribution |
| `filter_bam_by_end_damage.py` | Damage vs no-damage BAM split |
| `calculate_genus_ranking.py` | Genus rank vs Kraken report |
| `generate_pathogen_lists.py` | HOPS/Kraken pathogen lists |
| `create_hops_config.py` | Custom HOPS config |
| `compare_pathogens.py` | Kraken vs HOPS comparison TSV/HTML |
| `krakenuniq_abundance_matrix.R` | Cross-sample abundance matrix |
| `plot_krakenuniq_abundance_matrix.R` | Abundance heatmaps |
| `summarize_pathogen_data.py` | Per-sample pathogen summary CSV |
| `summarize_pathogen_data_pcr.py` | PCR-level summary (super_careful) |
| `summarize_host_mtdna.py` | Host/mtDNA Excel summary |
| `merge_summaries.py` | Merge CSVs → Excel |
| `create_pathogen_heatmap.py` | Detection score heatmap |
| `create_comprehensive_summary.py` | Combined Excel |
| `create_comprehensive_summary_pathogen_only.py` | Screening-only Excel |
| `generate_sample_report.R` | Per-sample PDF |
| `generate_pathogen_report.R` | Per-pathogen PDF |
| `calculate_detection_scores.py` | 9-point score matrix (HOPS on) |
| `create_multi_qc_dashboard.py` | HTML QC dashboard |
| `generate_pipeline_report.py` | Run monitoring HTML |
| `create_interactive_dashboard.py` | Interactive pathogen dashboard |
| `generate_publication_figures.py` | Extra publication plots |
| `write_run_manifest.py` | Reproducibility manifest |
| `check_indices.py` | Index validation marker |
| `render_workflow_diagram.py` | README workflow figure (manual/CI) |

## Utilities (not in `rule all`)

| Script | Notes |
|--------|------|
| `ena_to_pigsti_samplesheet.py` | ENA download + samples.tsv |
| `download_ena_fastq.py` | Standalone ENA fetch |
| `make_samples_from_fastq.py` | Build samples.tsv from FASTQ dir |
| `parse_samples_tsv.py` | Samplesheet helper |
| `check_samples_loaded.py` | Debug sample loading |
| `setup_diamond_databases.py` | Diamond DB build (no Snakemake rule) |
| `parse_diamond_hits.py` | Diamond parser (no Snakemake rule) |
| `prepare_bwa_targets.py` | Legacy helper |
| `bowtie2_align.py` | Legacy helper |
| `merge_pathogen_plots.py` | Legacy plotting |
| `alternative_escore_functions.py` | R&D |
| `evaluate_escore_functions.py` | R&D (uses random data) |
| `aDNA-BAMPlotter.py` | Removed from DAG (resource heavy) |
| `remove_duplicate_samples.sh` | Maintenance |
| `add_ll155_ll166_to_samplesheet.sh` | One-off |

## Shared modules (new)

| Module | Purpose |
|--------|---------|
| `pigsti_naming.py` | Canonical `safe_pathogen_name()` |
| `pigsti_paths.py` | Path builders matching Snakefile |
