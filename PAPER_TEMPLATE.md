# PIGSTI paper write-up template

Use this as a structured checklist for your Methods and Reproducibility sections.

## 1. Pipeline identity
- Pipeline name: PIGSTI (Pathogen anImal Genome Sequence ToolkIt)
- Repository snapshot used: commit hash or release tag: `<INSERT>`
- Pipeline command: `<INSERT snakemake command>`
- Snakemake version: `<INSERT>`

## 2. Input data
1. Sequencing data source
   - Input manifest: `config/samples.tsv`
   - Sample identifiers: `sample` (bio sample), `pcr` (PCR/library), `sequencing_run`
2. Read preprocessing
   - Adapter removal: adapterremoval++ parameters: `<INSERT>`
   - QC: FastQ Screen + Prinseq/Q30 filtering: `<INSERT>`

## 3. References and databases
1. Host reference
   - BWA indices: `config/bwa_indices`
2. mtDNA references
   - mtDNA indices: `config/mtDNA_indices`
3. Pathogen references
   - Pathogen spreadsheet: `config/Pathogen_spreadsheet.csv`
   - KrakenUniq DB: `kraken_db` from `config/config.yaml`
4. Optional HOPS
   - Enabled: `enable_hops: true/false`

## 4. Detection logic and thresholds
1. Pathogen mapping mode
   - `pathogen_mapping_mode`: `<INSERT default/super_careful>`
   - `host_aligner`: `<INSERT bwa/bowtie2>`
   - `pathogen_aligner`: `<INSERT bwa/bowtie2>`
2. Primary detection criteria
   - Criteria file: `config/pathogen_detection_criteria.yaml`
   - Report the final thresholds you used, e.g.:
     - `ani_threshold`
     - `entropy_threshold`
     - `breadth_ratio_threshold`
     - `damage_5p_ct_threshold`
     - `hops_detection_threshold` (if HOPS enabled)
3. Per-pathogen spreadsheet overrides
   - Mention that `Pathogen_spreadsheet.csv` can override thresholds per taxon.

## 5. Output suite
At minimum, reference the outputs below:
- Pathogen targets (checkpoint-derived): `results/pathogen_targets.txt`
- Reproducibility manifest (file hashes + configs): `results/run_manifest.json`
- Final summaries:
  - `results/final/host_mtdna_summary_all_samples.xlsx`
  - `results/final/pathogen_summary_all_samples.xlsx`
  - `results/final/pathogen_detection_scores_heatmap.pdf`
- Per-pathogen reports:
  - `results/{sample}/summary/{sample}_{pathogen_safe}_pathogen_report.pdf`

## 6. Reproducibility statement (recommended)
“All analyses were performed using PIGSTI with a fully specified configuration. We provide the exact workflow snapshot and a run manifest (`results/run_manifest.json`) containing file hashes of the pipeline and configuration used to generate the results.”

