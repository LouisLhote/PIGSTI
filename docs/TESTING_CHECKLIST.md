# PIGSTI testing checklist (post-audit)

Use after pulling checkpoint changes.

## 1. Dry-run DAG

```bash
cd PIGSTI   # or PIGSTI_publication
snakemake -n -p results/workflow/pathogen_mapping_complete.txt 2>&1 | tee logs/dryrun.txt
```

Check:

- [ ] `metagenomics_screening_cohort_ready` appears **after** all `evalue` / `krakenuniq` jobs
- [ ] `generate_pathogen_targets` depends on `metagenomics_screening_cohort_ready.json`
- [ ] `krakenuniq` input is `results/pools/unaligned_fastq/{bio}_unaligned.fastq.gz` (merged; no `prinseq_after_merge_unaligned` job)
- [ ] Per-PCR `prinseq` still runs under `results/libraries/{pcr}/prinseq/` before `bowtie2_unaligned`

## 2. Pathogen checkpoint (`--rerun-incomplete`)

1. Delete gate files (if recovering a bad partial run):

   ```bash
   rm -f results/workflow/metagenomics_screening_cohort_ready.json \
         results/workflow/pathogen_targets.txt \
         results/workflow/pathogen_targets.manifest.json \
         results/workflow/pathogen_mapping_complete.txt
   ```

2. Run:

   ```bash
   snakemake --rerun-incomplete -c 4
   ```

3. Confirm:

   - [ ] Checkpoint runs only when **all** `BIO_SAMPLES` have `*_pathogen.csv`
   - [ ] `pathogen_targets.manifest.json` → `cohort_samples` lists every bio sample
   - [ ] No pathogen BAM jobs start while any sample lacks E-score

## 3. Merged unaligned FASTQ → Kraken

Pick one bio sample `S`:

- [ ] `results/pools/unaligned_fastq/S_unaligned.fastq.gz` exists (merge)
- [ ] `krakenuniq` log for `S` shows `S_unaligned.fastq.gz` as input (not `.prinseq_passed.fq.gz`)

## 4. Scripts smoke test

```bash
python scripts/build_pathogen_target_pairs.py --help  # must list --evalue (and --escore alias)
python scripts/validate_pigsti_setup.py --config config/config.yaml
```
