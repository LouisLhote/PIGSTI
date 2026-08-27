# Contributing to PIGSTI

Thanks for improving PIGSTI. Short guidelines:

## Before you open a PR

1. Copy templates — never commit local paths:
   ```bash
   cp config/config.example.yaml config/config.yaml
   cp config/samples.example.tsv config/samples.tsv
   cp config/Pathogen_spreadsheet.example.csv config/Pathogen_spreadsheet.csv
   ```
2. Dry-run after changes:
   ```bash
   snakemake -n -p --use-conda --cores 4
   ```
3. Keep secrets and machine paths out of git (`config/config.yaml`, `samples.tsv`, real spreadsheets, `results/`, `.snakemake/`).

## Style

- Prefer small, focused PRs (one bugfix or one feature).
- Document new config keys in `docs/CONFIG.md` and `config/config.example.yaml`.
- If you change authentication metrics, update `docs/METRICS_DEFINITIONS.md`.

## Reporting issues

Include: PIGSTI version/commit (`results/final/run_manifest.json` if available), Snakemake version, the failing rule name, and the relevant log under `logs/`.
