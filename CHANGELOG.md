# Changelog

All notable changes to PIGSTI are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/).

## [1.0.0] — 2026-08-26

### Added
- Full Snakemake workflow for host/mtDNA, metagenomics, and pathogen authentication
- Interactive HTML config generator (`config/pigsti_config_generator.html`)
- Startup validation (`scripts/validate_pigsti_setup.py`)
- Guellil E-value as default detection filter (optional E-score)
- Always-on edit-distance damage vs no-damage split
- On-demand pathogen BWA index build before `bwa aln`
- Optional HOPS (including parallel / mmap MALT), decOM, genetic sexing, Multi-QC dashboard
- Metro-map and Function+Tool workflow figures
- Zenodo / GitHub release metadata (`CITATION.cff`, `.zenodo.json`, `LICENSE`)

### Changed
- Portable `config/config.example.yaml` (no machine-specific paths)
- Picard and FastQ Screen taken from conda environments only

[1.0.0]: https://github.com/LouisLhote/PIGSTI/releases/tag/v1.0.0
