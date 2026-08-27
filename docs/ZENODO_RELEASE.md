# Zenodo release guide — PIGSTI v1.0.0

This document describes how to publish **PIGSTI** (pipeline code) and **paper supplementary materials** on Zenodo for the F1000Research software article.

**Paper title (current):**  
*PIGSTI: a modular, reproducible pipeline for detecting species identity, pathogens, and microbes from animal palaeogenomic data*

**Repository:** [github.com/LouisLhote/PIGSTI](https://github.com/LouisLhote/PIGSTI)

---

## Recommended strategy (two related deposits)

For F1000 software papers it is usually clearest to use **two Zenodo records** that cross-link:

| Record | Contents | Typical resource type | License |
|--------|----------|----------------------|---------|
| **A — Software** | GitHub release `v1.0.0` (Snakemake pipeline) | Software | MIT (see `LICENSE`) |
| **B — Supplementary data** | Tables, configs, manifests, figures, simulation inputs, etc. | Dataset / Other | *TBD — often CC-BY 4.0* |

Link them in metadata (“Related identifiers” / description) and cite both in the paper (Software availability + Supplementary material).

**Alternative:** one manual Zenodo upload with a zip containing `PIGSTI/` + `supplementary/` — simpler for reviewers, but the software DOI will not auto-update on future GitHub releases.

---

## Part A — Pipeline (GitHub → Zenodo)

### 1. Prerequisites

- [ ] GitHub repo public: `LouisLhote/PIGSTI`
- [ ] `CITATION.cff` and `LICENSE` on `main` (already present)
- [ ] Tag matches version in `CITATION.cff` (`1.0.0`)

### 2. Enable GitHub–Zenodo integration (once)

1. Sign in at [zenodo.org](https://zenodo.org) with GitHub.
2. **Account → GitHub → Enable** for `LouisLhote/PIGSTI`.
3. Zenodo will create a draft deposit on each **GitHub Release** (not on every push).

### 3. Create GitHub release `v1.0.0`

From the publication repo root (after all doc changes are committed):

```bash
git tag -a v1.0.0 -m "PIGSTI v1.0.0 — public release for F1000 software article"
git push origin v1.0.0
```

On GitHub: **Releases → Draft a new release** → choose tag `v1.0.0`, title `v1.0.0`, paste release notes (see template below).

### 4. Publish Zenodo software record

1. Open the new draft on Zenodo (linked from GitHub release or Zenodo uploads).
2. Verify metadata (edit if needed):

| Field | Suggested value |
|-------|-----------------|
| **Upload type** | Software |
| **Title** | Same as `CITATION.cff` title |
| **Authors** | *Fill from paper author list — see § Authors below* |
| **Description** | Abstract from `CITATION.cff` + one sentence on Snakemake/conda |
| **License** | MIT |
| **Keywords** | ancient DNA, palaeogenomics, Snakemake, pathogen authentication, … |
| **Version** | 1.0.0 |
| **GitHub repository** | `https://github.com/LouisLhote/PIGSTI` |

3. **Publish** → copy the **concept DOI** (e.g. `10.5281/zenodo.XXXXXXX`) and version DOI.
4. Update `CITATION.cff` (`doi:` field) and README citation line when the DOI is assigned.

### Software release notes template

```markdown
## PIGSTI v1.0.0

First public release accompanying the F1000Research software article.

### Includes
- Snakemake workflow (host/mtDNA, metagenomics, pathogen authentication)
- Interactive config generator (`config/pigsti_config_generator.html`)
- Documentation (`docs/CONFIG.md`, `docs/METRICS_DEFINITIONS.md`)

### Requirements
- Linux/macOS, Conda/Mamba, Snakemake ≥ 7.32
- External databases: KrakenUniq (aMETA NT), host/chimera indices, pathogen references

### Citation
See `CITATION.cff`. Record the git tag or commit SHA in `results/final/run_manifest.json` for each run.
```

---

## Part B — Paper supplementary materials

### Folder layout (suggested)

Prepare a directory **`PIGSTI_paper_supplementary/`** (can live outside the git repo or in a separate private repo until upload):

```
PIGSTI_paper_supplementary/
├── README.md                 # Manifest: file list + one-line description each
├── METADATA.csv              # Optional: sample × pathogen hit table for paper
├── configs/                  # Frozen configs used for 953-sample screen (if shareable)
├── manifests/                # samples.tsv subsets, pathogen spreadsheet snapshot
├── simulation/               # Mock metagenome parameters / spike-in design (Figure 2)
├── phylogeny/                # PathPhynder inputs or tree summaries (if not too large)
├── figures/                  # Source data for Fig 2–3 (CSV/plot scripts)
└── tables/                   # Supplementary tables (SX) as CSV/XLSX
```

> **Do not** include raw FASTQ, full Kraken DB, or private sample paths unless ethics/data policy allows.

### Supplementary README manifest (template)

Fill this when files are final:

| File | Description | Paper reference |
|------|-------------|-----------------|
| `tables/SX_pathogen_panel.csv` | Pathogen reference panel used for empirical screen | Supp. Table SX |
| `tables/SX_sample_metadata.csv` | 953-sample metadata (published + new) | Supp. Table SX |
| `simulation/spike_in_design.tsv` | Ten pathogens × six abundance tiers (12M reads) | Methods — simulation |
| … | … | … |

### Manual Zenodo upload (supplementary record)

1. Zip `PIGSTI_paper_supplementary/` or upload files individually.
2. **Upload type:** Dataset (or “Other” if mixed).
3. **Title:** e.g. *Supplementary data for: PIGSTI — … (F1000Research)*
4. **Related identifier:** DOI of software record (Part A), relation “Is supplement to” / “Is documented by”.
5. **License:** agree with journal policy (often **CC-BY 4.0** for data tables).
6. Publish and add DOI to the paper supplementary-materials section.

---

## Authors (complete for Zenodo)

*Replace placeholders when `PIGSTI_paper_Metadatas` / author list is final.*

| # | Name | Affiliation | ORCID | Creator on Zenodo? |
|---|------|-------------|-------|-------------------|
| 1 | Louis Lhote | Trinity College Dublin | | Yes |
| 2 | | | | |
| … | | | | |

**Zenodo tip:** order creators as on the paper; mark corresponding author in description if Zenodo author order is fixed.

---

## Funding & grants (optional Zenodo field)

| Funder | Award / grant ID |
|--------|------------------|
| | |

---

## Post-release checklist

- [ ] Software DOI in `CITATION.cff` (`doi:` key)
- [ ] Software DOI badge/line in `README.md` (optional)
- [ ] Supplementary DOI in paper “Data availability”
- [ ] F1000 “Software availability” box: GitHub URL + version tag + software DOI
- [ ] F1000 supplementary: supplementary Zenodo DOI or file list
- [ ] Confirm `run_manifest.json` in pipeline documents commit/tag used for paper reruns

---

## Optional: `.zenodo.json` (GitHub-only deposit)

If you use **GitHub → Zenodo** for Part A only, you can add `.zenodo.json` at the repo root to pre-fill Zenodo fields on each release. Do **not** put large supplementary files in the git repo for this path.

Example skeleton (authors to be completed):

```json
{
  "title": "PIGSTI: a modular, reproducible pipeline for detecting species identity, pathogens, and microbes from animal palaeogenomic data",
  "description": "<paste CITATION.cff abstract>",
  "upload_type": "software",
  "publication_date": "2026-08-26",
  "creators": [
    { "name": "Lhote, Louis", "affiliation": "Trinity College Dublin", "orcid": "" }
  ],
  "keywords": [
    "ancient DNA",
    "palaeogenomics",
    "pathogen detection",
    "Snakemake",
    "KrakenUniq"
  ],
  "license": "mit",
  "version": "1.0.0",
  "notes": "Companion supplementary dataset: <Supplementary Zenodo DOI when available>"
}
```

---

## What to send next (so the manifest can be filled in)

Point the assistant (or edit § Part B yourself) to:

1. **Path** to supplementary folder or spreadsheet (`PIGSTI_paper_Metadatas`, Supp. tables, etc.)
2. **Final author list** + ORCIDs
3. **Which files are public** vs. embargoed / metadata-only
4. **Single Zenodo record vs. two** (preference)
5. **Supplementary license** (CC-BY 4.0 or other)

Once those are available, the manifest table, `.zenodo.json` creators block, and GitHub release text can be completed automatically.
