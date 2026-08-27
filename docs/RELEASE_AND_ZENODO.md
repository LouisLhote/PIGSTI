# GitHub release & Zenodo deposit

Checklist for publishing **PIGSTI** on GitHub and archiving a citable version on Zenodo.

## 0. Preflight (do this first)

- [ ] No machine paths in tracked files (`/raid_md0/…`, home directories, passwords).
- [ ] `config/config.yaml`, `config/samples.tsv`, and `config/Pathogen_spreadsheet.csv` are **gitignored** (local only).
- [ ] Templates exist and use placeholders: `config.example.yaml`, `samples.example.tsv`, `Pathogen_spreadsheet.example.csv`.
- [ ] `LICENSE` is MIT and matches `CITATION.cff`.
- [ ] `CHANGELOG.md` has an entry for this version.
- [ ] Dry-run succeeds on a clean clone with example configs filled for a tiny test set.

```bash
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/Pathogen_spreadsheet.example.csv config/Pathogen_spreadsheet.csv
# edit paths, then:
snakemake -n -p --use-conda --cores 4
```

## 1. GitHub repository

1. Push the publication branch to `https://github.com/LouisLhote/PIGSTI` (or your fork).
2. Set repository description: *Ancient DNA pathogen screening with Snakemake — host/mtDNA, KrakenUniq, authentication, cohort reports.*
3. Topics: `ancient-dna`, `snakemake`, `metagenomics`, `pathogen`, `krakenuniq`, `bioinformatics`.
4. Enable Issues for bug reports and questions.
5. Confirm README figures render (`docs/images/*.png`).

## 2. Create a GitHub Release (triggers Zenodo DOI if linked)

1. Update version strings in sync:
   - `CITATION.cff` → `version`
   - `CHANGELOG.md` → new section
   - optional: git tag message
2. Commit and push `main` (or release branch).
3. Tag and release:

```bash
git tag -a v1.0.0 -m "PIGSTI v1.0.0"
git push origin v1.0.0
```

4. On GitHub → **Releases** → **Draft a new release** from `v1.0.0`.
   - Title: `PIGSTI v1.0.0`
   - Paste the matching `CHANGELOG.md` section
   - Attach nothing huge (no FASTQs, no Kraken DB)

## 3. Link Zenodo ↔ GitHub (one-time)

1. Sign in at [zenodo.org](https://zenodo.org) (or [sandbox.zenodo.org](https://sandbox.zenodo.org) for a practice deposit).
2. **GitHub** → enable the `LouisLhote/PIGSTI` repository.
3. The next GitHub Release creates a Zenodo deposit and assigns a **DOI**.
4. After the first DOI exists, paste it into `CITATION.cff` as:

```yaml
identifiers:
  - type: doi
    value: 10.5281/zenodo.XXXXXXX
```

and update the README citation block.

## 4. What users should cite

Prefer, in order:

1. The paper (when published)
2. The Zenodo DOI for the **exact release** used
3. The GitHub commit SHA recorded in `results/final/run_manifest.json`

## 5. Do **not** upload to Zenodo / GitHub

| Content | Why |
|---------|-----|
| Raw FASTQ / BAM / large intermediates | Size + ethics / controlled access |
| KrakenUniq / MALT databases | Huge; distribute separately |
| Local `config.yaml` with absolute paths | Not portable; may leak structure |
| Real `samples.tsv` with private sample IDs / paths | Privacy |

Ship **code + docs + example configs + figures** only.

## 6. After publication

- [ ] Pin the release DOI in the manuscript Data/Code availability section
- [ ] Archive the paper’s exact `run_manifest.json` as supplementary material
- [ ] Bump version in `CITATION.cff` / `CHANGELOG.md` for the next development cycle
