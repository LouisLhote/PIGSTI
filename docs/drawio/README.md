# PIGSTI Draw.io diagram maintainer guide

## Primary editable diagrams

| File | Style |
|------|--------|
| **[`pigsti_workflow.drawio`](pigsti_workflow.drawio)** | Box / flowchart (stations as rounded rectangles) |
| **[`pigsti_subway_workflow.drawio`](pigsti_subway_workflow.drawio)** | Subway / metro map (**nf-metro light** aesthetic: soft gray panels, dark station markers, coloured tracks) |

Palette variants (same **light** style + icons): `pigsti_subway_{egypt,johnson,derain,hokusai3,java,archambault,hiroshige}.{drawio,png,pdf,svg}` and `*_light.svg` / `*_light.pdf` (high-quality PDF from the light SVG).

Open in [diagrams.net](https://app.diagrams.net/) or the Draw.io desktop app.

Regenerate from code after pipeline changes:

```bash
python scripts/rasterize_icons.py          # SVG → PNG (needed for icons in PDF/SVG)
python scripts/build_drawio_workflow.py     # box diagram
python scripts/build_drawio_subway.py       # subway draw.io (default Hiroshige)
python scripts/build_drawio_subway.py --all-palettes
python scripts/render_subway_workflow.py --all-palettes   # PNG/PDF/SVG with icons
```

### Load nf-core shapes

1. **File → Open Library from → Device**
2. Select [`../nf-core_components.xml`](../nf-core_components.xml)
3. Replace rounded boxes with official nf-core station / track shapes

### Host identification icons

The **Host identification** box has a sticky note placeholder (`🐄 🐷 🐑 🐕`). Replace with your species icons (cow, pig, sheep, dog, etc.) once designed.

### Correct route logic (do not regress)

| Route | Colour | Steps |
|-------|--------|-------|
| Preprocessing | slate `#475569` | Raw → Trim |
| Host route | blue `#376795` | Host ID → map → merge → Qualimap · DamageProfiler · **Endogenous DNA** · **Sexing** (opt) |
| Metagenomics screening | indigo `#5B4A8A` | PRINSEQ → … → KrakenUniq + optional MALT / decOM |
| Pathogen detection | coral `#D4483B` | **Pathogen detection** (E-value ∪ MaltExtract) |
| Pathogen mapping | rose `#8C1C2B` | Reference mapping (BWA / Bowtie2) |
| Authentication | deep rose `#6B1420` | Composite score panel (ANI · Entropy · Breadth · Damage · Edit dist · Mapping · Genus) |
| Cohort outputs | green `#059669` | Host QC (Qualimap, DamageProfiler) + pathogen auth **converge only here** |
| Optional | gray dashed `#94A3B8` | MALT, MaltExtract, decOM |
| Label style | purpose = route colour (bold) · tool = teal `#0D9488` (italic) | In every station box |

---

## Curated views (use these in papers / README)

| Asset | Role |
|-------|------|
| [`docs/pipeline_overview.html`](../pipeline_overview.html) | Single page: branded layout, embedded SVG, interactive Mermaid |
| [`docs/images/pipeline_workflow.svg`](../images/pipeline_workflow.svg) | Publication-style schematic (matplotlib) |
| [`docs/images/pipeline_workflow.png`](../images/pipeline_workflow.png) | Raster version of the same |
| [`docs/images/pipeline_workflow.mmd`](../images/pipeline_workflow.mmd) | Mermaid source for GitHub/GitLab render |

Editable sources: `pigsti_workflow.drawio` and `pigsti_subway_*.drawio` in this folder.

---

## Draw.io checklist (mirror the real pipeline)

Use this when adding or rearranging boxes so the drawing does not drift from Snakemake.

### Per PCR (`results/libraries/{pcr}/`)

- [ ] **AdapterRemoval** (trim) → **collapsed** reads
- [ ] **FastQ Screen** runs on **AdapterRemoval collapsed** reads (before PRINSEQ++)
- [ ] **Host** and **mtDNA** mapping inputs = collapsed + species from FastQ Screen
- [ ] **PRINSEQ++** (`prinseq/` `*-passed.fq.gz`) feeds the **pathogen / metagenomics** branch only
- [ ] Optional per-PCR QC: Qualimap on mapping, **DamageProfiler** on host BAM

### Biological sample merges

- [ ] **`results/host/`** merged host / mtDNA BAMs (`merge_*` rules)
- [ ] **`results/samples/{bio}/`** merged Qualimap (and **sexing** under `samples/.../sexing/` when enabled)

### Pathogen / metagenomics fork

- [ ] **PRINSEQ-passed** reads → Bowtie2 to host index → **unaligned** BAM/FASTQ (pathogen branch only)
- [ ] Merge unaligned **`results/pools/unaligned_fastq/`** (per cohort sample)
- [ ] **KrakenUniq** → **E-value / dExp** cohort step → **checkpoint** (`pathogen_targets.txt`) → pathogen **BWA/Bowtie2** mapping

### Optional branches (config / tooling)

- [ ] **HOPS** (`enable_hops`)
- [ ] **decOM** (`enable_decom`)
- [ ] **Genetic sexing** (`enable_sexing`, Cow / Goat / Sheep / Dog after merged host BAM)
- [ ] **`pathogen_screening_only`** (host/mtDNA lane off)
- [ ] **`pathogen_mapping_mode`** (default Bowtie2-unaligned vs host-unmapped “super_careful”)
- [ ] **`results_root`** / symlinked `results/` (outputs on large disk)

### Layout path reminder (single source of truth)

See [`OUTPUT_SCHEMA.md`](../OUTPUT_SCHEMA.md). Key folders: `libraries/`, `pools/`, `host/`, `samples/`, `metagenomics/`, `pathogen/`, `final/`, `workflow/`.

---

## Related docs

- [`CONFIG.md`](../CONFIG.md) — feature switches and paths
- [`METRICS_DEFINITIONS.md`](../METRICS_DEFINITIONS.md) — authentication scoring
- `Snakefile` — definitive rule ordering
