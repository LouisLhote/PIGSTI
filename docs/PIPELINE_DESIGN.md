# PIGSTI pipeline design (v3)

Visual reference for the Snakemake workflow and `results/` layout. Regenerate the figure with:

```bash
python scripts/build_drawio_workflow.py   # editable draw.io source
python scripts/render_subway_workflow.py   # nf-metro + matplotlib figures
python scripts/render_workflow_diagram.py  # box-style + Mermaid
```

Outputs:

- [`images/pigsti_subway_workflow.svg`](images/pigsti_subway_workflow.svg) / [`.png`](images/pigsti_subway_workflow.png) / [`.html`](images/pigsti_subway_workflow.html) — **paper figure:** official **nf-core metro-map** via [nf-metro](https://seqeralabs.github.io/nf-metro/) with four sections (Preprocessing · Host · Metagenomics · Pathogen)
- [`images/pigsti_subway_workflow.mmd`](images/pigsti_subway_workflow.mmd) — editable nf-metro source
- [`NF_CORE_WORKFLOW_FIGURE.md`](NF_CORE_WORKFLOW_FIGURE.md) — how to use nf-core draw.io shapes + Google Docs paste
- [`nf-core_components.xml`](nf-core_components.xml) — official nf-core station/track library for draw.io
- [`pipeline_overview.html`](pipeline_overview.html) — hero layout, static SVG, full Mermaid graph, Snakemake DAG how-to
- [`images/pipeline_workflow.png`](images/pipeline_workflow.png) — box-style schematic (PNG)
- [`images/pipeline_workflow.svg`](images/pipeline_workflow.svg) — same (SVG, zoom in docs)
- [`images/pipeline_workflow.mmd`](images/pipeline_workflow.mmd) — Mermaid source (renders on GitHub like nf-core pipelines)
- [`images/pipeline_workflow.html`](images/pipeline_workflow.html) — legacy multi-panel Mermaid (see overview for the full graph)
- [`images/snakemake_dag.svg`](images/snakemake_dag.svg) — **optional:** rule dependency graph from `snakemake --dag` (run `bash scripts/render_snakemake_dag.sh` after Graphviz is installed)
- Editable Draw.io: [`drawio/pigsti_workflow.drawio`](drawio/pigsti_workflow.drawio) — open in [diagrams.net](https://app.diagrams.net/) (see [`drawio/README.md`](drawio/README.md))
- Legacy Draw.io export: [`PIGSTI.drawio.html`](../PIGSTI.drawio.html)

---

## Results tree (single source of truth)

```text
results/
├── libraries/{pcr}/     adapter_removal, fastq_screen, host_mapping/, mtdna_mapping/, pathogen_mapping/
├── pools/                unaligned_fastq/ (metagenomics input)
├── host/                 host_mapping/, mtdna_mapping/ (merged per bio sample)
├── samples/{bio}/        merged host/mtDNA qualimap
├── metagenomics/         krakenuniq/, hops/, decOM/, kraken_abundance/
├── pathogen/{bio}/       evalue/, pathogen_mapping/, summary/, comparison/
├── final/                cohort Excel, heatmaps, pipeline report, manifest
└── workflow/             checkpoints, pathogen_mapping_complete.txt
```

Directory names are **aligner-neutral**. Tool choice is only in `config.yaml` (`host_aligner`, `pathogen_aligner`).

---

## End-to-end flow (nf-core / rnaseq style)

Like [nf-core/rnaseq](https://github.com/nf-core/rnaseq), each step shows **Function** (what the step does) and **Tool** (software). Layout:

1. **Spine (all libraries):** raw FASTQ → adapter trimming  
2. **Fork after trim:** host route (collapsed reads) **and** metagenomics route (PRINSEQ++) in parallel  
3. **Host route:** host identification (FastQ Screen) → host & mtDNA mapping → merge → QC  
4. **Metagenomics:** PRINSEQ++ → Bowtie2 unaligned → pool → KrakenUniq  
5. **Pathogen detection:** Guellil E-value (from Kraken) **∪** HOPS (from pooled FASTQ) → checkpoint targets  
6. **Pathogen mapping:** reference mapping → multi-criteria authentication → cohort reports in `results/final/`

| Spine step | Function | Tool |
|------------|----------|------|
| Input | Raw sequencing data | FASTQ · `samples.tsv` |
| Trim | Adapter & quality trimming | **AdapterRemoval** (PE) · **cutadapt** (SE) |

Regenerate: `python scripts/render_workflow_diagram.py` → PNG/SVG + `pipeline_workflow.mmd`.

### Host & mtDNA route

| Function | Tool | Results path |
|----------|------|----------------|
| Host identification | **FastQ Screen** (collapsed reads; species icons TBD) | `results/libraries/{pcr}/fastq_screen/` |
| Host genome mapping | BWA / Bowtie2 (`host_aligner`) | `results/libraries/{pcr}/host_mapping/` (collapsed reads) |
| Mitochondrial mapping | BWA / Bowtie2 | `results/libraries/{pcr}/mtdna_mapping/` (collapsed reads) |
| Merge host / mtDNA libraries | samtools | `results/host/host_mapping/`, `results/host/mtdna_mapping/` |
| Host / mtDNA mapping QC | Qualimap | `results/samples/{bio}/qualimap/` |

### Metagenomics route

| Function | Tool | Results path |
|----------|------|----------------|
| Complexity filter & dedup | **PRINSEQ++** | `results/libraries/{pcr}/prinseq/` |
| Separate unaligned reads | Bowtie2 · samtools | `unaligned_fastq/` (PRINSEQ-passed input) |
| Merge unaligned FASTQs | cat · pigz | `results/pools/unaligned_fastq/` |
| Taxonomic classification | KrakenUniq 1.0.4 | `results/metagenomics/krakenuniq/{bio}/` |

Optional: **decOM** (source screening).

### Pathogen route

| Function | Tool | Results path |
|----------|------|----------------|
| Kraken E-value scoring | `dExp_Escore.py` (Guellil E-value) | `results/pathogen/{bio}/evalue/` |
| aDNA screening (union with E-value) | **HOPS** (MALT · MaltExtract) | `results/metagenomics/hops/` |
| Select mapping targets | cohort checkpoint (E-value ∪ HOPS) | `results/workflow/pathogen_targets.txt` |
| Pathogen reference mapping | BWA / Bowtie2 (`pathogen_aligner`) | `results/pathogen/{bio}/pathogen_mapping/` |
| Pathogen mapping QC | ANI · DamageProfiler · Qualimap | `results/pathogen/{bio}/` |

Optional: **HOPS** (`enable_hops`) feeds pathogen **detection** together with E-value (union at checkpoint).

**Per-PCR order:** AdapterRemoval → fork: (A) FastQ Screen → host/mtDNA mapping on collapsed reads; (B) PRINSEQ++ → Bowtie2-unaligned → metagenomics/pathogen. Host species comes from FastQ Screen; host mapping does **not** use PRINSEQ-passed reads.

Open [`pipeline_overview.html`](pipeline_overview.html) in a browser for the main diagram page (static + interactive + DAG instructions).

---

## Optional branches

| Config flag | Effect |
|-------------|--------|
| `pathogen_screening_only: true` | Skips host/mtDNA lanes; metagenomics + pathogen only |
| `enable_hops: true` | HOPS + comparison plots; detection-score matrices |
| `enable_decom: true` | decOM source screening → `metagenomics/decOM/decOM_out/` |
| `pathogen_mapping_mode: super_careful` | Pathogen reads from BWA host-unmapped per PCR |
| `host_aligner` / `pathogen_aligner` | `bwa` or `bowtie2` (paths unchanged) |

---

## Checkpoint gate

Pathogen reference mapping starts only after:

1. All sample E-value pathogen CSVs exist  
2. Kraken reports exist  
3. Optional HOPS heatmap (if enabled)  

Barrier file: `results/workflow/pathogen_mapping_complete.txt`  
Rule: `pathogen_mapping_targets`

See [`CHECKPOINT_PATHOGEN_MAPPING.md`](CHECKPOINT_PATHOGEN_MAPPING.md).

---

## Interactive view

Open **`pigsti-pipeline-design.canvas.tsx`** in Cursor (Canvases) for a zoomable DAG of the same logic.
