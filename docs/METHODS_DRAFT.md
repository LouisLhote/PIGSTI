# Methods (detailed draft) — PIGSTI

Style modelled on nf-core/eager (Fellows Yates *et al.*, 2021).

- **Part A** — toolkit by analytical stage.  
- **Part B** — detailed steps and parameters (exact values from the Snakefile/scripts, with explanations).  
- **Part C** — this study, as paper-style narrative.

---

# Part A — Toolkit by analytical stage

## 1. Pre-processing

- **Adapter trimming:** AdapterRemoval v2 (Schubert, Lindgreen & Orlando, 2016); cutadapt (Martin, 2011) for single-end libraries  
- **Host species identification:** FastQ Screen (Wingett & Andrews, 2018) with Bowtie2 (Langmead & Salzberg, 2012)

## 2. Host / mtDNA mapping and post-processing

- **Alignment:** BWA *aln* / *samse* (Li & Durbin, 2009) or Bowtie2 (Langmead & Salzberg, 2012)  
- **Filtering, sorting, duplicate marking:** SAMtools (Li *et al.*, 2009; Danecek *et al.*, 2021) or Picard MarkDuplicates  
- **Terminal soft-clipping (before host Qualimap):** PIGSTI `softclip_mod.py`  
- **Mapping QC:** Qualimap2 (Okonechnikov, Conesa & García-Alcalde, 2016)  
- **Damage profiling:** DamageProfiler (Neukamm, Peltzer & Nieselt, 2021)  
- **Genetic sexing (optional):** residual sexing for cattle, goat, sheep, dog

## 3. Metagenomics screening

- **Dereplication:** PRINSEQ++ (Cantu, Sadikalay & Edwards, 2019)  
- **Mammalian host reads removal:** Bowtie2 vs multi-host index; retain unmapped (SAMtools `-f 4`)  
- **Pooling of PCR replicates:** `cat` + compression  
- **Metagenomics classification:** KrakenUniq (Breitwieser, Baker & Salzberg, 2018)  
- **Optional metagenomics classification:** HOPS / MALT / MaltExtract (Hübler *et al.*, 2019; Herbig *et al.*, 2016)  
- **Optional source tracking:** decOM (Duitama González *et al.*, 2023)

## 4. Pathogen detection, mapping and authentication

- **Candidate selection:** KrakenUniq counts ± Guellil E-value / E-score vs pathogen spreadsheet (± HOPS if enabled)  
- **Pathogen alignment:** BWA *aln*/*samse* or Bowtie2  
- **Post-map QC:** SAMtools MAPQ / sort / markdup; Qualimap2; DamageProfiler  
- **Authentication:** ANI, relative entropy, breadth ratio, 5′ C→T damage, edit-distance decay, mapping ratio, genus rank (± HOPS criteria)

## 5. Report generation

- Per-sample / per-pathogen PDFs; cohort Excel summaries; detection-score heatmaps; optional multi-QC dashboard

---

# Part B — Detailed steps and parameters

Values below are those implemented in the PIGSTI Snakefile and helper scripts. Which options were used in this study is summarised narratively in **Part C**.

---

## B1. Pre-processing

### Adapter trimming — AdapterRemoval v2 (paired-end)

```text
AdapterRemoval --collapse --minadapteroverlap 1 \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --minlength 30 --gzip --trimns --trimqualities
```

Overlapping PE mates are collapsed into a single consensus read (`--collapse`). `--minadapteroverlap 1` requires only 1 bp of adapter match to trigger clipping (sensitive). Reads shorter than 30 bp after trimming/collapse are discarded. `--trimns` and `--trimqualities` remove terminal Ns and low-quality bases (AdapterRemoval default Phred floor = 2 unless `--minquality` is set).

### Adapter trimming — cutadapt (single-end)

```text
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 -j 6 \
  -o <sample>.collapsed.gz <SE.fastq.gz>
```

Minimum overlap 1 (`-O 1`) and minimum length 30 (`-m 30`).

### Host species identification — FastQ Screen + Bowtie2

Adapter-trimmed collapsed reads are screened against the configured host genomes (from `bwa_indices` or `bowtie2_indices`), with FastQ Screen forced to Bowtie2 (`BOWTIE2` / `bowtie2` in the auto-written conf). Human is recorded for contamination warnings but excluded from best-host selection by default (`scripts/parse_fastq_screen.py`).

**Pass 1.** FastQ Screen’s default subset (~100 000 reads) is aligned to each database. For each genome the report includes, among other fields, `#One_hit_one_genome` (reads that hit exactly one genome once) and `%Unmapped` (fraction of the screened subset that did not map to that genome). Best-host selection uses **both** signals (`resolve_best_species`):

1. Identify the genome with the highest `#One_hit_one_genome`.  
2. Identify the genome with the lowest `%Unmapped`.  
3. If these two genomes agree, that species is written to `*_best_species.txt`.  
4. If they disagree, the genome with the **lowest %Unmapped** is selected (and a contamination warning is logged).  
5. An additional warning is raised if the top one-hit genome is only weakly ahead of the second-best (&lt; 1.5× one-hit ratio).

Nuclear and mtDNA references for mapping are taken from this best-species call.

**Pass 2 (optional full-dataset rescreen).** If `fastq_screen_full_dataset_rescreen` is true (pipeline default) and the top genome’s `#One_hit_one_genome` count is below `fastq_screen_full_dataset_min_one_hit` (default **50**), FastQ Screen is re-run with `--subset 0` (all reads). This catches low-endogenous or shallow libraries where the 100k subset is too noisy for a reliable species call. Best-species resolution again uses one-hit **and** `%Unmapped` as above. Threads: **4**.

---

## B2. Host / mtDNA mapping and post-processing

After FastQ Screen best-host identification, adapter-trimmed collapsed reads are aligned to that host’s nuclear and mtDNA references. Aligner choice: `host_aligner: bwa` or `bowtie2`. Host/mtDNA analysis can be skipped with `pathogen_screening_only: true` (FastQ Screen and unmapped generation for Kraken are retained).

### BWA *aln* / *samse*

```text
bwa aln -l 1024 -n 0.01 -o 2 -t 6 <index> <reads.fq.gz> > out.sai
bwa samse -r '@RG\tID:<sample>_host\tSM:<sample>\tPL:ILLUMINA' \
  <index> out.sai <reads.fq.gz> | samtools view -@ 6 -F 4 -b -o mapped.bam
```

These are classical ancient-DNA *bwa aln* settings (Schubert *et al.*, 2012): `-l 1024` effectively disables seeding so damaged terminal bases are not excluded from alignment; `-n 0.01` relaxes the allowed edit distance (~1% of read length); `-o 2` allows up to two gap opens. Mapped reads only are kept (`samtools view -F 4`). The same *aln* flags are used for mtDNA and for pathogen BWA mapping when `pathogen_aligner: bwa`. With `pathogen_mapping_mode: super_careful`, host-unmapped reads are also exported (`-f 4`) for pathogen mapping; with `default`, pathogen mapping later uses chimera-depleted unmapped reads.

### Bowtie2

```text
bowtie2 --end-to-end --sensitive -x <prefix> -U <reads> -p 6 …
```

End-to-end sensitive mode is a common Bowtie2 setting for aDNA (full-read alignment without local soft-clipping; sensitive preset). Comparative aDNA benchmarking has supported end-to-end Bowtie2 modes, including the sensitive preset, as competitive alternatives to BWA for short damaged reads (e.g. Poullet & Orlando, 2020).

### Post-mapping filters and QC

After mapping, reads are sorted (`samtools sort`), alignments with mapping quality below 30 are discarded (`samtools view -q 30`), and duplicates are removed with `samtools markdup -r` (or Picard MarkDuplicates if configured). BAMs are then soft-clipped for the first/last **4** bases with `softclip_mod.py` and written as CRAM. Coverage and mapping QC use Qualimap2; DamageProfiler is run to profile terminal substitutions. Endogenous DNA % = (deduplicated host-mapped reads / raw reads) × 100. Optional residual sexing is available on the merged host BAM for Cow, Goat, Sheep and Dog when `enable_sexing: true`.

---

## B3. Metagenomics screening

### Dereplication — PRINSEQ++

```text
prinseq++ -fastq <collapsed.gz> -derep 14 -VERBOSE 2 -threads 6 \
  -out_good <sample>-passed.fq -out_bad <sample>-bad.fq
```

`-derep 14` removes exact duplicate / exact-sequence reads (metagenomics branch only; not applied to the host-mapping read stream).

### Mammalian host-read removal

```text
bowtie2 -x <multi_host_chimera_index> -U <prinseq_passed.fq.gz> -p 6 | \
samtools view -Sb - -f4 > unaligned.bam
```

This step uses Bowtie2 **default** presets (not `--end-to-end --sensitive`). `-f4` keeps unmapped reads only (depleted of mammal/chimera-like hits). Per-PCR unaligned FASTQs are produced with `samtools bam2fq | pigz`, then PCR replicates are pooled with `cat` into one biological-sample FASTQ for KrakenUniq / optional HOPS / decOM.

### KrakenUniq

```text
krakenuniq --db <DB> --fastq-input <pooled_unaligned.fq.gz> \
  --threads 8 --output <out> --report-file <report> \
  --gzip-compressed --only-classified-out
```

Optional abundance heatmaps keep taxa with ≥ **1000** KrakenUniq reads and plot the top **25** (absolute and normalised).

### Optional HOPS / MALT / MaltExtract

When `enable_hops: true`, MALT runs in BlastN semi-global mode with percent identity / minimum percent identity **90.0**, `topMalt 1`, `mq 100`, then MaltExtract with `filter=def_anc`, `topMaltEx=0.01`, `minComp=0.6`, `minPIdent=85` (plus memory/thread settings from `hops_parameters`).

### Optional decOM

When `enable_decom: true`, pooled unaligned reads are compared to a pre-built k-mer source matrix (`decOM_sources`) for microbial source tracking.

---

## B4. Pathogen detection, mapping and authentication

### Candidate selection

A pathogen spreadsheet maps KrakenUniq names → reference indices (± HOPS names).

- **`use_evalue_for_detection: true`** — Guellil *et al.* E-value \(E = (K/R) \times C\), where *K* = unique k-mers assigned to the taxon, *R* = assigned reads, and *C* = estimated k-mer coverage from the KrakenUniq report; pass if **E &gt; 0.001** (`guellil_evalue_threshold`).  
- **`use_evalue_for_detection: false`** — E-score ≥ **5** and/or ≥ **50** KrakenUniq reads.  
- **`map_all_escore_pathogens`** — if true, map all spreadsheet taxa present in the detection table; if false, map only taxa passing thresholds.  
- If HOPS is enabled, mapping targets are the **union** of Kraken/E-value and HOPS-positive taxa.

### Pathogen reference mapping

Same aligner as configured (`pathogen_aligner: bwa` → `-l 1024 -n 0.01 -o 2`; or Bowtie2 `--end-to-end --sensitive`). Read source: chimera-unmapped pooled FASTQ (`pathogen_mapping_mode: default`) or host-unmapped FASTQ (`super_careful`). Post-map: keep mapped, MAPQ ≥ 30, sort, `samtools markdup -r`, Qualimap2, DamageProfiler.

### Authentication metrics

Each sample–pathogen pair is evaluated against the following criteria (thresholds from `pathogen_detection_criteria` / config). The reported **Score** is `criteria_passed / criteria_evaluated`.

| Metric | What it measures | Pass rule |
|--------|------------------|-----------|
| KrakenUniq reads | Taxon abundance (unique-kmer aware) | ≥ 50 (or spreadsheet) |
| Guellil E-value / E-score | Enrichment × complexity style signal | E &gt; 0.001 or E-score ≥ 5 |
| ANI | Average nucleotide identity of mapped reads vs reference | &gt; 96.5% |
| Relative entropy | Evenness of coverage along the genome (100 bp / 1000 bp windows) | ≥ 0.9 bacteria; ≥ 0.7 virus |
| Breadth ratio | Observed breadth / expected breadth given read count and length | ≥ 0.8 |
| 5′ C→T damage | Terminal deamination rate (aDNA signal); optional 3′ G→A | ≥ 0.01 |
| Edit-distance decay quality | How well the mismatch histogram (0, 1, 2, … edits) shows a descending / exponential pattern | ≥ 0.65 |
| Edit-distance (damage split) | Same metric on deaminated vs non-deaminated read subsets (window 5 bp) | no-damage ≥ 0.55; damage − no-damage ≥ 0.10 |
| Mapping ratio | Pathogen-mapped reads / KrakenUniq reads for that taxon | ≥ 0.5 |
| Genus rank | Rank of the target among taxa in its genus by Kraken counts | = 1 |
| HOPS detection / edit / damage | HOPS authentication levels (if HOPS on) | ≥ 2 / 3 / 4 |

Damage vs no-damage edit-distance splitting (`edit_distance_damage_split: true`) follows HOPS-like logic: reads with terminal C→T / G→A in a 5 bp window are scored separately so authentic ancient signal can be distinguished from modern contamination-like edit profiles.

---

## B5. Report generation

Per-sample / per-pathogen PDFs; cohort Excel summaries; detection-score heatmaps; optional multi-QC dashboard.

---

# Part C — This study (paper-style narrative)

In this study, shotgun libraries were adapter- and quality-trimmed with AdapterRemoval v2 for paired-end data (`--collapse`, `--minadapteroverlap 1`, `--minlength 30`, TruSeq adapters, `--trimns --trimqualities`) or cutadapt for single-end data (`-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30`). Host species was identified with FastQ Screen and Bowtie2 against Cow, Dog, Human, Horse, Sheep, Pig, Goat, Deer, Cat and Camel genomes, selecting the best host from `#One_hit_one_genome` and `%Unmapped` (preferring lowest unmapped when the two disagreed), with full-dataset rescreen when one-hit support was below 50 reads. Collapsed reads used for metagenomics were dereplicated with PRINSEQ++ (`-derep 14`).

Host nuclear and mitochondrial genomes were mapped with Bowtie2 in end-to-end sensitive mode (`--end-to-end --sensitive`). Alignments were filtered to mapping quality ≥ 30, deduplicated with SAMtools (`markdup -r`), soft-clipped by four terminal bases before Qualimap, and assessed with Qualimap2 and DamageProfiler. Endogenous content and optional residual sexing (cattle, goat, sheep, dog) were computed from the host BAMs.

For metagenomic screening, PRINSEQ-passed reads were depleted against a multi-host chimera Bowtie2 index, unmapped reads retained, PCR replicates pooled, and classified with KrakenUniq (`--gzip-compressed --only-classified-out`) against the NT_2020_microbes database. HOPS/MALT was not used. Microbial source tracking with decOM was enabled.

Pathogen candidates were selected from KrakenUniq using E-score and read thresholds (E-score ≥ 5; ≥ 50 reads), mapping spreadsheet taxa present in the detection table. Pathogen references were aligned with Bowtie2 (`--end-to-end --sensitive`) using chimera-unmapped pooled reads, then filtered (MAPQ ≥ 30), deduplicated, and profiled with Qualimap2 and DamageProfiler. Authentication required ANI &gt; 96.5%, relative entropy ≥ 0.9 (bacteria) or ≥ 0.7 (viruses), breadth ratio ≥ 0.8, 5′ C→T ≥ 0.01, edit-distance decay quality ≥ 0.65 (with damage/no-damage split; window 5 bp), mapping ratio ≥ 0.5 and genus rank = 1. Detection confidence was reported as criteria passed over criteria evaluated.

Analyses were run with Snakemake and conda environments (record pinned software versions before submission).

---

## Notes for editing

Confirm PE vs SE for the published cohort; paste pinned tool versions; convert Part A citations to the journal’s reference style.
