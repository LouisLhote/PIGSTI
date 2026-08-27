# PIGSTI pathogen authentication metrics

Definitions used by PIGSTI for **detection** (KrakenUniq / E-value ± HOPS) and **post-mapping authentication**. Each mapped candidate receives a score **out of 10**, or **out of 13** if HOPS is enabled.

These definitions match the PIGSTI software paper (F1000 Research, in preparation).

---

## Scoring overview

| Block | Criteria | Points |
|-------|----------|--------|
| Screening | KrakenUniq clade read count; Guellil E-value | 2 |
| Mapping authentication | Relative entropy; edit-distance decay (damaged); edit-distance decay (non-damaged); postmortem damage; ANI; breadth ratio; mapping ratio; genus rank | 8 |
| Optional HOPS | MaltExtract edit-distance decline; terminal aDNA damage; damaged reads among edit distance 0 | +3 |
| **Total** | | **10** (or **13** with HOPS) |

Default pass thresholds are listed with each metric. Per-pathogen overrides for reads / E-value can be set in `Pathogen_spreadsheet.csv`.

---

## 1. KrakenUniq read count (`Krakenuniq_reads`)

**Definition.** Clade-level read count from the KrakenUniq report (`reads`): reads assigned to a taxon **plus all descendant taxa**. This is distinct from `taxReads`, which counts only reads assigned to the terminal node.

**Interpretation.** A minimum clade read count filters weak or sporadic k-mer hits. Default: **≥ 50**. More stringent or permissive values can be set per pathogen (e.g. higher for opportunistic commensals such as *Escherichia coli*; lower for poxviruses).

**Source.** Breitwieser, F. P., Baker, D. N. & Salzberg, S. L. KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. *Genome Biology* 19, 198 (2018).

---

## 2. E-value (`Guellil_et_al_Evalue`)

**Definition.** Following Guellil and colleagues:

$$E = \frac{K}{R} \times C$$

where **K** = unique k-mers assigned to the taxon, **R** = reads assigned to the taxon, and **C** = estimated k-mer coverage (genome breadth) from the KrakenUniq report.

**Interpretation.** Elevated *E* reflects reads spread across a larger genomic fraction with proportionally unique k-mer support. Default detection threshold: **E > 0.001**. This is the **default screening filter** in PIGSTI (`use_evalue_for_detection: true`). Set `false` to use KrakenUniq E-score instead.

**Source.** Guellil, M. et al. An invasive *Haemophilus influenzae* serotype b infection in an Anglo-Saxon plague victim. *Genome Biology* 23, 22 (2022).

---

## 3. Average nucleotide identity — ANI (`ANI`)

**Definition.** Mean percentage of identical nucleotides, approximated from mapping statistics:

$$\text{ANI} \approx \left(1 - \frac{\text{mismatches}}{\text{bases mapped}}\right) \times 100$$

Mismatch and mapped-base counts come from `samtools stats`.

**Interpretation.** Reduced ANI may reflect damage-related substitutions, mapping to a distantly related reference, or misclassification of reads from a related species. Default: **> 96.5%**.

**Source.** Li, H. et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics* 25, 2078–2079 (2009).

---

## 4. Relative entropy (`Relative_entropy`)

**Definition.** Evenness of mapped **read start** positions along the reference, adapted from Sikora et al. (2025). Starts are binned into non-overlapping windows of **100 bp** and **1,000 bp**. Shannon entropy of the start distribution is normalized to the maximum under a uniform distribution:

$$H_{\mathrm{rel}} = \frac{H}{H_{\max}} = \frac{-\sum_i p_i \log_2 p_i}{\log_2(N_{\mathrm{windows}})}$$

where \(p_i\) is the proportion of read starts in window *i*.

**Interpretation.** Values approaching 1 indicate coverage consistent with a genuine, evenly distributed pathogen signal. Low values are consistent with localized clustering (contamination or mismapping). Default: **≥ 0.9** for bacteria/archaea; **≥ 0.7** for viruses.

**Source.** Sikora, M. et al. The spatiotemporal distribution of human pathogens in ancient Eurasia. *Nature* 643, 1011–1019 (2025).

---

## 5. Breadth ratio (`Breadth_ratio`)

**Definition.** Observed versus expected genome breadth, adapted from Sikora et al. (2025):

$$\text{Breadth ratio} = \frac{B_{\mathrm{observed}}}{B_{\mathrm{expected}}},\qquad B_{\mathrm{expected}} = 1 - e^{-\mathrm{mean\ depth}}$$

\(B_{\mathrm{expected}}\) is the Poisson expectation for uniform coverage at a given mean depth.

**Interpretation.** Values approaching 1 indicate even mapping; values well below 1 indicate patchy or clustered mapping (e.g. cross-mapping from a related taxon). Default: **≥ 0.8**.

**Source.** Sikora, M. et al. *Nature* 643, 1011–1019 (2025).

---

## 6–7. Edit-distance decay (`Edit_distance_decay_quality`, `Edit_distance_decay_quality_default`)

**Definition.** Edit distance is the number of mismatches between a read and the mapped reference. Reads are split into **damaged** and **non-damaged** subsets based on terminal deamination (5′ C→T and/or 3′ G→A within the terminal **five** aligned bases), following HOPS logic (Hübler et al. 2019).

For each subset, PIGSTI computes a composite **Decay Quality Score** (0–1):

| Component | Weight |
|-----------|--------|
| Monotonicity of read-count decline across mismatch bins | 0.35 |
| Ratio of reads at edit distance 0 vs 1 | 0.25 |
| Proportion of reads at edit distance 0 | 0.15 |
| Penalty if the modal edit distance is not 0 | 0.10 |
| Exponential decay rate | 0.10 |
| R² of the exponential fit | 0.05 |

**Interpretation.** Higher scores indicate a clean descending mismatch profile. Default pass: **≥ 0.65** (damaged) and **≥ 0.55** (non-damaged). This split is always on in the pipeline.

**Source.** Hübler, R. et al. HOPS: automated detection and authentication of pathogen DNA in archaeological remains. *Genome Biology* 20, 280 (2019). Composite score: PIGSTI `scripts/calculate_edit_distance_r2.py`.

---

## 8. Postmortem damage (`Damage_5p_CtoT`)

**Definition.** Frequency of C→T substitutions at the terminal 5′ position of aligned reads, estimated with DamageProfiler.

**Interpretation.** Supports ancient origin of mapped reads. Default: **≥ 0.01**.

**Source.** Neukamm, J., Peltzer, A. & Nieselt, K. DamageProfiler: fast damage pattern calculation for ancient DNA. *Bioinformatics* 37, 3652–3653 (2021).

---

## 9. Mapping ratio (`Read_mapping_ratio`)

**Definition.** Concordance between KrakenUniq classification and reference alignment:

$$\text{Mapping ratio} = \frac{\text{mapped reads (BWA/Bowtie2)}}{\text{KrakenUniq-assigned reads}}$$

**Interpretation.** Default: **≥ 0.5**. Values ≪ 0.5 suggest misclassification, a mismatched reference, or heavy mapping filters.

**Source.** PIGSTI cross-validation metric.

---

## 10. Genus rank (`Genus_ranking`)

**Definition.** Rank of the candidate taxon among all species of the **same genus** in the KrakenUniq report, ordered by descending clade read count.

**Interpretation.** Rank **1** (dominant genus-level hit) is required to award a detection point, reducing ambiguous assignment among closely related taxa.

**Source.** PIGSTI implementation (genus-level disambiguation as used in large-scale ancient pathogen screens; cf. Sikora et al. 2025).

---

## 11–13. HOPS / MaltExtract (optional)

When `enable_hops: true`, PIGSTI additionally evaluates three MaltExtract criteria (Hübler et al. 2019), contributing **up to 3 points** (maximum score 13):

1. Decline of the edit-distance distribution  
2. Characteristic terminal aDNA damage patterns  
3. Proportion of damaged reads among those with edit distance zero  

Default HOPS score floors in config: detection ≥ 2, edit distance ≥ 3, damage ≥ 4 (see `docs/CONFIG.md`).

---

## Default thresholds (software defaults)

| Metric | Threshold | Direction |
|--------|-----------|-----------|
| KrakenUniq clade reads | ≥ 50 | higher better |
| Guellil E-value | > 0.001 | higher better |
| ANI | > 96.5% | higher better |
| Relative entropy | ≥ 0.9 (bacteria/archaea); ≥ 0.7 (virus) | higher better |
| Breadth ratio | ≥ 0.8 | higher better |
| Edit-distance decay (damaged) | ≥ 0.65 | higher better |
| Edit-distance decay (non-damaged) | ≥ 0.55 | higher better |
| 5′ C→T damage | ≥ 0.01 | higher better |
| Mapping ratio | ≥ 0.5 | higher better |
| Genus ranking | = 1 | rank 1 required |

Empirical studies may relax evenness criteria (entropy, breadth) for very low-abundance authentic hits; that is a **study-specific** choice, not the pipeline default.

---

## Key references

1. Breitwieser, F. P., Baker, D. N. & Salzberg, S. L. *Genome Biology* 19, 198 (2018).
2. Guellil, M. et al. *Genome Biology* 23, 22 (2022).
3. Sikora, M. et al. *Nature* 643, 1011–1019 (2025).
4. Hübler, R. et al. *Genome Biology* 20, 280 (2019).
5. Neukamm, J. et al. *Bioinformatics* 37, 3652–3653 (2021).
6. Li, H. et al. *Bioinformatics* 25, 2078–2079 (2009).
