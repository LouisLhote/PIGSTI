configfile: "config/config.yaml"

ruleorder: adapter_removal_pe > adapter_removal_se

global_wildcard_constraints = {
    "krakenuniq_name": r"[A-Za-z0-9_.-]+",
    "sample": r"[A-Za-z0-9_.-]+",
    "ref_name_safe": r"[A-Za-z0-9_]+"
}

import os
import pandas as pd
import re
import csv
from pathlib import Path


# -------------------- Load sample info --------------------
SAMPLES_DICT = {}
with open("config/samples.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample = row["sample"]
        if sample not in SAMPLES_DICT:
            SAMPLES_DICT[sample] = {"r1": [], "r2": []}
        SAMPLES_DICT[sample]["r1"].append(row["r1"])
        if row["r2"] and row["r2"].strip():
            SAMPLES_DICT[sample]["r2"].append(row["r2"])

SAMPLES = list(SAMPLES_DICT.keys())
READ_GROUPS = {}
with open("config/samples.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample = row["sample"]
        rglb = row["RGLB"]
        sequencing_run = row["sequencing_run"]
        rg_id = f"{rglb}_{sequencing_run}"
        read_group = f"@RG\\tID:{rg_id}\\tPL:ILLUMINA\\tLB:{rglb}\\tSM:{sample}"
        READ_GROUPS[sample] = read_group



# Load spreadsheet once globally and strip columns once
spreadsheet_df = pd.read_csv("config/Pathogen_spreadsheet.csv")
spreadsheet_df.columns = spreadsheet_df.columns.str.strip()

def get_sample_ref_pairs():
    """Dynamically generate SAMPLE_REF_PAIRS based on available escore files"""
    pairs = []
    for sample in SAMPLES:
        escore_path = f"results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
        if not os.path.exists(escore_path):
            continue
        escore_df = pd.read_csv(escore_path)
        escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()

        # Only add pathogens from escore that are in master spreadsheet
        for pathogen in escore_df["taxonomy"]:
            if pathogen in spreadsheet_df["Krakenuniq name"].values:
                pairs.append((sample, pathogen))
    return pairs




def get_all_possible_pathogens():
    """Get all possible pathogens from the spreadsheet for all samples"""
    pairs = []
    for sample in SAMPLES:
        for pathogen in spreadsheet_df["Krakenuniq name"].values:
            pairs.append((sample, pathogen))
    return pairs







def get_reference_path(wc):
    row = spreadsheet_df[spreadsheet_df["Krakenuniq name"] == wc.krakenuniq_name]
    if row.empty:
        raise ValueError(f"No reference path found for {wc.krakenuniq_name}")
    return row.iloc[0]["bwa index"]


def get_bwa_targets(wildcards):
    escore_path = f"results/{wildcards.sample}/Escore/pathogen/{wildcards.sample}_pathogen.csv"
    if not os.path.exists(escore_path):
        return []
    escore_df = pd.read_csv(escore_path)
    # Use the globally loaded spreadsheet_df here
    targets = []
    for pathogen in escore_df["name"]:
        row = spreadsheet_df[spreadsheet_df["Krakenuniq name"] == pathogen]
        if not row.empty:
            krakenuniq_name = row.iloc[0]["Krakenuniq name"]
            bam_path = f"results/{wildcards.sample}/bwa_final/q30/{krakenuniq_name}_F4_q30_sort.bam"
            targets.append(bam_path)
    return targets

# rest of your code remains unchanged...



def safe_name(name):
    return name.replace(" ", "_")

def get_sample_ref_pairs_safe():
    """Get SAMPLE_REF_PAIRS_SAFE dynamically"""
    pairs = get_sample_ref_pairs()
    return [(sample, safe_name(pathogen)) for sample, pathogen in pairs]




print("SAMPLE_REF_PAIRS for all samples:")
for s, p in get_sample_ref_pairs():
    print(s, p)


def get_reference_from_safe_name(ref_name_safe):
    # Clean the ref_name_safe to remove any path components
    clean_ref_name = ref_name_safe.split("/")[0]
    pathogen_name = clean_ref_name.replace("_", " ")
    row = spreadsheet_df.loc[spreadsheet_df["Krakenuniq name"] == pathogen_name]
    if row.empty:
        raise ValueError(f"[ERROR] No bwa index found for pathogen name: '{pathogen_name}'")
    return row.iloc[0]["bwa index"]

def get_zipped_pairs():
    """Get zipped_pairs dynamically"""
    pairs = get_sample_ref_pairs()
    return [(s, safe_name(r)) for s, r in pairs]

def get_pathogen_bwa_targets(wildcards):
    """Get pathogen BWA targets after checkpoint completion"""
    try:
        # After checkpoint completion, get the actual targets
        targets = []
        for sample, pathogen in get_sample_ref_pairs():
            pathogen_safe = safe_name(pathogen)
            targets.append(f"results/{sample}/bwa_pathogen/{sample}_{pathogen_safe}.dedup.bam")
        return targets
    except:
        return []


###--------------------------wrappers-----------------------------------------------------


###-----------------------------------------------rules-------------------------------------------------------------------------------
rule all:
    input:
        expand("results/{sample}/adapter_removal/{sample}.collapsed.gz", sample=SAMPLES),
        expand("results/{sample}/krakenuniq/kraken-report.txt", sample=SAMPLES),
        expand("results/{sample}/fastq_screen/{sample}.collapsed_screen.html", sample=SAMPLES),
        expand("results/{sample}/fastq_screen/{sample}_best_species.txt", sample=SAMPLES),
        expand("results/{sample}/krona/{sample}.html", sample=SAMPLES),

        # Host alignments
        expand("results/{sample}/bwa_host/{sample}_F4.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}_F4_q30.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}_F4_q30.sorted.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}.dedup.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}.dedup.metrics.txt", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}_F4_q30.bam.bai", sample=SAMPLES),
        expand("results/{sample}/bwa_host/{sample}.dedup.bam.bai", sample=SAMPLES),

        # Host QC and damage
        expand("results/{sample}/qualimap/genome_results.txt", sample=SAMPLES),
        expand("results/{sample}/qualimap_mtdna/genome_results.txt", sample=SAMPLES),
        expand("results/{sample}/damageprofiler_host", sample=SAMPLES),

        # mtDNA alignments
        expand("results/{sample}/bwa_mtdna/{sample}_F4.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_mtdna/{sample}_F4_q30.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_mtdna/{sample}.dedup.bam", sample=SAMPLES),
        expand("results/{sample}/bwa_mtdna/{sample}.dedup.metrics.txt", sample=SAMPLES),
        expand("results/{sample}/bwa_mtdna/{sample}_F4_q30.bam.bai", sample=SAMPLES),
        expand("results/{sample}/bwa_mtdna/{sample}.dedup.bam.bai", sample=SAMPLES),

        # mtDNA QC and damage
        expand("results/{sample}/damageprofiler_mtdna", sample=SAMPLES),

        # E-Score outputs
        expand("results/{sample}/Escore/genus/{sample}_genus.csv", sample=SAMPLES),
        expand("results/{sample}/Escore/species/{sample}_species.csv", sample=SAMPLES),
        expand("results/{sample}/Escore/pathogen/{sample}_pathogen.csv", sample=SAMPLES),

        # Comparison outputs
        expand("results/comparison/{sample}_comparison.tsv", sample=SAMPLES),
        expand("results/comparison/{sample}_heatmap.html", sample=SAMPLES),

        #decom
        "results/p_sink.txt",
        expand("results/p_keys/{sample}.fof", sample=SAMPLES),
        "results/decOM",

        # Shared/global files
        expand("results/{sample}/krakenuniq/output.txt", sample=SAMPLES),
        "results/hops/maltExtract/heatmap_overview_Wevid.tsv",
        "lists/krakenuniq_pathogen_list.txt",
        "lists/hops_pathogen_list.txt",
        "config/config_hops_custom.txt",
        "results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix_absolute.csv",
        "results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix_normalized.csv",
        "results/KRAKENUNIQ_ABUNDANCE_MATRIX/heatmap_absolute.pdf",
        "results/KRAKENUNIQ_ABUNDANCE_MATRIX/heatmap_normalized.pdf",

        # All BWA pathogen targets - will be empty initially, populated after escore
        "results/pathogen_targets.txt",



        # Final summary reports

        "results/final/host_mtdna_summary_all_samples.xlsx",
        "results/pathogen_bwa_complete.txt",
        "results/pathogen_detection/detection_scores_heatmap.pdf",
        "results/pathogen_detection/detection_scores_matrix.csv"

# Generate pathogen targets after escore completes
checkpoint generate_pathogen_targets:
    input:
        escore_files = expand("results/{sample}/Escore/pathogen/{sample}_pathogen.csv", sample=SAMPLES)
    output:
        targets_file = "results/pathogen_targets.txt"
    run:
        targets = []
        for sample in SAMPLES:
            escore_path = f"results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
            if os.path.exists(escore_path):
                escore_df = pd.read_csv(escore_path)
                escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()
                
                # Only add pathogens from escore that are in master spreadsheet
                for pathogen in escore_df["taxonomy"]:
                    if pathogen in spreadsheet_df["Krakenuniq name"].values:
                        pathogen_safe = safe_name(pathogen)
                        targets.append(f"results/{sample}/bwa_pathogen/{sample}_{pathogen_safe}.dedup.bam")
        
        with open(output.targets_file, "w") as f:
            for target in targets:
                f.write(target + "\n")

def expand_downstream_targets(wildcards):
    # Get checkpoint output path
    ckpt_output = checkpoints.generate_pathogen_targets.get(**wildcards).output.targets_file
    targets = []
    with open(ckpt_output) as f:
        for bam in f:
            bam = bam.strip()
            if not bam:
                continue
            sample, pathogen_file = Path(bam).stem.split("_", 1)
            pathogen = pathogen_file.replace(".dedup", "")
            targets.extend([
                bam,
                f"results/{sample}/bwa_pathogen/qualimap_{pathogen}",
                f"results/{sample}/bwa_pathogen/damageprofiler_{pathogen}",
                f"results/{sample}/bwa_pathogen/adnaplotter_{pathogen}.pdf",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.ani.txt",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.depth.txt",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.breadth.txt",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.expected_breadth.txt",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.breadth_ratio.txt",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.entropy_plot.png",
                f"results/{sample}/bwa_pathogen/{sample}_{pathogen}.mean_entropy.txt"
            ])
    return targets

checkpoint generate_pathogen_targets:
    input:
        escore_files = expand("results/{sample}/Escore/pathogen/{sample}_pathogen.csv", sample=SAMPLES)
    output:
        targets_file = "results/pathogen_targets.txt"
    run:
        targets = []
        for sample in SAMPLES:
            escore_path = f"results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
            if os.path.exists(escore_path):
                escore_df = pd.read_csv(escore_path)
                escore_df["taxonomy"] = escore_df["taxonomy"].str.strip()
                for pathogen in escore_df["taxonomy"]:
                    if pathogen in spreadsheet_df["Krakenuniq name"].values:
                        pathogen_safe = safe_name(pathogen)
                        targets.append(f"results/{sample}/bwa_pathogen/{sample}_{pathogen_safe}.dedup.bam")
        with open(output.targets_file, "w") as f:
            for target in targets:
                f.write(target + "\n")

rule pathogen_bwa_targets:
    input:
        expand_downstream_targets
    output:
        touch("results/pathogen_bwa_complete.txt")
    shell:
        "echo 'Pathogen BWA alignment and downstream analysis completed for all samples' > {output}"

#--------------------list and configs files -------------------------------------------

rule generate_pathogen_lists:
    input:
        spreadsheet = "config/Pathogen_spreadsheet.csv"
    output:
        kraken_list = "lists/krakenuniq_pathogen_list.txt",
        hops_list = "lists/hops_pathogen_list.txt"
    conda:
        "workflow/envs/python.yaml"  # or your python env
    shell:
        """
        python scripts/generate_pathogen_lists.py {input.spreadsheet}
        """

rule create_hops_config:
    input:
        original_config = "config/config_hops2.0.txt",
        hops_list = "lists/hops_pathogen_list.txt"
    output:
        new_config = "config/config_hops_custom.txt"
    conda:
        "workflow/envs/python.yaml"
    shell:
        """
        python scripts/create_hops_config.py {input.original_config} {output.new_config} {input.hops_list}
        """

rule select_best_pathogens:
    input:
        spreadsheet="config/Pathogen_spreadsheet.csv",
        escore="results/{sample}/sample_pathogen.csv",
        taxdump_dir="ncbi_taxdump/nodes.dmp"
    output:
        selected="results/{sample}/selected_pathogens.txt"
    params:
        sample="{sample}"
    shell:
        "python scripts/group_by_genus.py --sample {params.sample}"
#-------------------bwa pathogen mapping----------------------------------------------

rule bwa_aln:
    input:
        reads = "results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz",
        pathogen ="results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
    output:
        sai = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.sai"
    conda:
        "workflow/envs/bwa.yaml"
    threads: 6
    params:
        reference = lambda wc: get_reference_from_safe_name(wc.ref_name_safe)
    shell:
        "bwa aln -l 1024 -n 0.01 -o 2 -t {threads} {params.reference} {input.reads} > {output.sai}"

rule bwa_samse:
    input:
        sai = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.sai",
        reads = "results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz"
    output:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4.bam"
    conda:
        "workflow/envs/bwa.yaml"
    params:
        reference =  lambda wc: get_reference_from_safe_name(wc.ref_name_safe),
        read_group = lambda wc: READ_GROUPS[wc.sample]
    shell:
        """
        bwa samse -r '{params.read_group}' {params.reference} {input.sai} {input.reads} | \
        samtools view -F 4 -Sb - > {output.bam}
        """

rule bwa_sort:
    input: "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4.bam"
    output: "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4.sorted.bam"
    conda: "workflow/envs/bwa.yaml"
    shell: "samtools sort {input} -o {output}"

rule bwa_q30:
    input: "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4.sorted.bam"
    output: "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4_q30.bam"
    conda: "workflow/envs/bwa.yaml"
    shell: "samtools view -q 30 -o {output} {input}"

rule bwa_q30_sort:
    input:
        "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4_q30.bam"
    output:
        "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4_q30_sort.bam"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule mark_duplicates:
    input:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}_F4_q30_sort.bam"
    output:
        dedup = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam",
        metrics = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.metrics.txt"
    conda: "workflow/envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.dedup} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT
        """

rule index_dedup_bam:
    input:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam"
    output:
        bai = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam.bai"
    conda: "workflow/envs/samtools.yaml"
    shell: "samtools index {input.bam}"

rule qualimap_bamqc_bwa:
    input:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam"
    output:
        report = directory("results/{sample}/bwa_pathogen/qualimap_{ref_name_safe}")
    log:
        "logs/qualimap/{sample}_{ref_name_safe}.log"
    conda: "workflow/envs/qualimap.yaml"
    shell:
        """
        qualimap bamqc \
            -bam {input.bam} \
            -outdir {output.report} \
            -outformat html > {log} 2>&1
        """
def get_ref_path(wc):
    ref_name_clean = wc.ref_name_safe.split("/")[0]  # strip anything after slash
    matches = spreadsheet_df.loc[
        spreadsheet_df["Krakenuniq name"] == ref_name_clean.replace("_", " "),
        "bwa index"
    ].values
    if len(matches) == 0:
        raise ValueError(f"No bwa index found for reference: {ref_name_clean}")
    return matches[0]


rule damageprofiler:
    input:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam",
        ref = get_ref_path
    output:
        directory("results/{sample}/bwa_pathogen/damageprofiler_{ref_name_safe}")
    conda: "workflow/envs/damageprofiler.yaml"
    shell:
        "damageprofiler -i {input.bam} -o {output} -r {input.ref}"
rule adna_bamplotter:
    input:
        bam = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam",
        bai = "results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam.bai",
        dpdir = "results/{sample}/bwa_pathogen/damageprofiler_{ref_name_safe}/"
    output:
        pdf = "results/{sample}/bwa_pathogen/adnaplotter_{ref_name_safe}.pdf"
    conda: "workflow/envs/adna_plotter.yaml"
    shell:
        """
        python scripts/aDNA-BAMPlotter.py \
            -b {input.bam} \
            -d {input.dpdir}/misincorporation.txt \
            -o {output.pdf}
        """
rule Compute_ANI:
    input:
        bam="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam"
    output:
        ani="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.ani.txt"
    log:
        "logs/ANI/{sample}_{ref_name_safe}.log"
    threads: 1
    conda:
        "workflow/envs/samtools.yaml"  # ensure samtools is included
    message:
        "Compute_ANI: Calculating ANI for {input.bam}"
    shell:
        """
        samtools stats {input.bam} 2>> {log} | \
        awk '/^SN/ && /mismatches:/ {{mis=$3}} /^SN/ && /bases mapped:/ {{map=$4}} \
        END {{if (map > 0) printf("ANI ≈ %.2f%%\\n", (1 - mis/map)*100); else print "No mapped bases!"}}' \
        > {output.ani}
        """

rule MappingStats:
    input:
        bam="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam"
    output:
        depth="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.depth.txt",
        breadth="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.breadth.txt",
        expected_breadth="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.expected_breadth.txt",
        ratio="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.breadth_ratio.txt"
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_breadth_stats.py"


rule EntropyProfile:
    input:
        bam="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.dedup.bam"
    output:
        plot="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.entropy_plot.png",
        summary="results/{sample}/bwa_pathogen/{sample}_{ref_name_safe}.mean_entropy.txt"
    conda:
        "workflow/envs/qc.yaml"
    script:
        "scripts/calculate_entropy_profile.py"


##plus---------------




# -------------------- Pathogen Tracker Preprocessing -----------------------

rule adapter_removal_pe:
    input:
        r1 = lambda wc: SAMPLES_DICT[wc.sample]["r1"][0],
        r2 = lambda wc: SAMPLES_DICT[wc.sample]["r2"][0]
    output:
        collapsed = "results/{sample}/adapter_removal/{sample}.collapsed.gz"
    conda:
        "workflow/envs/adapterremoval.yaml"
    threads: 6
    shell:
        """
        AdapterRemoval --file1 {input.r1} --file2 {input.r2} \
        --basename results/{wildcards.sample}/adapter_removal/{wildcards.sample} \
        --threads {threads} --collapse --minadapteroverlap 1 \
        --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --minlength 30 --gzip --trimns --trimqualities
        """

rule adapter_removal_se:
    input:
        r1 = lambda wc: SAMPLES_DICT[wc.sample]["r1"][0]
    output:
        collapsed = "results/{sample}/adapter_removal/{sample}.collapsed.gz"
    conda:
        "workflow/envs/adapterremoval.yaml"
    threads: 6
    shell:
        """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                 -O 1 -m 30 \
                 {input.r1} -o {output.collapsed} -j {threads}
        """


rule fastq_screen:
    input:
        collapsed="results/{sample}/adapter_removal/{sample}.collapsed.gz"
    output:
        html="results/{sample}/fastq_screen/{sample}.collapsed_screen.html",
        txt="results/{sample}/fastq_screen/{sample}.collapsed_screen.txt"
    threads: 4
    shell:
        """
        /raid_md0/Software/FastQ-Screen-0.15.2/fastq_screen \
        -conf /home/kdaly/fastq_screen.conf \
        {input.collapsed} \
        --outdir results/{wildcards.sample}/fastq_screen
        """

rule prinseq:
    input:
        collapsed="results/{sample}/adapter_removal/{sample}.collapsed.gz"
    output:
        passed="results/{sample}/prinseq/{sample}-passed.fq.gz"
    threads: 6
    conda:
        "workflow/envs/prinseq.yaml"
    shell:
        """
        prinseq++ -fastq {input.collapsed} -derep 14 \
        -out_good results/{wildcards.sample}/prinseq/{wildcards.sample}-passed.fq \
        -VERBOSE 2 -threads {threads}
        pigz results/{wildcards.sample}/prinseq/{wildcards.sample}-passed.fq
        """

rule bowtie2_unaligned:
    input:
        passed="results/{sample}/prinseq/{sample}-passed.fq.gz"
    output:
        bam="results/{sample}/bowtie2/{sample}_unaligned.bam"
    threads: 6
    conda:
        "workflow/envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x {config[host_index]} -U {input.passed} -p {threads} | \
        samtools view -Sb - -f4 > {output.bam}
        """

rule bam_to_fastq:
    input:
        bam="results/{sample}/bowtie2/{sample}_unaligned.bam"
    output:
        fastq="results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz"
    threads: 2
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools bam2fq {input.bam} -@{threads} | pigz - > {output.fastq}"

rule krakenuniq:
    input:
        fastq="results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz"
    output:
        report="results/{sample}/krakenuniq/kraken-report.txt",
        output="results/{sample}/krakenuniq/output.txt"
    threads: 8
    conda:
        "workflow/envs/krakenuniq.yaml"
    shell:
        """
        krakenuniq --preload --db {config[kraken_db]} --fastq-input {input.fastq} \
        --threads {threads} --output {output.output} --report-file {output.report} \
        --gzip-compressed --only-classified-out
        """
rule escore:
    input:
        report="results/{sample}/krakenuniq/kraken-report.txt",
        config="config/config.yaml"
    output:
        genus="results/{sample}/Escore/genus/{sample}_genus.csv",
        species="results/{sample}/Escore/species/{sample}_species.csv",
        pathogen="results/{sample}/Escore/pathogen/{sample}_pathogen.csv"
    conda:
        "workflow/envs/escore.yaml"
    params:
        script="scripts/dExp_Escore.py"
    shell:
        """
        python {params.script} {input.report} {output.genus} {output.species} {output.pathogen} {input.config}
        """


rule parse_fastq_screen:
    input:
        "results/{sample}/fastq_screen/{sample}.collapsed_screen.txt"
    output:
        "results/{sample}/fastq_screen/{sample}_best_species.txt"
    params:
        exclude_human=True
    script:
        "scripts/parse_fastq_screen.py"
    

rule hops:
    input:
        genus=expand("results/{sample}/Escore/genus/{sample}_genus.csv", sample=SAMPLES),
        species=expand("results/{sample}/Escore/species/{sample}_species.csv", sample=SAMPLES),
        pathogen=expand("results/{sample}/Escore/pathogen/{sample}_pathogen.csv", sample=SAMPLES),
        fq=expand("results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz", sample=SAMPLES)
    output:
        heatmap="results/hops/maltExtract/heatmap_overview_Wevid.tsv"
    params:
        config="config/config_hops_custom.txt"
    conda:"workflow/envs/hops.yaml"
    threads: 15
        
    shell:
        """
        hops -Xmx800G -input {input.fq} -output results/hops -m full -c {params.config}
        """

rule compare_pathogens:
    input:
        escore="results/{sample}/Escore/pathogen/{sample}_pathogen.csv",
        hops="results/hops/maltExtract/heatmap_overview_Wevid.tsv",
        spreadsheet="config/Pathogen_spreadsheet.csv"
    output:
        tsv="results/comparison/{sample}_comparison.tsv",
        html="results/comparison/{sample}_heatmap.html"
    params:
        sample="{sample}"
    conda:
        "workflow/envs/compare_pathogens.yaml"
    shell:
        "python scripts/compare_pathogens.py {input.escore} {input.hops} {params.sample} {input.spreadsheet}"


##------------------------------------------Microbiome--------------------#######################

rule gather_sinks:
    input:
        fastqs=expand("results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz", sample=SAMPLES)
    output:
        sink_txt="results/p_sink.txt",
        fof_files=expand("results/p_keys/{sample}.fof", sample=SAMPLES)
    run:
        import os
        os.makedirs("results/p_keys", exist_ok=True)
        with open(output.sink_txt, "w") as f_sink:
            for sample in SAMPLES:
                f_sink.write(sample + "\n")
                with open(f"results/p_keys/{sample}.fof", "w") as f_fof:
                    f_fof.write(f"{sample}: results/{sample}/unaligned_fastq/{sample}_unaligned.fastq.gz\n")


rule decom_run:
    input:
        p_sink="results/p_sink.txt",
        fof_files=expand("results/p_keys/{sample}.fof", sample=SAMPLES)
    output:
        directory("results/decOM")
    params:
        p_sources=config["decOM_sources"],
        memory="64GB",
        threads=8
    conda:
        "workflow/envs/decom.yaml"
    log:
        "logs/decom_run.log"
    shell:
        """
    rm -rf {output}
    decOM -p_sinks {input.p_sink} \
          -p_sources {params.p_sources} \
          -p_keys results/p_keys/ \
          -mem {params.memory} \
          -t {params.threads} \
          -o {output} || true

        """


rule update_krona_taxonomy:
    conda:
        "workflow/envs/krona.yaml"
    output:
        touch("taxonomy/.updated")  # a dummy file to mark completion
    shell:
        """
        # Make sure taxonomy dir exists
        mkdir -p taxonomy
        # Run taxonomy update script; adjust path if needed
        ktUpdateTaxonomy.sh
        # Touch the marker file to signal completion
        touch {output}
        """

rule krakenuniq_to_krona:
    input:
        kraken_report="results/{sample}/krakenuniq/kraken-report.txt",
        taxonomy_update="taxonomy/.updated"  # dependency on taxonomy update
    output:
        "results/{sample}/krona/{sample}.html"
    conda:
        "workflow/envs/krona.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/krona
        ktImportTaxonomy -t 5 -m 3 -o {output} {input.kraken_report}
        """


rule krakenuniq_abundance_matrix:
    input:
        reports = expand("results/{sample}/krakenuniq/kraken-report.txt", sample=SAMPLES),
        script1 = "scripts/krakenuniq_abundance_matrix.R",
        script2 = "scripts/plot_krakenuniq_abundance_matrix.R"
    output:
        matrix = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix_absolute.csv",
        matrix_norm = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix_normalized.csv",
        abs_plot = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/heatmap_absolute.pdf",
        norm_plot = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/heatmap_normalized.pdf"
    conda:
        "workflow/envs/r.yaml"
    shell:
        """
        Rscript {input.script1} results results/KRAKENUNIQ_ABUNDANCE_MATRIX 1000 25
        Rscript {input.script2} results/KRAKENUNIQ_ABUNDANCE_MATRIX {output.abs_plot} {output.norm_plot}
        """

# -------------------- Host bwa aln Mapping -----------------------------
rule bwa_aln_host:
    input:
        reads = "results/{sample}/adapter_removal/{sample}.collapsed.gz",
        species = "results/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        bam = "results/{sample}/bwa_host/{sample}_F4.bam"
    threads: 6
    conda:
        "workflow/envs/bwa.yaml"
    script:
        "scripts/bwa_aln_host.py"


rule sort_bam_initial_bwa_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}_F4.bam"
    output:
        bam = "results/{sample}/bwa_host/{sample}_F4.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"

rule filter_q30_bwa_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}_F4.sorted.bam"
    output:
        bam = "results/{sample}/bwa_host/{sample}_F4_q30.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools view -q 30 -@ {threads} -b {input.bam} -o {output.bam}"

rule index_q30_bam_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}_F4_q30.bam"
    output:
        bai = "results/{sample}/bwa_host/{sample}_F4_q30.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule sort_q30_bam_bwa_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}_F4_q30.bam"
    output:
        bam = "results/{sample}/bwa_host/{sample}_F4_q30.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"

rule mark_duplicates_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}_F4_q30.sorted.bam"
    output:
        bam = "results/{sample}/bwa_host/{sample}.dedup.bam",
        metrics = "results/{sample}/bwa_host/{sample}.dedup.metrics.txt"
    threads: 4
    conda: "workflow/envs/picard.yaml"
    shell:
        """
        java -jar /raid_md0/Software/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT
        """
rule index_dedup_bam_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}.dedup.bam"
    output:
        bai = "results/{sample}/bwa_host/{sample}.dedup.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule qualimap_bamqc_bwa_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}.dedup.bam"
    output:
        txt = "results/{sample}/qualimap/genome_results.txt"
    log:
        "logs/qualimap/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/qualimap
        qualimap --java-mem-size=9G \
                 bamqc \
                 -bam {input.bam} \
                 -outdir results/{wildcards.sample}/qualimap \
                 > {log} 2>&1
        """

rule damage_profiler_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}.dedup.bam",
        species_file = "results/{sample}/fastq_screen/{sample}_best_species.txt",
        config = "config/config.yaml"
    output:
        dir = directory("results/{sample}/damageprofiler_host/")
    log:
        "logs/damageprofiler/{sample}.log"
    conda:
        "workflow/envs/damageprofiler.yaml"
    shell:
        """
        mkdir -p {output.dir}

        species=$(cat {input.species_file})

        ref=$(python -c "import yaml; cfg = yaml.safe_load(open('{input.config}')); print(cfg['bwa_indices'].get('${{species}}', ''))")
        if [ -z "$ref" ]; then
            echo "Reference for species '$species' not found in {input.config}" >&2
            exit 1
        fi

        damageprofiler -Xmx12G \
            -i {input.bam} \
            -o {output.dir} \
            -r "$ref" > {log} 2>&1
        """

rule softclip_cram_host:
    input:
        bam = "results/{sample}/bwa_host/{sample}.dedup.bam",
        ref = lambda wc: config["bwa_indices"][
            open(f"results/{wc.sample}/fastq_screen/{wc.sample}_best_species.txt").read().strip()
        ]
    output:
        cram = "results/{sample}/bwa_host/{sample}.dedup_q30_softclipped.cram"
    threads: 4
    conda: "workflow/envs/samtools.yaml"  # Or a custom env with Python2 + samtools if needed
    shell:
        """
        samtools view -h {input.bam} | \
        python2 /home/kdaly/programs/scripts_for_goat_project/softclip_mod.py - 4 | \
        samtools view -@ {threads} -T {input.ref} -O CRAM -o {output.cram}
        """


####-----------------------------------------------------########mtDNA mapping#####----------------------------------------------------------------------------
rule bwa_aln_mtdna:
    input:
        reads = "results/{sample}/adapter_removal/{sample}.collapsed.gz",
        species = "results/{sample}/fastq_screen/{sample}_best_species.txt"
    output:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4.bam"
    threads: 6
    conda:
        "workflow/envs/bwa.yaml"
    script:
        "scripts/bwa_aln_mtdna.py"


rule sort_bam_initial_bwa_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4.bam"
    output:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"


rule filter_q30_bwa_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4.sorted.bam"
    output:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4_q30.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools view -q 30 -@ {threads} -b {input.bam} -o {output.bam}"

rule index_q30_bam_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4_q30.bam"
    output:
        bai = "results/{sample}/bwa_mtdna/{sample}_F4_q30.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule sort_q30_bam_bwa_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4_q30.bam"
    output:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4_q30.sorted.bam"
    threads: 4
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam}"


rule mark_duplicates_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}_F4_q30.sorted.bam"
    output:
        bam = "results/{sample}/bwa_mtdna/{sample}.dedup.bam",
        metrics = "results/{sample}/bwa_mtdna/{sample}.dedup.metrics.txt"
    threads: 4
    conda: "workflow/envs/picard.yaml"
    shell:
        """
        java -jar /raid_md0/Software/picard.jar MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT
        """
rule index_dedup_bam_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}.dedup.bam"
    output:
        bai = "results/{sample}/bwa_mtdna/{sample}.dedup.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"

rule qualimap_bamqc_bwa_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}.dedup.bam"
    output:
        txt = "results/{sample}/qualimap_mtdna/genome_results.txt"
    log:
        "logs/qualimap_mtdna/{sample}.log"
    conda:
        "workflow/envs/qualimap.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/qualimap_mtdna
        qualimap --java-mem-size=9G \
                 bamqc \
                 -bam {input.bam} \
                 -outdir results/{wildcards.sample}/qualimap_mtdna \
                 > {log} 2>&1
        """


rule damage_profiler_mtdna:
    input:
        bam = "results/{sample}/bwa_mtdna/{sample}.dedup.bam",
        species_file = "results/{sample}/fastq_screen/{sample}_best_species.txt",
        config = "config/config.yaml"  # Correct path here
    output:
        dir = directory("results/{sample}/damageprofiler_mtdna/")
    log:
        "logs/damageprofiler_mtdna/{sample}.log"
    conda:
        "workflow/envs/damageprofiler.yaml"
    shell:
        """
        mkdir -p {output.dir}

        species=$(cat {input.species_file})

        ref=$(python -c "import yaml; cfg = yaml.safe_load(open('{input.config}')); print(cfg['mtDNA_indices'].get('${{species}}', ''))")
        if [ -z "$ref" ]; then
            echo "Reference for species '$species' not found in {input.config}" >&2
            exit 1
        fi

        damageprofiler -Xmx12G \
            -i {input.bam} \
            -o {output.dir} \
            -r "$ref" > {log} 2>&1
        """



rule softclip_cram_mtdna:
    input:
        cram = "results/{sample}/bwa_mtdna/{sample}.dedup.cram",
        ref = lambda wc: config["mtDNA_indices"][
            open(f"results/{wc.sample}/fastq_screen/{wc.sample}_best_species.txt").read().strip()
        ]
    output:
        cram = "results/{sample}/bwa_mtdna/{sample}.dedup_q30_softclipped.cram"
    threads: 4
    conda: "workflow/envs/samtools.yaml"
    shell:
        """
        samtools view -h {input.cram} | \
        python2 /home/kdaly/programs/scripts_for_goat_project/softclip_mod.py - 4 | \
        samtools view -@ {threads} -T {input.ref} -O CRAM -o {output.cram}
        """

###----------------------------------------wrappers------------------------------------------------######################
rule merge_pathogen_summaries:
    input:
        summaries=expand("results/{sample}/summary/{sample}_pathogen_summary.csv", sample=SAMPLES)
    output:
        excel="results/final/pathogen_summary_all_samples.xlsx"
    conda:
        "workflow/envs/summary.yaml"
    script:
        "scripts/merge_summaries.py"
def refs_for_sample(sample):
    pairs = get_sample_ref_pairs()
    return [ref.replace(" ", "_") for s, ref in pairs if s == sample]

def summarize_inputs(wc):
    refs = refs_for_sample(wc.sample)
    # Instead of specific files inside these directories, just input the directories:
    damage_dirs = [f"results/{wc.sample}/bwa_pathogen/damageprofiler_{ref}" for ref in refs]
    qualimap_dirs = [f"results/{wc.sample}/bwa_pathogen/qualimap_{ref}" for ref in refs]

    return [
        f"results/{wc.sample}/Escore/pathogen/{wc.sample}_pathogen.csv",
        "results/hops/maltExtract/heatmap_overview_Wevid.tsv",
        "config/Pathogen_spreadsheet.csv",
    ] + damage_dirs + qualimap_dirs

def safe_name(pathogen):
    return pathogen.replace(" ", "_").replace("/", "_") 

rule summarize_pathogen_data:
    input:
        escore = "results/{sample}/Escore/pathogen/{sample}_pathogen.csv",
        hops = "results/hops/maltExtract/heatmap_overview_Wevid.tsv",
        spreadsheet = "config/Pathogen_spreadsheet.csv",
        qualimaps = lambda wildcards: [
            f"results/{wildcards.sample}/bwa_pathogen/qualimap_{safe_name(p)}"
            for s, p in get_sample_ref_pairs() if s == wildcards.sample
        ],
        damage = lambda wildcards: [
            f"results/{wildcards.sample}/bwa_pathogen/damageprofiler_{safe_name(p)}"
            for s, p in get_sample_ref_pairs() if s == wildcards.sample
        ],
    output:
        "results/{sample}/summary/{sample}_pathogen_summary.csv"
    conda:
        "workflow/envs/summary.yaml"
    shell:
        """
        python scripts/summarize_pathogen_data.py \
            --sample {wildcards.sample} \
            --escore {input.escore} \
            --hops {input.hops} \
            --spreadsheet {input.spreadsheet} \
            --bam_dir results/{wildcards.sample}/bwa_pathogen \
            --qualimap_dir results/{wildcards.sample}/bwa_pathogen \
            --damage_dir results/{wildcards.sample}/bwa_pathogen \
            --output {output}
        """


rule summarize_host_mtdna:
    input:
        host_bams = expand("results/{sample}/bwa_host/{sample}.dedup.bam", sample=SAMPLES),
        mtdna_bams = expand("results/{sample}/bwa_mtdna/{sample}.dedup.bam", sample=SAMPLES),
        qualimap_host = expand("results/{sample}/qualimap/genome_results.txt", sample=SAMPLES),
        qualimap_mtdna = expand("results/{sample}/qualimap_mtdna/genome_results.txt", sample=SAMPLES),
        samples_tsv = "config/samples.tsv",
        species = expand("results/{sample}/fastq_screen/{sample}_best_species.txt", sample=SAMPLES),
        collapsed = expand("results/{sample}/adapter_removal/{sample}.collapsed.gz", sample=SAMPLES)
    output:
        "results/final/host_mtdna_summary_all_samples.xlsx"
    params:
        samples = ",".join(SAMPLES)
    conda:
        "workflow/envs/summary.yaml"
    script:
        "scripts/summarize_host_mtdna.py"



rule pathogen_report:
    input:
        escore="results/{sample}/Escore/pathogen/{sample}_pathogen.csv",
        spreadsheet="config/Pathogen_spreadsheet.csv",
        summary="results/{sample}/summary/{sample}_pathogen_summary.csv"
    output:
        pdf="results/{sample}/summary/{sample}_{ref_name_safe}_pathogen_report_merged.pdf"
    params:
        pathogen=lambda wc: wc.ref_name_safe.replace("_", " ")
    conda:
        "workflow/envs/report_env.yaml"
    shell:
        """
        Rscript scripts/generate_pathogen_report.R {wildcards.sample} '{params.pathogen}' {input.spreadsheet} {output.pdf}
        """

# Pathogen Detection Scoring System
rule calculate_pathogen_detection_scores:
    input:
        # E-Score results (already filtered by user thresholds)
        escore_files = expand("results/{sample}/Escore/pathogen/{sample}_pathogen.csv", sample=SAMPLES),
        # Hops results
        hops_results = "results/hops/maltExtract/heatmap_overview_Wevid.tsv",
        # BWA results for detected pathogens
        bwa_results = lambda wildcards: [
            f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.ani.txt"
            for sample, pathogen in get_sample_ref_pairs()
        ],
        # Damageprofiler results
        damage_results = lambda wildcards: [
            f"results/{sample}/bwa_pathogen/damageprofiler_{safe_name(pathogen)}"
            for sample, pathogen in get_sample_ref_pairs()
        ],
        # Breadth and entropy results
        breadth_results = lambda wildcards: [
            f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.breadth_ratio.txt"
            for sample, pathogen in get_sample_ref_pairs()
        ],
        entropy_results = lambda wildcards: [
            f"results/{sample}/bwa_pathogen/{sample}_{safe_name(pathogen)}.mean_entropy.txt"
            for sample, pathogen in get_sample_ref_pairs()
        ],
        # Comparison results (for k-mer ranking)
        comparison_results = expand("results/comparison/{sample}_comparison.tsv", sample=SAMPLES)
    output:
        scores_matrix = "results/pathogen_detection/detection_scores_matrix.csv",
        scores_heatmap = "results/pathogen_detection/detection_scores_heatmap.pdf",
        detailed_scores = "results/pathogen_detection/detailed_scores.csv"
    conda:
        "workflow/envs/python.yaml"
    script:
        "scripts/calculate_detection_scores.py"
