"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import yaml
import pathlib
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'bismark_base.smk')

# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # once summary is generated, snakemake will delete these stats
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",cell_id=CELL_IDS,mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS,mcg_context=mcg_context),
        expand("fastq/{cell_id}-{read_type}.trimmed.stats.tsv", cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS,read_type=['R1','R2']),
    output:
        "MappingSummary.csv.gz"
    params:
        outdir="./" if not config["gcp"] else workflow.default_remote_prefix,
    shell:
        """
        yap-internal summary --output_dir {params.outdir} --fastq_dir {fastq_dir} \
                    --mode {mode} --barcode_version {barcode_version} \
                    --mc_stat_feature "{mc_stat_feature}" --mc_stat_alias "{mc_stat_alias}" \
                    --num_upstr_bases {num_upstr_bases}
        """

# Trim reads
rule trim:
    input:
        fq=get_fastq_path(),
    output:
        fq=local(temp("fastq/{cell_id}-{read_type}.trimmed.fq.gz")),
        stats="fastq/{cell_id}-{read_type}.trimmed.stats.tsv"
    params:
        adapter=lambda wildcards: r1_adapter if wildcards.read_type=='R1' else r2_adapter, #r1_adapter, r2_adapter and other config_str will be written into the header
        left_cut=lambda wildcards: r1_left_cut if wildcards.read_type=='R1' else r2_left_cut,
        right_cut= lambda wildcards: r1_right_cut if wildcards.read_type == 'R1' else r2_right_cut,
    threads:
        2
    shell:
        """
        cutadapt --report=minimal -a {params.adapter} {input.fq} 2> {output.stats} | cutadapt --report=minimal -O 6 -q 20 -u {params.left_cut} -u -{params.right_cut} -m 30 -o {output.fq} - >> {output.stats}
        """

# bismark mapping
rule bismark:
    input:
        local("fastq/{cell_id}-{read_type}.trimmed.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.bam")),
        stats=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt"))
    params:
        mode=lambda wildcards: "--pbat" if wildcards.read_type == "R1" else ""
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode
        """
        bismark {bismark_reference} {unmapped_param_str} --bowtie2 {input} {params.mode} -o {bam_dir} --temp_dir {bam_dir}
        """

# filter bam
rule filter_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.bam")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.filter.bam"))
    shell:
        "samtools view -b -h -q 10 -o {output} {input}"

# sort bam by position
rule sort_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.filter.bam")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.sorted.bam"))
    resources:
        mem_mb=1000
    shell:
        """
        samtools sort -o {output} {input}
        """

# remove PCR duplicates
rule dedup_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.sorted.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.bam")),
        stats=bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt"
    params:
        tmp_dir="bam/temp" if not config["gcp"] else workflow.default_remote_prefix+"/bam/temp"
    resources:
        mem_mb=3000
    shell:
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR {params.tmp_dir}
        """

# merge R1 and R2, get final bam
rule merge_bam:
    input:
        local(bam_dir+"/{cell_id}-R1.trimmed_bismark_bt2.deduped.bam"),
        local(bam_dir+"/{cell_id}-R2.trimmed_bismark_bt2.deduped.bam")
    output:
        "bam/{cell_id}.final.bam"
    shell:
        "samtools merge -f {output} {input}"

# generate ALLC
rule allc:
    input:
        bam="bam/{cell_id}.final.bam"
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        tbi= "allc/{cell_id}.allc.tsv.gz.tbi",
        stats="allc/{cell_id}.allc.tsv.gz.count.csv"
    threads:
        2
    resources:
        mem_mb=500
    shell:
        """
        mkdir -p {allc_dir}
        allcools bam-to-allc \
                --bam_path {input.bam} \
                --reference_fasta {reference_fasta} \
                --output_path {output.allc} \
                --cpu 1 \
                --num_upstr_bases {num_upstr_bases} \
                --num_downstr_bases {num_downstr_bases} \
                --compress_level {compress_level} \
                --save_count_df
        """

# CGN extraction from ALLC
rule cgn_extraction:
    input:
        allc="allc/{cell_id}.allc.tsv.gz",
        tbi="allc/{cell_id}.allc.tsv.gz.tbi"
    output:
        allc="allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",
        tbi="allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",
    params:
        prefix=allc_mcg_dir+"/{cell_id}",
    threads:
        1
    resources:
        mem_mb=100
    shell:
        """
        mkdir -p {allc_mcg_dir}
        allcools extract-allc --strandness merge \
--allc_path  {input.allc} --output_prefix {params.prefix} \
--mc_contexts {mcg_context} --chrom_size_path {chrom_size_path}
        """
