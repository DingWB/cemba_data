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

mhap_dir=os.path.abspath(workflow.default_remote_prefix+"/mhap") if config["gcp"] else "mhap"
if not os.path.exists(mhap_dir):
    os.mkdir(mhap_dir)
# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here, so they won't be deleted until summary is generated
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        # mhap
        expand("mhap/{cell_id}.mhap.gz",cell_id=CELL_IDS),
        expand("mhap/{cell_id}.mhap.gz.tbi",cell_id=CELL_IDS),
        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",cell_id=CELL_IDS,mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS,mcg_context=mcg_context),

        expand("fastq/{cell_id}-{read_type}.trimmed.stats.tsv",cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}-{read_type}.merged.deduped.matrix.txt",cell_id=CELL_IDS,read_type=['R1','R2']),
        local(expand(bam_dir+"/{cell_id}-{read_type}.merged.filter.bam",cell_id=CELL_IDS,read_type=['R1','R2'])),
        local(expand(bam_dir+"/{cell_id}-{read_type}.merged.deduped.bam",cell_id=CELL_IDS,read_type=['R1','R2'])),
        expand("hic/{cell_id}.3C.contact.tsv.gz", cell_id=CELL_IDS),
        expand("hic/{cell_id}.3C.contact.tsv.counts.txt", cell_id=CELL_IDS)
    output:
        "MappingSummary.csv.gz"
    params:
        outdir="./" if not config["gcp"] else workflow.default_remote_prefix,
    shell:
        """
        yap-internal summary --output_dir {params.outdir} --fastq_dir {fastq_dir} --mode {mode} --barcode_version {barcode_version} \
--mc_stat_feature "{mc_stat_feature}" --mc_stat_alias "{mc_stat_alias}" \
--num_upstr_bases {num_upstr_bases}
        """

# Trim reads
rule trim:
    input:
        fq=get_fastq_path()
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
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark.bam")),
        um=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.fq.gz")),
        stats=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_SE_report.txt"))
    params:
        mode=lambda wildcards: "--pbat" if wildcards.read_type=="R1" else ""
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode; map R2 with normal SE mode
        """
        mkdir -p {bam_dir}
        bismark {bismark_reference} -un --bowtie1 {input} {params.mode} -o {bam_dir} --temp_dir {bam_dir}
        """

# split unmapped fastq
rule split_um_fastq:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.fq.gz")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split.fq.gz"))
    threads:
        1
    shell:
        """
        yap-internal m3c-split-reads --fastq_path {input} --output_path {output} --size_l {split_left_size} --size_r {split_right_size} --size_m {split_middle_min_size} --trim_b {trim_on_both_end}
        """

# map split fastq again
rule bismark_split:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark.bam")),
        stats=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark_SE_report.txt"))
    params:
        mode=lambda wildcards: "--pbat" if wildcards.read_type=="R1" else ""
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        """
        bismark {bismark_reference} --bowtie1 {input} {params.mode} -o {bam_dir} --temp_dir {bam_dir}
        """

# merge two bam files
rule merge_raw_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark.bam"),
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark.bam")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.merged.bam"))
    shell:
        """
        samtools merge -f {output} {input}
        """

# filter bam
rule filter_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.merged.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.merged.filter.bam"))
    shell:
        """
        samtools view -b -h -q 10 -o {output.bam} {input}
        """

# sort bam by coords
rule sort_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.merged.filter.bam")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.merged.sorted.bam"))
    resources:
        mem_mb=1000
    shell:
        """
        samtools sort -o {output} {input}
        """

# remove PCR duplicates
rule dedup_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.merged.sorted.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.merged.deduped.bam")),
        stats=bam_dir+"/{cell_id}-{read_type}.merged.deduped.matrix.txt"
    params:
        tmp_dir="bam/temp" if not config["gcp"] else workflow.default_remote_prefix+"/bam/temp"
    resources:
        mem_mb=3000
    shell:
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR {params.tmp_dir}
        """

# merge R1 and R2, get final bam for mC calling
rule merge_mc_bam:
    input:
        local(bam_dir+"/{cell_id}-R1.merged.deduped.bam"),
        local(bam_dir+"/{cell_id}-R2.merged.deduped.bam")
    output:
        bam="bam/{cell_id}.mC.bam",
        bai="bam/{cell_id}.mC.bam.bai"
    shell:
        """
        samtools merge -f {output.bam} {input} && samtools index {output.bam}
        """

# generate ALLC
rule allc:
    input:
        bam="bam/{cell_id}.mC.bam",
        index="bam/{cell_id}.mC.bam.bai"
    output:
        allc = "allc/{cell_id}.allc.tsv.gz",
        tbi = "allc/{cell_id}.allc.tsv.gz.tbi",
        stats = "allc/{cell_id}.allc.tsv.gz.count.csv"
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

# Convert bam to mhap
rule bam_to_mhap:
    input: #sorted bam
        bam="bam/{cell_id}.mC.bam",
        bai="bam/{cell_id}.mC.bam.bai"
    output:
        mhap="mhap/{cell_id}.mhap.gz",
        tbi="mhap/{cell_id}.mhap.gz.tbi"
    params:
        cpgPath=os.path.expanduser(config['cpgPath']),
    resources:
        mem_mb=500
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        outfile=output.mhap[:-3]
        bam2mhap(bam_path=input.bam,cpg_path=params.cpgPath,output=outfile)

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

# merge and sort (by read name) bam before dedup for generating contact
# contact dedup happen within generate contact
rule merge_3c_bam_for_contact:
    input:
        local(bam_dir+"/{cell_id}-R1.merged.sorted.bam"),
        local(bam_dir+"/{cell_id}-R2.merged.sorted.bam")
    output:
        local(temp(bam_dir+"/{cell_id}.3C.bam"))
    shell:
        """
        samtools merge -f {output} {input}
        """

# sort by name
rule sort_bam_for_contact:
    input:
        local(bam_dir+"/{cell_id}.3C.bam")
    output:
        "bam/{cell_id}.3C.sorted.bam"
    resources:
        mem_mb=1000
    shell:
        """
        samtools sort -n -o {output} {input}
        """

rule generate_contact:
    input:
        "bam/{cell_id}.3C.sorted.bam"
    output:
        contact="hic/{cell_id}.3C.contact.tsv.gz",
        stats="hic/{cell_id}.3C.contact.tsv.counts.txt"
    resources:
        mem_mb=300
    shell:
        """
        mkdir -p {hic_dir}
        yap-internal generate-contacts --bam_path {input} --output_path {output.contact} \
--chrom_size_path {chrom_size_path} --min_gap {min_gap}
        """
