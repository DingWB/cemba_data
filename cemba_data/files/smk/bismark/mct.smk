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

rna_bam_dir=os.path.abspath(workflow.default_remote_prefix+"/rna_bam") if config["gcp"] else "rna_bam"

if not os.path.exists(rna_bam_dir):
    os.mkdir(rna_bam_dir)
# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here,
        # once summary is generated, snakemake will delete these stats
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",cell_id=CELL_IDS,mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS,mcg_context=mcg_context),
        expand("fastq/{cell_id}-{read_type}.trimmed.stats.txt", cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}.dna_reads.bam.reads_profile.csv", cell_id=CELL_IDS),
        'rna_bam/TotalRNAAligned.filtered.bam',  # needed for count star mapped reads by RG
        'rna_bam/TotalRNALog.final.out',
        'rna_bam/TotalRNALog.out',
        'rna_bam/TotalRNALog.progress.out',
        # 'rna_bam/TotalRNAAligned.rna_reads.bam',
        'rna_bam/TotalRNAAligned.rna_reads.bam.reads_profile.csv',
        'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv',
        'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary'
    output:
        "MappingSummary.csv.gz"
    params:
        outdir="./" if not config["gcp"] else workflow.default_remote_prefix,
    shell:
        """
        yap-internal summary --output_dir {params.outdir} --fastq_dir {fastq_dir} --mode {mode} --barcode_version {barcode_version} \
--mc_stat_feature "{mc_stat_feature}" --mc_stat_alias "{mc_stat_alias}" \
--num_upstr_bases {num_upstr_bases} --mc_rate_max_threshold {mc_rate_max_threshold} \
--dna_cov_min_threshold {dna_cov_min_threshold}
        """

# Trim reads
rule trim:
    input:
        fq=get_fastq_path()
    output:
        fq=local(temp("fastq/{cell_id}-{read_type}.trimmed.fq.gz")),
        stats="fastq/{cell_id}-{read_type}.trimmed.stats.txt"
    params:
        adapter=lambda wildcards: r1_adapter if wildcards.read_type == 'R1' else r2_adapter,
        left_cut= lambda wildcards: r1_left_cut if wildcards.read_type == 'R1' else r2_left_cut,
        right_cut=lambda wildcards: r1_right_cut if wildcards.read_type == 'R1' else r2_right_cut,
    threads:
        2
    shell:
        "cutadapt -a R1Adapter={params.adapter} "
        "-a TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG "
        "-a N6=AAGCAGTGGTATCAACGCAGAGTAC "
        "-a TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT "
        "-a N6_rc=GTACTCTGCGTTGATACCACTGCTT "
        "-a 3PpolyT=TTTTTTTTTTTTTTTX "
        "-g 5PpolyT=XTTTTTTTTTTTTTTT "
        "-a 3PpolyA=AAAAAAAAAAAAAAAX "
        "-g 5PpolyA=XAAAAAAAAAAAAAAA "
        "-a polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "
        "-a polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
        "-a ISPCR_F=AAGCAGTGGTATCAACGCAGAGT "
        "-a ISPCR_R=ACTCTGCGTTGATACCACTGCTT "
        "{input.fq} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {params.left_cut} -u -{params.right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

rule bismark:
    input:
        local("fastq/{cell_id}-{read_type}.trimmed.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.bam")),
        stats="bam/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt"
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
        stats="bam/{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt"
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

# select DNA reads from final bismark mapped bam
rule select_dna:
    input:
        "bam/{cell_id}.final.bam"
    output:
        bam="bam/{cell_id}.dna_reads.bam",
        stats='bam/{cell_id}.dna_reads.bam.reads_profile.csv'
    shell:
        """
        yap-internal select-dna-reads --input_bam {input} \
--output_bam {output.bam} --mc_rate_max_threshold {mc_rate_max_threshold} \
--cov_min_threshold {dna_cov_min_threshold} {nome_flag_str} --assay_type mc
        """

# generate ALLC using dna_reads.bam
rule allc:
    input:
        bam="bam/{cell_id}.dna_reads.bam"
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

# RNA mapping, also start from trimmed fastq
cell_ids_str = ' , ID:'.join(CELL_IDS)
# star separate multiple input by ,
star_input_str = ','.join([fastq_dir+f"/{cell_id}-R1.trimmed.fq.gz" for cell_id in CELL_IDS])

rule star:
    input:
        # here we only use R1 SE for RNA,
        # R2 SE or R1R2 PE is worse than R1 actually, due to R2's low quality
        # And we map all cells together, so the genome is only load once
        # each cell will have a different @RG tag
        local(expand("fastq/{cell_id}-R1.trimmed.fq.gz", cell_id = CELL_IDS))
    output:
        'rna_bam/TotalRNAAligned.out.bam',
        'rna_bam/TotalRNALog.final.out',
        'rna_bam/TotalRNALog.out',
        'rna_bam/TotalRNALog.progress.out',
        'rna_bam/TotalRNASJ.out.tab'
    params:
        prefix="rna_bam/TotalRNA" if not config["gcp"] else workflow.default_remote_prefix+"/rna_bam/TotalRNA"
    threads:
        workflow.cores * 0.8  # workflow.cores is user provided cores for snakemake
    resources:
        mem_mb=48000
    shell:
        'STAR --runThreadN {threads} '
        '--genomeDir {star_reference} '
        '--alignEndsType Local '
        '--genomeLoad NoSharedMemory '
        '--outSAMstrandField intronMotif '
        '--outSAMtype BAM Unsorted '
        '--outSAMunmapped None '
        '--outSAMattributes NH HI AS NM MD '
        '--sjdbOverhang 100 '
        '--outFilterType BySJout '  # ENCODE standard options
        '--outFilterMultimapNmax 20 '  # ENCODE standard options
        '--alignSJoverhangMin 8 '  # ENCODE standard options
        '--alignSJDBoverhangMin 1 '  # ENCODE standard options
        '--outFilterMismatchNmax 999 '  # ENCODE standard options
        '--outFilterMismatchNoverLmax 0.04 '  # ENCODE standard options
        '--alignIntronMin 20 '  # ENCODE standard options
        '--alignIntronMax 1000000 '  # ENCODE standard options
        '--alignMatesGapMax 1000000 '  # ENCODE standard options
        '--outFileNamePrefix {params.prefix} '
        '--readFilesIn {star_input_str} '
        '--readFilesCommand gzip -cd '
        '--outSAMattrRGline ID:{cell_ids_str}'

rule filter_rna_bam:
    input:
        'rna_bam/TotalRNAAligned.out.bam'
    output:
        'rna_bam/TotalRNAAligned.filtered.bam'
    threads:
        min(workflow.cores * 0.8, 10)
    shell:
        "samtools sort -@ {threads} -m 2G {input} | samtools view -bh -q 10 -o {output} -"

rule select_rna:
    input:
        'rna_bam/TotalRNAAligned.filtered.bam'
    output:
        bam='rna_bam/TotalRNAAligned.rna_reads.bam',
        stats='rna_bam/TotalRNAAligned.rna_reads.bam.reads_profile.csv'
    shell:
        """
        yap-internal select-rna-reads  --input_bam {input}  --output_bam {output.bam}  \
--mc_rate_min_threshold {mc_rate_min_threshold} --cov_min_threshold {rna_cov_min_threshold} {nome_flag_str} --assay_type mc
        """

rule feature_count:
    input:
        'rna_bam/TotalRNAAligned.rna_reads.bam'
    output:
        count_tsv='rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv',
        stats='rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary'
    threads:
        min(workflow.cores * 0.8, 10)
    resources:
        mem_mb=1000
    shell:
        """
        featureCounts -t {feature_type} -g {id_type} -a {gtf_path} -o {output.count_tsv} --byReadGroup -T {threads} {input}
        """
