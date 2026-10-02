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
        # methylation output (allc and/or cz, controlled by mc_format)
        get_methylation_targets(CELL_IDS),
        # mhap (optional, controlled by generate_mhap)
        get_mhap_targets(CELL_IDS),
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
                    --num_upstr_bases {num_upstr_bases} --mc_format {mc_format}
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
# no --convert_bam_strandness: bismark SE flags already match the converted strand
# (identical calls to bam_to_cz(convert_bam_strandness=True), without a temp BAM)
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
                --chroms {chrom_size_path} \
                --save_count_df
        """

# index the final bam (needed by cz / mhap conversion)
rule index_final_bam:
    input:
        bam="bam/{cell_id}.final.bam"
    output:
        bai="bam/{cell_id}.final.bam.bai"
    shell:
        "samtools index {input.bam}"

# ==================================================
# Generate CZ (cytozip), alternative/addition to ALLC
# ==================================================
rule bam_to_cz:
    input:
        bam="bam/{cell_id}.final.bam",
        bai="bam/{cell_id}.final.bam.bai"
    output:
        cz="cz/{cell_id}.cz"
    threads:
        2
    resources:
        mem_mb=1300
    run:
        from cytozip import bam_to_cz
        os.makedirs(cz_dir, exist_ok=True)
        if reference_cz in (None, '', 'None') or \
                not os.path.exists(os.path.expanduser(str(reference_cz))):
            raise FileNotFoundError(
                "reference_cz is required to generate .cz files. Build one with:\n"
                f"    czip build_ref -g {reference_fasta} "
                f"-O <output.allc.cz> -s {chrom_size_path} -j 20\n"
                "then set 'reference_cz = <output.allc.cz>' in the mapping config.")
        bam_to_cz(
            bam_path=input.bam,
            genome=os.path.expanduser(str(reference_fasta)),
            output=output.cz,
            reference=os.path.expanduser(str(reference_cz)),
            num_upstr_bases=int(num_upstr_bases),
            num_downstr_bases=int(num_downstr_bases),
            convert_bam_strandness=True,
            chroms=os.path.expanduser(str(chrom_size_path)),
            save_count_df=True)

# ==================================================
# Convert bam to mhap (optional, enabled by generate_mhap)
# ==================================================
rule bam_to_mhap:
    input:
        bam="bam/{cell_id}.final.bam",
        bai="bam/{cell_id}.final.bam.bai"
    output:
        mhap_cg="mhap/{cell_id}.CG.mhap.gz",
        tbi_cg="mhap/{cell_id}.CG.mhap.gz.tbi",
        mhap_ch="mhap/{cell_id}.CH.mhap.gz",
        tbi_ch="mhap/{cell_id}.CH.mhap.gz.tbi"
    resources:
        mem_mb=500
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        os.makedirs(mhap_dir, exist_ok=True)
        if annotation_path in (None, '', 'None'):
            raise ValueError(
                "generate_mhap=True requires 'annotation_path' (path to the "
                "*_allc.gz annotation) in the mapping config.")
        annotation = os.path.expanduser(str(annotation_path))
        outfile_cg = output.mhap_cg[:-3]  # strip ".gz"; bgzipped + tabixed by bam2mhap
        bam2mhap(bam_path=input.bam, annotation=annotation,
                 output=outfile_cg, pattern="CGN")
        outfile_ch = output.mhap_ch[:-3]
        bam2mhap(bam_path=input.bam, annotation=annotation,
                 output=outfile_ch, pattern="CHN")

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
