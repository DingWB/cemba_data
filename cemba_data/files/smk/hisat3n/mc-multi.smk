import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n_base.smk')

# the summary rule is the final target
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt",cell_id=CELL_IDS),

        # bam dir
        expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt", cell_id=CELL_IDS),

        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz",cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi",cell_id=CELL_IDS),

        #allc-multi
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv",cell_id=CELL_IDS),
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz",cell_id=CELL_IDS),
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz.tbi",cell_id=CELL_IDS),

        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi", cell_id=CELL_IDS, mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS,mcg_context=mcg_context)
    output:
        csv="MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        # generate the final summary
        indir='.' if not gcp else workflow.default_remote_prefix
        snmc_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

module hisat3n:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
    config: config

use rule * from hisat3n exclude unique_reads_allc,hisat_3n_pair_end_mapping_dna_mode,index_bam as hisat3n_*

rule hisat_3n_pair_end_mapping_dna_mode:
    input:
        R1=local("fastq/{cell_id}-R1.trimmed.fq.gz"),
        R2=local("fastq/{cell_id}-R2.trimmed.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.unsort.bam")),
        stats="bam/{cell_id}.hisat3n_dna_summary.txt",
    threads:
        config['hisat3n_threads']
    resources:
        mem_mb=14000
    shell: # -q 1
        """
        hisat-3n {config[hisat3n_dna_reference]} -q  -1 {input.R1} -2 {input.R2} \
--directional-mapping-reverse --base-change C,T {repeat_index_flag} \
--no-spliced-alignment --no-temp-splicesite -t  --new-summary \
--summary-file {output.stats} --threads {threads} | samtools view -b -q 1 -o {output.bam}
        """

rule mc_multi_sort_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unsort.bam") #output of rule hisat_3n_pair_end_mapping_dna_mode from hisat3n.smk
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

# Separate unique aligned reads and multi-aligned reads with length > 30
rule split_unique_and_multi_align_bam_dna:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.bam")
    output:
        unique=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.unique_align.bam")),
        multi=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.multi_align.bam"))
    run:
        separate_unique_and_multi_align_reads(
            in_bam_path=input.bam,
            out_unique_path=output.unique,
            out_multi_path=output.multi,
            out_unmappable_path=None,
            mapq_cutoff=10,
            qlen_cutoff=30
        )

rule mc_dedup_unique_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unique_align.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        stats="bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt"
    resources:
        mem_mb=2000
    threads:
        2
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR bam/temp/
        """

# ==================================================
# Generate ALLC
# ==================================================
rule unique_reads_allc:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        tbi="allc/{cell_id}.allc.tsv.gz.tbi",
        stats="allc/{cell_id}.allc.tsv.gz.count.csv"
    threads:
        1.5
    resources:
        mem_mb=500
    shell:
        """
        mkdir -p {allc_dir}
        allcools bam-to-allc --bam_path {input.bam} \
--reference_fasta {config[reference_fasta]} --output_path {output.allc} \
--num_upstr_bases {config[num_upstr_bases]} \
--num_downstr_bases {config[num_downstr_bases]} \
--compress_level {config[compress_level]} --save_count_df \
--convert_bam_strandness
        """

rule dedup_multi_bam: #dedup_unique_bam is included in mc.smk
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.multi_align.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        stats="bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt"
    resources:
        mem_mb=3000
    threads:
        2
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR bam/temp/
        """

# ==================================================
# Generate ALLC
# ==================================================
rule index_bam:
    input:
        bam="{input_name}.bam"
    output:
        bai="{input_name}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule multi_reads_allc: #unique reads allc is included in rule: mc_unique_reads_allc from mc.smk
    input:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam.bai"
    output:
        allc="allc-multi/{cell_id}.allc_multi.tsv.gz",
        tbi="allc-multi/{cell_id}.allc_multi.tsv.gz.tbi",
        stats="allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv"
    threads:
        1.5
    resources:
        mem_mb=500
    shell:
        """
        mkdir -p {allc_multi_dir}
        allcools bam-to-allc --bam_path {input.bam} \
--reference_fasta {config[reference_fasta]} --output_path {output.allc} \
--num_upstr_bases {config[num_upstr_bases]} \
--num_downstr_bases {config[num_downstr_bases]} \
--compress_level {config[compress_level]} --save_count_df \
--min_mapq 0 --convert_bam_strandness
        """
