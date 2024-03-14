import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'mct_base.smk')

# the summary rule is the final target
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt",cell_id=CELL_IDS),

        # bam dir
        expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv",
                        cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.reads_mch_frac.csv", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.bam.bai",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam.bai",
                    cell_id=CELL_IDS),

        # rna mapping
        expand("bam/{cell_id}.hisat3n_rna_summary.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam.bai",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv.summary",cell_id=CELL_IDS),

        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz",cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi",cell_id=CELL_IDS),

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
        aggregate_feature_counts(indir=indir)
        snmct_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

module hisat3n:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
    config: config

# use rule * from hisat3n exclude trim as hisat3n_*
use rule sort_fq,unique_reads_cgn_extraction,index_bam from hisat3n

module mct:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            os.path.join(PACKAGE_DIR,"files","smk",'hisat3n',"mct.smk")

use rule * from mct exclude sort_dna_bam as mct_*

rule sort_dna_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unsort.bam"),
    output:
        bam = local(temp(bam_dir + "/{cell_id}.hisat3n_dna.bam")),
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

# Separate unique aligned reads and multi-aligned reads with length > 30
# TODO: make sure how to separate multi-align reads? or reads map to repeat regions in the genome?
# TODO right now, we are just using mapq == 1 as multi-align reads, but this might not be right
rule split_unique_and_multi_align_bam_dna:
    input:
        bam = local(bam_dir + "/{cell_id}.hisat3n_dna.bam"),
    output:
        unique=local(temp(bam_dir + "/{cell_id}.hisat3n_dna.unique_align.bam")),
        multi=local(temp(bam_dir + "/{cell_id}.hisat3n_dna.multi_align.bam"))
    run:
        separate_unique_and_multi_align_reads(
            in_bam_path=input.bam,
            out_unique_path=output.unique,
            out_multi_path=output.multi,
            out_unmappable_path=None,
            mapq_cutoff=10,
            qlen_cutoff=30
        )

# remove PCR duplicates
rule dedup_multi_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.multi_align.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.multi_align.deduped.bam")),
        stats="bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt"
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR bam/temp/
        """

rule select_multi_bam_dna_reads:
    input:
        bam=local(bam_dir + "/{cell_id}.hisat3n_dna.multi_align.deduped.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.bam",
        stats="bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.reads_mch_frac.csv"

    resources:
        mem_mb=100
    run:
        select_mct_reads(
            input_bam=input.bam,
            output_bam=output.bam,
            mode='dna',
            mc_rate_max_threshold=0.5,
            cov_min_threshold=3,
            nome=False
        )


# ==================================================
# Generate ALLC
# ==================================================

rule multi_reads_allc:
    input:
        bam = "bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.bam",
        bai = "bam/{cell_id}.hisat3n_dna.multi_align.deduped.dna_reads.bam.bai"
    output:
        allc="allc-multi/{cell_id}.allc_multi.tsv.gz",
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
--cpu {threads} --num_upstr_bases {config[num_upstr_bases]} \
--num_downstr_bases {config[num_downstr_bases]} \
--compress_level {config[compress_level]} --save_count_df \
--min_mapq 0 --convert_bam_strandness
        """
