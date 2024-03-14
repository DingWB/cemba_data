import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'base.smk')

include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')

# ==================================================
# Mapping summary
# ==================================================
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt",cell_id=CELL_IDS),

        # bam dir
        expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna_split_reads_summary.{read_type}.txt",
                        cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam",cell_id=CELL_IDS),

        # 3C contacts
        expand("hic/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv", cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.3C.contact.tsv.gz", cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.dedup_contacts.tsv.gz", cell_id=CELL_IDS),

        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi", cell_id=CELL_IDS),

        #allc-multi
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv",cell_id=CELL_IDS),
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz", cell_id=CELL_IDS),
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz.tbi", cell_id=CELL_IDS),

        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",cell_id=CELL_IDS, mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS, mcg_context=mcg_context),
    output:
        csv="MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        # generate the final summary
        indir = '.' if not gcp else workflow.default_remote_prefix
        snm3c_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")


#=====================================================================
# Processing multi-alignment reads and generate allc-multi
#=====================================================================
rule sort_multi_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.multi_aligned.bam"), #"bam/{cell_id}.hisat3n_dna.multi_aligned.bam"
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna_sorted.multi_align.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

rule dedup_multi_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna_sorted.multi_align.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        stats="bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt"
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        """
        picard MarkDuplicates -I {input} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR bam/temp/
        """

# generate ALLC
rule multi_reads_allc:
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
--cpu {threads} --num_upstr_bases {config[num_upstr_bases]} \
--num_downstr_bases {config[num_downstr_bases]} \
--compress_level {config[compress_level]} --save_count_df \
--min_mapq 0 --convert_bam_strandness
        """

