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

        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz",cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi",cell_id=CELL_IDS),

        # mhap
        expand("mhap/{cell_id}.mhap.gz",cell_id=CELL_IDS),
        expand("mhap/{cell_id}.mhap.gz.tbi",cell_id=CELL_IDS),

        # allc-CGN
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi", cell_id=CELL_IDS, mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",cell_id=CELL_IDS,mcg_context=mcg_context)
    output:
        csv="MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        # generate the final summary
        indir='.' if not config["gcp"] else workflow.default_remote_prefix
        snmc_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

module hisat3n:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
    config: config

# use rule * from hisat3n exclude dedup,unique_reads_allc,hisat_3n_pair_end_mapping_dna_mode,index_bam as hisat3n_*
use rule sort_fq,trim,unique_reads_cgn_extraction from hisat3n as hisat3n_*

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
    shell: # -q 10 will filter out multi-aligned reads
        """
        mkdir -p {bam_dir}
        hisat-3n {config[hisat3n_dna_reference]} -q  -1 {input.R1} -2 {input.R2} \
--directional-mapping-reverse --base-change C,T {repeat_index_flag} \
--no-spliced-alignment --no-temp-splicesite -t  --new-summary \
--summary-file {output.stats} --threads {threads} | samtools view -b -q 10 -o {output.bam}
        """

rule mc_sort_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unsort.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.unique_align.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

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

rule index_bam:
    input:
        bam="{input_name}.bam"
    output:
        bai="{input_name}.bam.bai"
    shell:
        """
        samtools index {input.bam}
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

# Convert bam to mhap
rule bam_to_mhap:
    input: #sorted bam
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
    output:
        mhap="mhap/{cell_id}.mhap.gz",
        tbi="mhap/{cell_id}.mhap.gz.tbi"
    params:
        cpgPath=os.path.expanduser(config['cpg_path']),
    resources:
        mem_mb=400
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        outfile=output.mhap[:-3] #"allc/{cell_id}.mhap", will be bgzipped and tabix indexed in mhap
        bam2mhap(bam_path=input.bam,cpg_path=params.cpgPath,output=outfile)

# rule stat_mhap_gene:
#     input: #sorted .mhap.gz
#         mhap="mhap/{cell_id}.mhap.gz",
#         tbi="mhap/{cell_id}.mhap.gz.tbi"
#     output:
#         gene_stat="mhap/{cell_id}.mhap.gene.stat.tsv.gz",
#     params:
#         cpgPath=os.path.expanduser(config['cpg_path']),
#         geneBedPath=os.path.expanduser(config['gene_bed_path']),
#     resources:
#         mem_mb=400
#     run:
#         from cemba_data.mapping.pipelines import stat_mhap
#         if not os.path.exists(mhap_dir):
#             os.mkdir(mhap_dir)
#         stat_mhap(mhap_path=input.mhap,cpg_path=params.cpgPath,
#                 region=None,bed=params.geneBedPath,output=output.gene_stat)