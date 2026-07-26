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

        # methylation output (allc and/or cz, controlled by config['mc_format'])
        get_methylation_targets(CELL_IDS),

        # mhap (optional, controlled by config['generate_mhap'])
        get_mhap_targets(CELL_IDS)
    output:
        csv="MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        # generate the final summary
        indir='.' if not config["gcp"] else workflow.default_remote_prefix
        snmc_summary(outname=output.csv,indir=indir,mc_format=config['mc_format'])

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

module hisat3n:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
    config: config

# use rule * from hisat3n exclude unique_reads_allc,hisat_3n_pair_end_mapping_dna_mode,index_bam as hisat3n_*
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
        mem_mb=4000
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

# ==================================================
# Generate CZ (cytozip), alternative/addition to ALLC
# ==================================================
rule unique_reads_cz:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
    output:
        cz="cz/{cell_id}.cz"
    threads:
        1.5
    resources:
        mem_mb=1300
    run:
        from cytozip import bam_to_cz
        os.makedirs(cz_dir, exist_ok=True)
        reference_cz = config.get('reference_cz', None)
        if reference_cz in (None, '', 'None') or \
                not os.path.exists(os.path.expanduser(str(reference_cz))):
            raise FileNotFoundError(
                "reference_cz is required to generate .cz files. Build one with:\n"
                f"    czip build_ref -g {config['reference_fasta']} "
                f"-O <output.allc.cz> -s {config['chrom_size_path']} -j 20\n"
                "then set 'reference_cz = <output.allc.cz>' in the mapping config.")
        bam_to_cz(
            bam_path=input.bam,
            genome=os.path.expanduser(str(config['reference_fasta'])),
            output=output.cz,
            reference=os.path.expanduser(str(reference_cz)),
            num_upstr_bases=int(config['num_upstr_bases']),
            num_downstr_bases=int(config['num_downstr_bases']),
            convert_bam_strandness=True,
            save_count_df=True)

# ==================================================
# Generate mhap (optional, enabled by generate_mhap)
# ==================================================
rule unique_reads_mhap:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
    output:
        mhap_cg="mhap/{cell_id}.CG.mhap.gz",
        tbi_cg="mhap/{cell_id}.CG.mhap.gz.tbi",
        mhap_ch="mhap/{cell_id}.CH.mhap.gz",
        tbi_ch="mhap/{cell_id}.CH.mhap.gz.tbi"
    resources:
        mem_mb=400
    threads:
        1
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        os.makedirs(mhap_dir, exist_ok=True)
        annotation = config.get('annotation_path', None)
        if annotation in (None, '', 'None'):
            raise ValueError(
                "generate_mhap=True requires 'annotation_path' (path to the "
                "*_allc.gz annotation) in the mapping config.")
        annotation = os.path.expanduser(str(annotation))
        outfile_cg = output.mhap_cg[:-3]  # strip ".gz"; bgzipped + tabixed by bam2mhap
        bam2mhap(bam_path=input.bam, annotation=annotation,
                 output=outfile_cg, pattern="CGN")
        outfile_ch = output.mhap_ch[:-3]
        bam2mhap(bam_path=input.bam, annotation=annotation,
                 output=outfile_ch, pattern="CHN")

