import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'mct_base.smk')

# ==================================================
# Mapping summary
# ==================================================
# the summary rule is the final target
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt",cell_id=CELL_IDS),

        # bam dir
        expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv",
                        cell_id=CELL_IDS), # local(expand("bam/{cell_id}.*", cell_id=CELL_IDS)),
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
        indir='.' if not config["gcp"] else workflow.default_remote_prefix
        aggregate_feature_counts(indir=indir)
        snmct_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

module hisat3n:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
    config: config
    # skip_validation: True

# use rule * from hisat3n exclude trim as hisat3n_*
use rule sort_fq,unique_reads_cgn_extraction from hisat3n

rule trim:
    input:
        # change to sort_R1 and sort_R2 output if the FASTQ name is disordered
        R1=local("fastq/{cell_id}-R1_sort.fq"),
        R2=local("fastq/{cell_id}-R2_sort.fq")
    output:
        R1=local(temp("fastq/{cell_id}-R1.trimmed.fq.gz")),
        R2=local(temp("fastq/{cell_id}-R2.trimmed.fq.gz")),
        stats="fastq/{cell_id}.trimmed.stats.txt"
    threads:
        1
    shell:
        """
        cutadapt -a R1Adapter={config[r1_adapter]} -a TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG \
-a TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT -a N6=AAGCAGTGGTATCAACGCAGAGTAC \
-a N6_rc=GTACTCTGCGTTGATACCACTGCTT -a 3PpolyT=TTTTTTTTTTTTTTTX \
-a 3PpolyA=AAAAAAAAAAAAAAAX -a polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
-a polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -a ISPCR_F=AAGCAGTGGTATCAACGCAGAGT \
-a ISPCR_R=ACTCTGCGTTGATACCACTGCTT \
-A R2Adapter={config[r2_adapter]} -A TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG \
-A TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT -A N6=AAGCAGTGGTATCAACGCAGAGTAC \
-A N6_rc=GTACTCTGCGTTGATACCACTGCTT -A 3PpolyT=TTTTTTTTTTTTTTTX \
-A 3PpolyA=AAAAAAAAAAAAAAAX -A polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
-A polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -A ISPCR_F=AAGCAGTGGTATCAACGCAGAGT \
-A ISPCR_R=ACTCTGCGTTGATACCACTGCTT \
-g 5PpolyT=XTTTTTTTTTTTTTTT \
-g 5PpolyA=XAAAAAAAAAAAAAAA -G 5PpolyT=XTTTTTTTTTTTTTTT \
-G 5PpolyA=XAAAAAAAAAAAAAAA --report=minimal \
-O 6 -q 20 -u {config[r1_left_cut]} -u -{config[r1_right_cut]} -U {config[r2_left_cut]} \
-U -{config[r2_right_cut]} -Z -m {config[min_read_length]}:{config[min_read_length]} \
--pair-filter 'both' -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {output.stats}
        """

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

rule sort_dna_bam:
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

# remove PCR duplicates
rule dedup_unique_bam:
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

rule select_unique_bam_dna_reads:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam"
    output:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam",
        stats="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv"
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
# HISAT-3N RNA Mapping
# ==================================================

# Paired-end Hisat3n mapping using RNA mode
rule hisat2_pairend_mapping_rna_mode:
    input:
        R1=local("fastq/{cell_id}-R1.trimmed.fq.gz"),
        R2=local("fastq/{cell_id}-R2.trimmed.fq.gz")
    output:
        bam= local(temp(bam_dir + "/{cell_id}.hisat3n_rna.unsort.bam")),
        stats="bam/{cell_id}.hisat3n_rna_summary.txt",
    threads:
        config["hisat_threads"]
    resources:
        mem_mb=8000
    shell: # add read group @RG to the reads in order to use featuerCounts
        """
        hisat2 -x {config[hisat_rna_reference]} -q -1 {input.R1} -2 {input.R2} -t --new-summary \
--summary-file {output.stats} --threads {threads} \
| samtools addreplacerg -r '@RG\tID:{wildcards.cell_id}' -u -o - - \
| samtools view -b -q 10 -o {output.bam}
        """

rule sort_rna_bam:
    input:
        bam=local(bam_dir + "/{cell_id}.hisat3n_rna.unsort.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_rna.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

# skip dedup step for RNA reads
rule select_unique_bam_rna_reads:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_rna.bam")
    output:
        bam="bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam",
        stats="bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv"
    resources:
        mem_mb=100
    run:
        select_mct_reads(
            input_bam=input.bam,
            output_bam=output.bam,
            mode='rna',
            mc_rate_min_threshold=0.9,
            cov_min_threshold=3,
            nome=False
        )

rule feature_count:
    input:
        bam="bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam",
    output:
        tsv='bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv',
        stats='bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv.summary'
    threads:
        1
    resources:
        mem_mb=1000
    shell: #version 2.0.1, if there is overlap between two records in gtf, then there will be no reads assigned to these two features.
        """
        featureCounts -t {config[feature_type]} -g {config[id_type]} \
-a {config[gtf_path]} -o {output.tsv} --byReadGroup -T {threads} {input}
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
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam.bai"
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
--cpu {threads} --num_upstr_bases {config[num_upstr_bases]} \
--num_downstr_bases {config[num_downstr_bases]} \
--compress_level {config[compress_level]} --save_count_df \
--convert_bam_strandness
        """

