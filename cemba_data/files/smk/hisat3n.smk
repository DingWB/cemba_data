import os
from cemba_data.hisat3n import *
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n_base.smk')
# ==================================================
# FASTQ Trimming
# ==================================================
# print(config)

# Trim reads
# sort the fastq files so that R1 and R2 are in the same order
rule sort_fq:
    input:
        fq=get_fastq_path(),
    output:
        fq=local(temp("fastq/{cell_id}-{read_type}_sort.fq")),
    threads:
        1.5
    resources:
        high_io_job=1,
        mem_mb=200
    shell:
        'zcat {input.fq} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output.fq} '

rule trim:
    input:
        # change to R1_sort and R2_sort output if the FASTQ name is disordered
        R1=local("fastq/{cell_id}-R1_sort.fq"),  #if local_fastq else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-R1.fq.gz"),
        R2=local("fastq/{cell_id}-R2_sort.fq")  #if local_fastq else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-R2.fq.gz")
    output:
        R1=local(temp("fastq/{cell_id}-R1.trimmed.fq.gz")),
        R2=local(temp("fastq/{cell_id}-R2.trimmed.fq.gz")),
        stats="fastq/{cell_id}.trimmed.stats.txt"
    resources:
        mem_mb=200
    threads:
        1
    shell:
        """
        cutadapt -a R1Adapter={config[r1_adapter]} \
-A R2Adapter={config[r2_adapter]} --report=minimal \
-O 6 -q 20 -u {config[r1_left_cut]} -u -{config[r1_right_cut]} \
-U {config[r2_left_cut]} -U -{config[r2_right_cut]} -Z \
-m {config[min_read_length]}:{config[min_read_length]} \
--pair-filter 'both' -o {output.R1} -p {output.R2} \
{input.R1} {input.R2} > {output.stats}
        """


# ==================================================
# HISAT-3N DNA Mapping
# ==================================================

# Paired-end Hisat3n mapping using DNA mode
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
        mem_mb=14000 # lambda wc, input: max(25 * input.size_mb, 14000); Request 4 GB of memory for this rule
    # benchmark: # https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake
    #         "fastq/{cell_id}.hisat_3n_pair_end_mapping_dna_mode.benchmark.txt"
    shell: # # do not filter any reads in this step
        """
        mkdir -p {bam_dir}
        hisat-3n {config[hisat3n_dna_reference]} -q  -1 {input.R1} -2 {input.R2} \
--directional-mapping-reverse --base-change C,T {repeat_index_flag} \
--no-spliced-alignment --no-temp-splicesite -t  --new-summary \
--summary-file {output.stats} --threads {threads} | samtools view -b -q 0 -o {output.bam}
        """


# separate hisat-3n unmapped reads
rule separate_unmapped_reads:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unsort.bam"),
    output:
        unique_bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.unique_aligned.bam")),
        multi_bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.multi_aligned.bam")),
        unmapped_fastq=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.unmapped.fastq"))
    threads:
        1
    run:
        separate_unique_and_multi_align_reads(in_bam_path=input.bam,
                                              out_unique_path=output.unique_bam,
                                              out_multi_path=output.multi_bam,
                                              out_unmappable_path=output.unmapped_fastq,
                                              unmappable_format='fastq',
                                              mapq_cutoff=10,
                                              qlen_cutoff=config['min_read_length'])


# split unmapped reads
rule split_unmapped_reads:
    input:
        unmapped_reads=local(bam_dir+"/{cell_id}.hisat3n_dna.unmapped.fastq")
    output:
        split_r1=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.R1.fastq")),
        split_r2=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.R2.fastq"))
    params:
        output_prefix=lambda wildcards: bam_dir+f"/{wildcards.cell_id}.hisat3n_dna.split_reads"
    threads:
        1
    run:
        split_hisat3n_unmapped_reads(fastq_path=input.unmapped_reads,
                                     output_prefix=params.output_prefix,
                                     min_length=config['min_read_length'])


# remap the split reads in SE mode
# Aligned reads FLAG and MAPQ possibilities:
# - [0, 60], uniquely mapped to forward strand
# - [16, 60], uniquely mapped to reverse strand
rule hisat_3n_single_end_mapping_dna_mode:
    input:
        fastq=local(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.{read_type}.fastq"), #"bam/{cell_id}.hisat3n_dna.split_reads.R1.fastq"
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.{read_type}.bam")),
        stats="bam/{cell_id}.hisat3n_dna_split_reads_summary.{read_type}.txt"
    params:
        direction=lambda wildcards: "--directional-mapping-reverse " if wildcards.read_type=="R1" else "--directional-mapping "
    threads:
        config['hisat3n_threads']
    resources:
        mem_mb=14000
    # benchmark:
    #         "fastq/{cell_id}-{read_type}.hisat_3n_single_end_mapping_dna_mode.benchmark.txt"
    shell:
        """
        hisat-3n {config[hisat3n_dna_reference]} -q -U {input.fastq} \
{params.direction} --base-change C,T {repeat_index_flag} \
--no-spliced-alignment --no-temp-splicesite -t --new-summary --summary-file {output.stats} \
--threads {threads} | samtools view -b -q 10 -o {output.bam}
        """
        # # only take the unique aligned reads

# sort split reads bam file by read name
rule merge_and_sort_split_reads_by_name:
    input:
        r1_bam=local(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.R1.bam"),
        r2_bam=local(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.R2.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.name_sort.bam"))
    threads:
        1
    shell:
        """
        samtools merge -o - {input.r1_bam} {input.r2_bam} | samtools sort -n -o {output.bam} -
        """


# remove overlap read parts from the split alignment bam file
rule remove_overlap_read_parts:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.name_sort.bam") #rules.merge_and_sort_split_reads_by_name.output.bam
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.no_overlap.bam"))
    threads:
        1
    run:
        remove_overlap_read_parts(in_bam_path=input.bam, out_bam_path=output.bam)


# merge all mapped reads
rule merge_original_and_split_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.unique_aligned.bam"),
        split_bam=local(bam_dir+"/{cell_id}.hisat3n_dna.split_reads.no_overlap.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.bam"))
    threads:
        1
    shell:
        """
        samtools merge -f {output.bam} {input.bam} {input.split_bam}
        """


# sort split reads bam file by read name
rule sort_all_reads_by_name:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam" #do not add local, upload to remote
    threads:
        1
    shell:
        """
        samtools sort -n -o {output.bam} {input.bam}
        """

# remove overlap parts and call contacts
rule call_chromatin_contacts:
    input:
        bam="bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam"
    output:
        stats="hic/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv",
        contact_tsv="hic/{cell_id}.hisat3n_dna.all_reads.3C.contact.tsv.gz",
        ded_contact="hic/{cell_id}.hisat3n_dna.all_reads.dedup_contacts.tsv.gz"
    params:
        contact_prefix=lambda wildcards: hic_dir+f"/{wildcards.cell_id}.hisat3n_dna.all_reads",
    threads:
        1
    run:
        if not os.path.exists(hic_dir):
            os.mkdir(hic_dir)
        call_chromatin_contacts(bam_path=input.bam,
                                contact_prefix=params.contact_prefix,
                                save_raw=False,
                                save_hic_format=True)


rule sort_bam_by_pos:
    input:
        bam="bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam"
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.pos_sort.bam"))
    resources:
        mem_mb=1000
    # benchmark:
    #         "fastq/{cell_id}.sort_bam_by_pos.benchmark.txt"
    threads:
        1
    shell:
        """
        samtools sort -O BAM -o {output.bam} {input.bam}
        """

# remove PCR duplicates
rule dedup:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.pos_sort.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam")), #to keep this bam, change to: "bam/{cell_id}.hisat3n_dna.all_reads.deduped.bam",
        stats="bam/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt"
    resources:
        mem_mb=3000
    # benchmark:
    #         "fastq/{cell_id}.dedup.benchmark.txt"
    threads:
        2
    shell:
        """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.stats} -REMOVE_DUPLICATES true -TMP_DIR bam/temp/
        """

# index the bam file
rule index_bam:
    input:
        bam=local(bam_dir+"/{input_name}.bam")
    output:
        bai=local(temp(bam_dir+"/{input_name}.bam.bai"))
    shell:
        """
        samtools index {input.bam}
        """
# ==================================================
# Generate ALLC
# ==================================================
rule unique_reads_allc:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam"),
        bai=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai")
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        tbi="allc/{cell_id}.allc.tsv.gz.tbi",
        stats="allc/{cell_id}.allc.tsv.gz.count.csv"
    threads:
        1.5
    resources:
        mem_mb=500
    # benchmark:
    #         "fastq/{cell_id}.unique_reads_allc.benchmark.txt"
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

# CGN extraction from ALLC
rule unique_reads_cgn_extraction:
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
--mc_contexts {mcg_context} --chrom_size_path {config[chrom_size_path]}
        """