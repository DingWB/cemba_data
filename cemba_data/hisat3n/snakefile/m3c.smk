"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import yaml
import pathlib
from cemba_data.hisat3n import *

# ==================================================
# Preparation
# ==================================================
# read mapping config and put all variables into the locals()
DEFAULT_CONFIG = {
    'hisat3n_repeat_index_type': '',
    'r1_adapter': 'AGATCGGAAGAGCACACGTCTGAAC',
    'r2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGA',
    'r1_right_cut': 10,
    'r2_right_cut': 10,
    'r1_left_cut': 10,
    'r2_left_cut': 10,
    'min_read_length': 30,
    'num_upstr_bases': 0,
    'num_downstr_bases': 2,
    'compress_level': 5,
    'hisat3n_threads': 11,
    # the post_mapping_script can be used to generate dataset, run other process etc.
    # it gets executed before the final summary function.
    # the default command is just a placeholder that has no effect
    'post_mapping_script': 'true',
}
REQUIRED_CONFIG = ['hisat3n_dna_reference', 'reference_fasta', 'chrom_size_path']

if "gcp" in config:
    gcp=config["gcp"] # if the fastq files stored in GCP cloud, set gcp=True in snakemake: --config gcp=True
else:
    gcp=False

if "local_fastq" in config and gcp:
    local_fastq=config["local_fastq"] # if the fastq files stored in GCP cloud, set local_fastq=False in snakemake: --config local_fastq=False
else:
    local_fastq=True

if not local_fastq or gcp:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')

bam_dir=os.path.abspath(workflow.default_remote_prefix+"/bam") if gcp else "bam"
allc_dir=os.path.abspath(workflow.default_remote_prefix+"/allc") if gcp else "allc"
hic_dir=os.path.abspath(workflow.default_remote_prefix+"/hic") if gcp else "hic"

local_config = read_mapping_config()
DEFAULT_CONFIG.update(local_config)

for k, v in DEFAULT_CONFIG.items():
    if k not in config:
        config[k] = v

missing_key = []
for k in REQUIRED_CONFIG:
    if k not in config:
        missing_key.append(k)
if len(missing_key) > 0:
    raise ValueError('Missing required config: {}'.format(missing_key))

if not gcp:
    # fastq table and cell IDs
    fastq_table = validate_cwd_fastq_paths()
    CELL_IDS = fastq_table.index.tolist() # CELL_IDS will be writen in the beginning of this snakemake file.

mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"
allc_mcg_dir=os.path.abspath(workflow.default_remote_prefix+f"/allc-{mcg_context}") if gcp else f"allc-{mcg_context}"
# print(f"bam_dir: {bam_dir}\n allc_dir: {allc_dir}\n hic_dir: {hic_dir} \n allc_mcg_dir: {allc_mcg_dir}")

for dir in [bam_dir,allc_dir,hic_dir,allc_mcg_dir]:
    if not os.path.exists(dir):
        os.mkdir(dir)

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
        expand("bam/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna_split_reads_summary.{read_type}.txt",
                        cell_id=CELL_IDS,read_type=['R1','R2']),
        expand("bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam", cell_id=CELL_IDS),

        # 3C contacts
        expand("hic/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv", cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.3C.contact.tsv.gz",cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.dedup_contacts.tsv.gz",cell_id=CELL_IDS),

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
        snm3c_summary(outname=output.csv)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")


# ==================================================
# FASTQ Trimming
# ==================================================


# Trim reads
# sort the fastq files so that R1 and R2 are in the same order
rule sort_fq:
    input:
        fq=local("fastq/{cell_id}-{read_type}.fq.gz") if local_fastq else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-{read_type}.fq.gz"),
    output:
        fq=local(temp("fastq/{cell_id}-{read_type}_sort.fq")),
    threads:
        1.5
    resources:
        high_io_job=1
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
    threads:
        1
    shell:
        "cutadapt "
        "-a R1Adapter={config[r1_adapter]} "
        "-A R2Adapter={config[r2_adapter]} "
        "--report=minimal "
        "-O 6 "
        "-q 20 "
        "-u {config[r1_left_cut]} "
        "-u -{config[r1_right_cut]} "
        "-U {config[r2_left_cut]} "
        "-U -{config[r2_right_cut]} "
        "-Z "
        "-m {config[min_read_length]}:{config[min_read_length]} "
        "--pair-filter 'both' "
        "-o {output.R1} "
        "-p {output.R2} "
        "{input.R1} {input.R2} "
        "> {output.stats}"


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
        mem_mb=14000
    shell:
        """
        hisat-3n {config[hisat3n_dna_reference]} -q  -1 {input.R1} -2 {input.R2} \
            --directional-mapping-reverse \
            --base-change C,T {repeat_index_flag} \
            --no-spliced-alignment \
            --no-temp-splicesite -t  --new-summary \
            --summary-file {output.stats} \
            --threads {threads} | samtools view -b -q 0 -o {output.bam}
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
    shell:
        "hisat-3n "
        "{config[hisat3n_dna_reference]} "
        "-q "
        "-U {input.fastq} "
        "--directional-mapping-reverse "  # map R1 in pbat mode
        "--base-change C,T "
        "{repeat_index_flag} "
        "--no-spliced-alignment "  # this is important for DNA mapping
        "--no-temp-splicesite "
        "-t "
        "--new-summary "
        "--summary-file {output.stats} "
        "--threads {threads} "
        "| "
        "samtools view "
        "-b -q 10 -o {output.bam}"  # only take the unique aligned reads

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
        "samtools merge -o - {input.r1_bam} {input.r2_bam} | samtools sort -n -o {output.bam} -"


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
        "samtools merge -f {output.bam} {input.bam} {input.split_bam}"


# sort split reads bam file by read name
rule sort_all_reads_by_name:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.bam")
    output:
        bam="bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam" #do not add local, upload to remote
    threads:
        1
    shell:
        "samtools sort -n -o {output.bam} {input.bam}"

# remove overlap parts and call contacts
rule call_chromatin_contacts:
    input:
        bam=rules.sort_all_reads_by_name.output.bam #"bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam"
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


rule sort_bam:
    input:
        bam=rules.sort_all_reads_by_name.output.bam #"bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam"
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.pos_sort.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        "samtools sort -O BAM -o {output.bam} {input.bam}"

# remove PCR duplicates
rule dedup_unique_bam:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.pos_sort.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam")),
        stats="bam/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt"
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        "picard MarkDuplicates I={input.bam} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"


# index the bam file
rule index_unique_bam_dna_reads:
    input:
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam")
    output:
        bai=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai"))
    shell:
        "samtools index {input.bam}"

# ==================================================
# Generate ALLC
# ==================================================
# generate ALLC
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
    shell:
        """
        mkdir -p {allc_dir}
        allcools bam-to-allc --bam_path {input.bam} \
            --reference_fasta {config[reference_fasta]} \
            --output_path {output.allc} \
            --num_upstr_bases {config[num_upstr_bases]} \
            --num_downstr_bases {config[num_downstr_bases]} \
            --compress_level {config[compress_level]} \
            --save_count_df \
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
            --allc_path  {input.allc} \
            --output_prefix {params.prefix} \
            --mc_contexts {mcg_context} \
            --chrom_size_path {config[chrom_size_path]}
        """