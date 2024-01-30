"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import yaml
import pathlib
from cemba_data.hisat3n import *

# if "gcp" in config:
#     gcp=config["gcp"] # if the fastq files stored in GCP cloud, set gcp=True in snakemake: --config gcp=True
# else:
#     gcp=False
#
# if gcp:
#     from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
#     GS = GSRemoteProvider()
#     os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
#
# bam_dir=workflow.default_remote_prefix+"/bam" if gcp else "bam"
# allc_dir=workflow.default_remote_prefix+"/allc" if gcp else "allc"
# hic_dir=workflow.default_remote_prefix+"/hic" if gcp else "hic"
# ==================================================
# Preparation
# ==================================================
# read mapping config and put all variables into the locals()
DEFAULT_CONFIG = {
    'hisat3n_repeat_index_type': 'repeat',
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

if "local_fastq" in config:
    local_fastq=config["local_fastq"] # if the fastq files stored in GCP cloud, set local_fastq=False in snakemake: --config local_fastq=False
else:
    local_fastq=True

if not local_fastq or gcp:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')

bam_dir=workflow.default_remote_prefix+"/bam" if gcp else "bam"
allc_dir=workflow.default_remote_prefix+"/allc" if gcp else "allc"
hic_dir=workflow.default_remote_prefix+"/hic" if gcp else "hic"

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

# fastq table and cell IDs
# fastq_table = validate_cwd_fastq_paths()
# CELL_IDS = fastq_table.index.tolist()

mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"
allc_mcg_dir=workflow.default_remote_prefix+f"allc-{mcg_context}" if gcp else f"allc-{mcg_context}"
allc_multi_dir=workflow.default_remote_prefix+"allc-multi" if gcp else "allc-multi"

module m3c:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        "https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/hisat3n/snakefile/m3c.Snakefile"

use rule * from m3c exclude summary

# ==================================================
# Mapping summary
# ==================================================
rule summary:
    input:
        # fastq trim
        local(expand("fastq/{cell_id}.trimmed.stats.tsv",
                        cell_id=CELL_IDS)),
        # dna mapping
        local(expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS)),

        # bam dir
        local(expand(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt",cell_id=CELL_IDS)),
        local(expand(bam_dir+"/{cell_id}.hisat3n_dna_split_reads_summary.{read_type}.txt",
                        cell_id=CELL_IDS,read_type=['R1','R2'])),

        # 3C contacts
        local(expand(hic_dir+"/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv"), cell_id=CELL_IDS),
        # allc
#         expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        local(expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS)),
        local(expand("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv",cell_id=CELL_IDS)),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",cell_id=CELL_IDS, mcg_context=mcg_context),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",
               cell_id=CELL_IDS, mcg_context=mcg_context),
    output:
        "MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        # generate the final summary
        snm3c_summary()

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

# include everything in m3c.Snakefile except summary (rewrite summary)
#=====================================================================
# Processing multi-alignment reads and generate allc-multi
#=====================================================================
rule sort_multi_bam:
    input:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.multi_aligned.bam")), #"bam/{cell_id}.hisat3n_dna.multi_aligned.bam"
    output:
        bam=local(temp(bam_dir+"/{cell_id}.hisat3n_dna_sorted.multi_align.bam"))
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        "samtools sort -O BAM -o {output.bam} {input.bam}"


rule dedup_multi_bam:
    input:
        bam=rules.sort_multi_bam.output.bam #"bam/{cell_id}.hisat3n_dna_sorted.multi_align.bam"
    output:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        stats=local(temp(bam_dir+"/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt"))
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"

rule index_multi_bam_dna_reads:
    input:
        bam=rules.dedup_multi_bam.output.bam #"bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam"
    output:
        bai="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam.bai"
    shell:
        "samtools index {input.bam}"


# generate ALLC
rule multi_reads_allc:
    input:
        bam=rules.dedup_multi_bam.output.bam, #"bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        bai=rules.index_multi_bam_dna_reads.output.bai #"bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam.bai"
    output:
        allc="allc-multi/{cell_id}.allc_multi.tsv.gz",
        stats=local(temp("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv"))
    threads:
        1.5
    resources:
        mem_mb=500
    shell:
        """
        mkdir -p {allc_multi_dir}
        allcools bam-to-allc --bam_path {input.bam} \
            --reference_fasta {config[reference_fasta]} \
            --output_path {output.allc} \
            --num_upstr_bases {config[num_upstr_bases]} \
            --num_downstr_bases {config[num_downstr_bases]} \
            --compress_level {config[compress_level]} \
            --save_count_df \
            --min_mapq 0 # for multi-mapped reads, skip mapq filter \
            --convert_bam_strandness
        """

