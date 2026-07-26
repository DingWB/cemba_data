import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')
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
        # expand("bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam", cell_id=CELL_IDS),

        # 3C contacts
        expand("hic/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv", cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.3C.contact.tsv.gz",cell_id=CELL_IDS),
        expand("hic/{cell_id}.hisat3n_dna.all_reads.dedup_contacts.tsv.gz",cell_id=CELL_IDS),

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
        snm3c_summary(outname=output.csv,indir=indir,mc_format=config['mc_format'])

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")
