import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n_base.smk')

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

        # allc & mhap
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz",cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi",cell_id=CELL_IDS),
        expand("allc/{cell_id}.mhap.gz", cell_id=CELL_IDS),
        expand("allc/{cell_id}.mhap.gz.tbi",cell_id=CELL_IDS),

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
        snm3c_summary(outname=output.csv,indir=indir)

        # cleanup
        shell(f"rm -rf {bam_dir}/temp")

# Convert bam to mhap
rule bam_to_mhap:
    input: #sorted bam
            bam="bam/{cell_id}.hisat3n_dna.all_reads.deduped.bam",
            bai="bam/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai"
    output:
        mhap="allc/{cell_id}.mhap.gz",
        tbi="allc/{cell_id}.mhap.gz.tbi"
    params:
        cpgPath=os.path.expanduser(config['cpg_path']),
    resources:
        mem_mb=500
    run:
        from pym3c import bam2mhap
        outfile=output.mhap[:-3] #"allc/{cell_id}.mhap", will be bgzipped and tabix indexed in mhap
        bam2mhap(bam_path=input.bam,cpg_path=params.cpgPath,output=outfile)
