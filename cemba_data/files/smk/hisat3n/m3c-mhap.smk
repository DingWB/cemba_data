import os,sys
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n_base.smk')

include:
    os.path.join(PACKAGE_DIR,"files","smk",'hisat3n.smk')

mhap_dir=config['mhap_dir']
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

        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz",cell_id=CELL_IDS),
        expand("allc/{cell_id}.allc.tsv.gz.tbi",cell_id=CELL_IDS),

        # mhap
        expand("mhap/{cell_id}.mhap.gz", cell_id=CELL_IDS),
        expand("mhap/{cell_id}.mhap.gz.tbi",cell_id=CELL_IDS),
        expand("mhap/{cell_id}.mhap.gene.stat.tsv.gz", cell_id=CELL_IDS),
        expand("mhap/{cell_id}.mhap.promoter.stat.tsv.gz",cell_id=CELL_IDS),

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
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam"),
        bai=local(bam_dir + "/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai")
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

rule stat_mhap_gene:
    input: #sorted .mhap.gz
        mhap="mhap/{cell_id}.mhap.gz",
        tbi="mhap/{cell_id}.mhap.gz.tbi"
    output:
        gene_stat="mhap/{cell_id}.mhap.gene.stat.tsv.gz",
    params:
        cpgPath=os.path.expanduser(config['cpg_path']),
        geneBedPath=os.path.expanduser(config['gene_bed_path']),
    resources:
        mem_mb=400
    run:
        from cemba_data.mapping.pipelines import stat_mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        stat_mhap(mhap_path=input.mhap,cpg_path=params.cpgPath,
                region=None,bed=params.geneBedPath,output=output.gene_stat)

rule stat_mhap_promoter:
    input: #sorted .mhap.gz
        mhap="mhap/{cell_id}.mhap.gz",
        tbi="mhap/{cell_id}.mhap.gz.tbi"
    output:
        promoter_stat="mhap/{cell_id}.mhap.promoter.stat.tsv.gz",
    params:
        cpgPath=os.path.expanduser(config['cpg_path']),
        promoterBedPath=os.path.expanduser(config['promoter_bed_path']),
    resources:
        mem_mb=400
    run:
        from cemba_data.mapping.pipelines import stat_mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        stat_mhap(mhap_path=input.mhap,cpg_path=params.cpgPath,
                region=None,bed=params.promoterBedPath,output=output.promoter_stat)