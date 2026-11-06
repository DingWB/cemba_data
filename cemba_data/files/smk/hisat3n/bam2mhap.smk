import os,sys
import glob

gcp=config.get("gcp",False)
bam_dir=os.path.abspath(workflow.default_remote_prefix+"/bam") if gcp else "bam"
mhap_dir=os.path.abspath(workflow.default_remote_prefix+"/mhap") if gcp else "mhap"

suffix=config.get("suffix",".hisat3n_dna.all_reads.name_sort.bam")
if "CELL_IDS" not in globals(): # CELL_IDS is written at the top of the Snakefile by start_from_cell_bam
    CELL_IDS=[os.path.basename(infile).replace(suffix,'') for infile in glob.glob(f"{bam_dir}/*{suffix}")]
print(len(CELL_IDS))

# CpG / CH sites come from the genome fasta: --config reference_fasta=<genome.fa> [chrom_size_path=<sizes>],
# otherwise from the mapping_config.ini next to the Snakefile or the working directory (or their parent)
if "reference_fasta" not in config:
    from cemba_data.utilities import get_configuration
    for _d in [workflow.basedir, os.path.dirname(workflow.basedir), os.getcwd(), os.path.dirname(os.getcwd())]:
        if os.path.exists(os.path.join(_d, "mapping_config.ini")):
            _mapping_config = get_configuration(os.path.join(_d, "mapping_config.ini"))
            for _k in ["reference_fasta", "chrom_size_path", "mhap_cache_dir"]:
                if _k in _mapping_config and _k not in config:
                    config[_k] = _mapping_config[_k]
            break
if "reference_fasta" not in config:
    raise ValueError("bam2mhap.smk needs reference_fasta: --config reference_fasta=<genome fasta the BAMs were aligned to>, "
                     "or a mapping_config.ini in the working directory")
reference_fasta=os.path.expanduser(str(config["reference_fasta"]))
chrom_size_path=config.get("chrom_size_path",None)
# optional on-disk cache of the site flags: --config mhap_cache_dir=<dir> (or mhap_cache_dir in mapping_config.ini)
mhap_cache_dir=config.get("mhap_cache_dir",None)
mhap_cache_dir=None if mhap_cache_dir in (None,"","None") else os.path.expanduser(str(mhap_cache_dir))

rule summary:
    input:
        # mhap
        expand("mhap/{cell_id}.CG.mhap.gz", cell_id=CELL_IDS),
        expand("mhap/{cell_id}.CG.mhap.gz.tbi",cell_id=CELL_IDS),
        expand("mhap/{cell_id}.CH.mhap.gz",cell_id=CELL_IDS),
        expand("mhap/{cell_id}.CH.mhap.gz.tbi",cell_id=CELL_IDS),

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

# Convert bam to mhap
rule bam_to_mhap:
    input: #sorted bam
        bam=local(bam_dir+"/{cell_id}.hisat3n_dna.all_reads.deduped.bam"),
        bai=local(bam_dir + "/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai")
    output:
        mhap1="mhap/{cell_id}.CG.mhap.gz",
        tbi1="mhap/{cell_id}.CG.mhap.gz.tbi",
        mhap2="mhap/{cell_id}.CH.mhap.gz",
        tbi2="mhap/{cell_id}.CH.mhap.gz.tbi"
    resources:
        mem_mb=1500
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        bam2mhap(bam_path=input.bam,
                 reference_fasta=reference_fasta,
                 chroms=os.path.expanduser(chrom_size_path) if chrom_size_path else None,
                 outputs={'CGN': output.mhap1, 'CHN': output.mhap2},
                 cache_dir=mhap_cache_dir)
