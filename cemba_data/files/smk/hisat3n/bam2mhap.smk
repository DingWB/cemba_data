import os,sys
import glob

gcp=config.get("gcp",False)
bam_dir=os.path.abspath(workflow.default_remote_prefix+"/bam") if gcp else "bam"
mhap_dir=os.path.abspath(workflow.default_remote_prefix+"/mhap") if gcp else "mhap"

suffix=config.get("suffix",".hisat3n_dna.all_reads.name_sort.bam")
CELL_IDS=[os.path.basename(infile).replace(suffix,'') for infile in glob.glob(f"{bam_dir}/*{suffix}")]
print(len(CELL_IDS))

annotation_dir=os.path.expanduser(config['annotation_dir'])
annotation_suffix=config.get("annotation_suffix","hg38_allc.gz")

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
    params:
        annotation_path=lambda wildcards: os.path.join(annotation_dir,f"{wildcards.cell_id.split('_')[0]}.{annotation_suffix}")
    resources:
        mem_mb=400
    run:
        from cemba_data.mapping.pipelines import bam2mhap
        if not os.path.exists(mhap_dir):
            os.mkdir(mhap_dir)
        outfile1=output.mhap1[:-3] #"allc/{cell_id}.mhap", will be bgzipped and tabix indexed in mhap
        bam2mhap(bam_path=input.bam,annotation=params.annotation_path,
            output=outfile1,pattern="CGN")
        outfile2 = output.mhap2[:-3]
        bam2mhap(bam_path=input.bam,annotation=params.annotation_path,
            output=outfile2,pattern="CHN")