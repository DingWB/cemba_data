"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import yaml
import pathlib
import pandas as pd

if "gcp" not in config:
    config["gcp"]=False #whether run on GCP (write output to GCP bucket)

if "fastq_server" not in config:
    config["fastq_server"]='local' # can be local, gcp, ftp

bam_dir=os.path.abspath(workflow.default_remote_prefix+"/bam") if config["gcp"] else "bam"
allc_dir=os.path.abspath(workflow.default_remote_prefix+"/allc") if config["gcp"] else "allc"
hic_dir=os.path.abspath(workflow.default_remote_prefix+"/hic") if config["gcp"] else "hic"
fastq_dir=os.path.abspath(workflow.default_remote_prefix+"/fastq") if config["gcp"] else "fastq"
mcg_context = 'CGN' if int(num_upstr_bases) == 0 else 'HCGN'
allc_mcg_dir=os.path.abspath(workflow.default_remote_prefix+f"/allc-{mcg_context}") if config["gcp"] else f"allc-{mcg_context}"

if config["fastq_server"]=='gcp' or config["gcp"]:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
elif config["fastq_server"]=='ftp':
    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()
    fastq_dir = os.path.abspath(workflow.default_remote_prefix + "/fastq") if config["gcp"] else "fastq"
    os.makedirs(fastq_dir,exist_ok=True)
    cell_id_path=os.path.abspath(workflow.default_remote_prefix + "/CELL_IDS") if config["gcp"] else "CELL_IDS"
    # instead of creating fastq directory, there should be a file names CELL_IDS, columns: cell_id,read_type and fastq_path should be present
    cell_dict=pd.read_csv(cell_id_path,sep='\t').set_index(['cell_id','read_type']).fastq_path.to_dict()

for dir in [bam_dir,allc_dir,hic_dir,allc_mcg_dir]:
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_fastq_path():
    if config["fastq_server"]=='ftp':
        # FTP.remote("ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/010/SRR24316310/SRR24316310_1.fastq.gz", keep_local=True)
        key=lambda wildcards: tuple([wildcards.cell_id,wildcards.read_type])
        print("ftp:",key,cell_dict[key])
        return FTP.remote(cell_dict[key])
    elif config["fastq_server"]=='gcp':
        return GS.remote("gs://" + workflow.default_remote_prefix + "/fastq/{cell_id}-{read_type}.fq.gz")
    else: # local
        return local("fastq/{cell_id}-{read_type}.fq.gz")