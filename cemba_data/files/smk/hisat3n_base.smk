"""
Snakemake pipeline for hisat-3n mapping of snm3C-seq data

hg38 normal index uses ~9 GB of memory
repeat index will use more memory
"""
import os,sys
import pandas as pd
import yaml
import pathlib
from cemba_data.hisat3n import *

# ==================================================
# Preparation
# ==================================================
# read mapping config and put all variables into the locals()
DEFAULT_CONFIG = {
    'hisat3n_repeat_index_type': 'no-repeat',
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

if "gcp" not in config:
    config["gcp"]=False #whether run on GCP (write output to GCP bucket)

if "fastq_server" not in config:
    config["fastq_server"]='local' # can be local, gcp, ftp

bam_dir=os.path.abspath(workflow.default_remote_prefix+"/bam") if config["gcp"] else "bam"
allc_dir=os.path.abspath(workflow.default_remote_prefix+"/allc") if config["gcp"] else "allc"
allc_multi_dir=os.path.abspath(workflow.default_remote_prefix+"/allc-multi") if config["gcp"] else "allc-multi"
hic_dir=os.path.abspath(workflow.default_remote_prefix+"/hic") if config["gcp"] else "hic"
mhap_dir=os.path.abspath(workflow.default_remote_prefix+"/mhap") if config["gcp"] else "mhap"

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

# if not config["gcp"]: # local
#     # fastq table and cell IDs
#     fastq_table = validate_cwd_fastq_paths() #get fastq path from pathlib.Path(f'{cwd}/fastq/').glob('*.[fq.gz][fastq.gz]')
#     CELL_IDS = fastq_table.index.tolist() # CELL_IDS will be writen in the beginning of this snakemake file.

mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
#repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"
repeat_index_flag="--no-repeat-index" #repeat would cause some randomness, get different output (mapping summary) even using the same input and parameters
allc_mcg_dir=os.path.abspath(workflow.default_remote_prefix+f"/allc-{mcg_context}") if config["gcp"] else f"allc-{mcg_context}"
# print(f"bam_dir: {bam_dir}\n allc_dir: {allc_dir}\n hic_dir: {hic_dir} \n allc_mcg_dir: {allc_mcg_dir}")

for dir in [bam_dir,allc_dir]:
    if not os.path.exists(dir):
        os.mkdir(dir)

print(f"bam_dir: {os.path.abspath(bam_dir)}")
print(f"allc_dir: {os.path.abspath(allc_dir)}")

if config["fastq_server"]=='gcp' or config["gcp"]:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
elif config["fastq_server"]=='ftp':
    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()
    fastq_dir = os.path.abspath(workflow.default_remote_prefix + "/fastq") if config["gcp"] else "fastq"
    os.makedirs(fastq_dir,exist_ok=True)
    # print(f"cwd: {os.getcwd()}") # the same as parameter: snakemake -d
    cell_id_path=os.path.abspath("gs://"+workflow.default_remote_prefix + "/CELL_IDS") if config["gcp"] else "CELL_IDS"
    # instead of creating fastq directory, there should be a file names CELL_IDS, columns: cell_id,read_type and fastq_path should be present
    df1=pd.read_csv(cell_id_path,sep='\t')
    # print(df1.head())
    cell_dict=df1.set_index(['cell_id','read_type']).fastq_path.to_dict()

def get_fastq_path():
    if config["fastq_server"]=='ftp':
        # FTP.remote("ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/010/SRR24316310/SRR24316310_1.fastq.gz", keep_local=True)
        return lambda wildcards: FTP.remote(cell_dict[tuple([wildcards.cell_id,wildcards.read_type])])
    elif config["fastq_server"]=='gcp':
        return GS.remote("gs://" + workflow.default_remote_prefix + "/fastq/{cell_id}-{read_type}.fq.gz")
    else: # local
        return local("fastq/{cell_id}-{read_type}.fq.gz")