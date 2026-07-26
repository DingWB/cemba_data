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
    # Methylation output format: 'allc' (ALLCools bam-to-allc),
    # 'cz' (cytozip bam_to_cz, default), or 'both'.
    'mc_format': 'cz',
    # Reference .cz file (built by `czip build_ref`), required when mc_format
    # is 'cz' or 'both'.
    'reference_cz': None,
    # Whether to also generate per-cell .mhap.gz files (bam -> mhap).
    'generate_mhap': False,
    # Annotation *_allc.gz path used by bam2mhap, required when generate_mhap.
    'annotation_path': None,
    # Whether to generate the CGN-merged ALLC (allc-CGN/*.CGN-Merge.allc.tsv.gz).
    # Off by default; set extract_mcg = True to produce it.
    'extract_mcg': False,
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
cz_dir=os.path.abspath(workflow.default_remote_prefix+"/cz") if config["gcp"] else "cz"

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

# ==================================================
# Methylation output format (allc / cz) and optional mhap generation
# ==================================================
def _coerce_bool(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ('true', '1', 'yes', 'y', 't', 'on')

mc_format = str(config.get('mc_format', 'allc')).lower()
if mc_format not in ('allc', 'cz', 'both'):
    raise ValueError(
        f"Unknown mc_format {mc_format!r}, choose from 'allc', 'cz', 'both'")
config['mc_format'] = mc_format

generate_mhap = _coerce_bool(config.get('generate_mhap', False))
config['generate_mhap'] = generate_mhap

# Whether to generate the CGN-merged ALLC (allc-CGN). Off by default.
extract_mcg = _coerce_bool(config.get('extract_mcg', False))
config['extract_mcg'] = extract_mcg

# When .cz output is requested, a reference .cz (built by `czip build_ref`)
# is required. Warn (and show how to build it) if it is missing so the user
# can generate it before/while running the pipeline.
if mc_format in ('cz', 'both'):
    reference_cz = config.get('reference_cz', None)
    genome_fasta = config.get('reference_fasta', 'GENOME.fa')
    chrom_size_path = config.get('chrom_size_path', 'CHROM.sizes')
    default_ref_cz = os.path.splitext(
        os.path.expanduser(str(genome_fasta)))[0] + '.allc.cz'
    if reference_cz in (None, '', 'None'):
        sys.stderr.write(
            "\n[WARNING] mc_format=%r requires a reference .cz file but "
            "'reference_cz' is not set in the mapping config.\n"
            "Build one from the genome fasta and chrom_size with:\n"
            "    czip build_ref -g %s -O %s -s %s -j 20\n"
            "then set 'reference_cz = %s' in the mapping config.\n\n"
            % (mc_format, genome_fasta, default_ref_cz,
               chrom_size_path, default_ref_cz))
    elif not os.path.exists(os.path.expanduser(str(reference_cz))):
        sys.stderr.write(
            "\n[WARNING] reference_cz %r does not exist. Build it with:\n"
            "    czip build_ref -g %s -O %s -s %s -j 20\n\n"
            % (reference_cz, genome_fasta,
               os.path.expanduser(str(reference_cz)), chrom_size_path))


def get_methylation_targets(cell_ids):
    """Methylation output targets for the summary rule, based on mc_format."""
    targets = []
    if mc_format in ('allc', 'both'):
        targets += expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=cell_ids)
        targets += expand("allc/{cell_id}.allc.tsv.gz", cell_id=cell_ids)
        targets += expand("allc/{cell_id}.allc.tsv.gz.tbi", cell_id=cell_ids)
        targets += get_mcg_targets(cell_ids)
    if mc_format in ('cz', 'both'):
        targets += expand("cz/{cell_id}.cz", cell_id=cell_ids)
    return targets


def get_mcg_targets(cell_ids):
    """CGN-merged ALLC targets, only when extract_mcg is enabled."""
    if not extract_mcg:
        return []
    return (expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",
                   cell_id=cell_ids, mcg_context=mcg_context)
            + expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",
                     cell_id=cell_ids, mcg_context=mcg_context))


def get_mhap_targets(cell_ids):
    """mhap output targets for the summary rule, enabled by generate_mhap."""
    if not generate_mhap:
        return []
    return (expand("mhap/{cell_id}.CG.mhap.gz", cell_id=cell_ids)
            + expand("mhap/{cell_id}.CG.mhap.gz.tbi", cell_id=cell_ids)
            + expand("mhap/{cell_id}.CH.mhap.gz", cell_id=cell_ids)
            + expand("mhap/{cell_id}.CH.mhap.gz.tbi", cell_id=cell_ids))

# print(f"bam_dir: {os.path.abspath(bam_dir)}")
# print(f"allc_dir: {os.path.abspath(allc_dir)}")

if config["fastq_server"]=='gcp' or config["gcp"]:
    print("gcp")
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
elif config["fastq_server"]=='ftp':
    # print("ftp")
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