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

# ==================================================
# Methylation output format (allc / cz) and optional mhap generation
# The bare variables mc_format / reference_cz / generate_mhap / annotation_path
# are written at the top of the generated Snakefile by the *_config_str helpers
# (see cemba_data/mapping/pipelines/{m3c,mc}.py). Fall back to defaults for
# modes (mct / 4m) that do not define them.
# ==================================================
mhap_dir=os.path.abspath(workflow.default_remote_prefix+"/mhap") if config["gcp"] else "mhap"
cz_dir=os.path.abspath(workflow.default_remote_prefix+"/cz") if config["gcp"] else "cz"

def _coerce_bool(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ('true', '1', 'yes', 'y', 't', 'on')

mc_format = str(globals().get('mc_format', 'allc')).lower()
if mc_format not in ('allc', 'cz', 'both'):
    raise ValueError(
        f"Unknown mc_format {mc_format!r}, choose from 'allc', 'cz', 'both'")
generate_mhap = _coerce_bool(globals().get('generate_mhap', False))
reference_cz = globals().get('reference_cz', None)
# Whether to generate the CGN-merged ALLC (allc-CGN). Off by default.
extract_mcg = _coerce_bool(globals().get('extract_mcg', False))

# When .cz output is requested, a reference .cz (built by `czip build_ref`)
# is required. Warn and show how to build it if it is missing.
if mc_format in ('cz', 'both'):
    _genome_fasta = globals().get('reference_fasta', 'GENOME.fa')
    _chrom_size = globals().get('chrom_size_path', 'CHROM.sizes')
    _default_ref_cz = os.path.splitext(
        os.path.expanduser(str(_genome_fasta)))[0] + '.allc.cz'
    if reference_cz in (None, '', 'None'):
        sys.stderr.write(
            "\n[WARNING] mc_format=%r requires a reference .cz file but "
            "'reference_cz' is not set in the mapping config.\n"
            "Build one from the genome fasta and chrom_size with:\n"
            "    czip build_ref -g %s -O %s -s %s -j 20\n"
            "then set 'reference_cz = %s' in the mapping config.\n\n"
            % (mc_format, _genome_fasta, _default_ref_cz,
               _chrom_size, _default_ref_cz))
    elif not os.path.exists(os.path.expanduser(str(reference_cz))):
        sys.stderr.write(
            "\n[WARNING] reference_cz %r does not exist. Build it with:\n"
            "    czip build_ref -g %s -O %s -s %s -j 20\n\n"
            % (reference_cz, _genome_fasta,
               os.path.expanduser(str(reference_cz)), _chrom_size))


def get_methylation_targets(cell_ids):
    """Methylation output targets for the summary rule, based on mc_format."""
    targets = []
    if mc_format in ('allc', 'both'):
        targets += expand("allc/{cell_id}.allc.tsv.gz", cell_id=cell_ids)
        targets += expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=cell_ids)
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

def get_fastq_path():
    if config["fastq_server"]=='ftp':
        # FTP.remote("ftp.sra.ebi.ac.uk/vol1/fastq/SRR243/010/SRR24316310/SRR24316310_1.fastq.gz", keep_local=True)
        return lambda wildcards: FTP.remote(cell_dict[tuple([wildcards.cell_id,wildcards.read_type])])
    elif config["fastq_server"]=='gcp':
        return GS.remote("gs://" + workflow.default_remote_prefix + "/fastq/{cell_id}-{read_type}.fq.gz")
    else: # local
        return local("fastq/{cell_id}-{read_type}.fq.gz")