import os,sys
import pandas as pd
import pathlib
import re
import glob
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
from cemba_data.gcp import *
from cemba_data.demultiplex import _parse_index_fasta,_read_cutadapt_result

default_config={
    'total_read_pairs_min':1,
    'total_read_pairs_max':10000000
}

# demultiplex can not be ran using spot mode, because in cutadapt step,
# {dir}/{uid}/lanes/{uid}-{lane}-{name}-R1.fq.gz the name is unknown, so
# it can not be upload onto cloud, if ran with spot, those files would lost.

if 'gcp' in config and config['gcp']:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
    fq_dir=config["fq_dir"]
    run_on_gcp=True
else:
    fq_dir=pathlib.Path(config["fq_dir"]).absolute()
    run_on_gcp=False

for key in default_config:
    if key not in config:
        config[key]=default_config[key]

outdir=config["outdir"] if 'outdir' in config else 'mapping'
local_outdir=outdir if not run_on_gcp else workflow.default_remote_prefix+"/"+outdir
barcode_version = config["barcode_version"] if 'barcode_version' in config else "V2"

df_fq=get_fastq_info(fq_dir,run_on_gcp)
# df_fq['stats_out']=df_fq.apply(lambda row: os.path.join(local_outdir, f"{row.uid}/demultiplex.stats.txt"),axis=1) #"{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt"
#for each pool, there are 16 old uid, 96 (16*6) new uid (6 multiplex groups)
uid_fastqs_dict=df_fq.loc[:,['uid','R1','R2']].groupby('uid').agg(lambda x:x.tolist()).to_dict(orient='index') #uid_fastqs_dict[uid]['R1'] and R2 are list

if barcode_version == 'V2' and df_fq['multiplex_group'].nunique() == 1:
    # print('Detect only single multiplex group in each plate, will use V2-single mode.')
    barcode_version = 'V2-single'

df_index=get_random_index(df_fq.uid.unique().tolist(),barcode_version) #old uids, 16 each pool

rule summary_demultiplex:
    input:
        stats=local(expand(local_outdir+"/{uid}/demultiplex.stats.txt",
                                    uid=df_index.old_uid.unique().tolist())) #old uids, 16
    output:
        csv=outdir+"/stats/demultiplex.stats.csv",
        fq_info=outdir+"/stats/fastq_info.tsv",
        index_info=outdir+"/stats/index_info.tsv"
    params:
        stat_dir=os.path.join(local_outdir,"stats")
    run:
        shell(f"mkdir -p {params.stat_dir}")
        df_fq.to_csv(output.fq_info,sep='\t',index=False)
        df_index.to_csv(output.index_info,sep='\t',index=False)
        # pathlib.Path(params.stat_dir).mkdir(exist_ok=True)
        random_index_fasta_path=os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') if barcode_version=='V1' else os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa')
        index_seq_dict = _parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        stat_list = []
        for path in glob.glob(local_outdir + "/*/demultiplex.stats.txt"):
            uid = path.split('/')[-2]
            single_df=_read_cutadapt_result(path)
            single_df['uid'] = uid
            single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
            assert single_df['index_name'].isna().sum() == 0
            stat_list.append(single_df)
        df_stats = pd.concat(stat_list)
        df_stats['multiplex_group']=df_stats['uid'].apply(lambda x:x.split('-')[1])
        df_stats['real_multiplex_group']=df_stats.index_name.apply(lambda x:((int(x[1:])-1) % 12) // 2 + 1 if 'unknow' not in x.lower() else 'NA')
        df_stats=df_stats.loc[df_stats.real_multiplex_group !='NA']
        df_stats['plate']=df_stats['uid'].apply(lambda x:x.split('-')[0])
        df_stats['primer_name']=df_stats['uid'].apply(lambda x:x.split('-')[-1])
        df_stats['uid']= df_stats.plate.map(str)+'-'+df_stats.real_multiplex_group.map(str)+'-'+df_stats.primer_name.map(str)
        df_stats['cell_id'] = df_stats['uid'] + '-' + df_stats['index_name']
        df_cell=df_stats.groupby('cell_id').agg({
                                        'Trimmed':'sum',
                                        'TotalPair':'sum',
                                        'index_name':lambda i: i.unique()[0],
                                        'uid':lambda i: i.unique()[0]})
        df_cell.rename(columns={'Trimmed': 'CellInputReadPairs',
                                'TotalPair': 'MultiplexedTotalReadPairs',
                                'index_name': 'IndexName',
                                'uid': 'UID'},inplace=True)
        df_cell['CellBarcodeRate'] = df_cell['CellInputReadPairs'] / df_cell['MultiplexedTotalReadPairs']
        df_cell['BarcodeVersion'] = barcode_version
        # print(type(output),output)
        df_cell.to_csv(output.csv)

rule download_from_gcp:
    output:
        fq=local("download/{path}.gz") #use download prefix to avoid conflict with merge_lanes
    run:
        dirname=os.path.dirname(output.fq)
        if not os.path.exists(dirname):
            os.makedirs(dirname,exist_ok=True)
        url='gs://'+wildcards.path+'.gz'
        if not os.path.exists(output.fq):
            os.system(f"gsutil -m cp -n {url} {dirname}")

rule merge_lanes:
    input:
        fqs=lambda wildcards: [local("download/"+fq.replace('gs://','')) for fq in uid_fastqs_dict[wildcards.uid][wildcards.read_type]] if run_on_gcp \
                                        else uid_fastqs_dict[wildcards.uid][wildcards.read_type]

    output:
        fq=local(temp(local_outdir+"/{uid}/{read_type}.fq.gz"))

    params:
        outdir=lambda wildcards: local_outdir+f"/{wildcards.uid}"

    run:
        shell(f"mkdir -p {params.outdir}")
        shell(f"cat {input.fqs} > {output.fq}")
        if run_on_gcp:
            for fq in input.fqs:
                print(f"Removing temporary raw fastq: {fq}")
                os.remove(fq)

rule run_demultiplex: #{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}.fastq.gz
    input: #uid = {plate}-{multiplex_group}-{primer_name} # primer_name is pcr index?
        R1 = lambda wildcards: local(local_outdir+f"/{wildcards.uid}/R1.fq.gz"),
        R2 = lambda wildcards: local(local_outdir+f"/{wildcards.uid}/R2.fq.gz")

    output: #uid, lane, index_name, read_type; dynamic: https://stackoverflow.com/questions/52598637/unknown-output-in-snakemake
        stats_out =local(local_outdir+"/{uid}/demultiplex.stats.txt"), # old_uid

    params:
        random_index_fa=lambda wildcards: os.path.join(PACKAGE_DIR, 'files', 'random_index_v1.fa') if barcode_version == "V1" \
		                            else os.path.join(PACKAGE_DIR, 'files', 'random_index_v2', 'random_index_v2.multiplex_group_' + wildcards.uid.split('-')[-2] + '.fa') if barcode_version == "V2" \
		                            else os.path.join(PACKAGE_DIR, 'files', 'random_index_v2', 'random_index_v2.fa'),
        outdir=lambda wildcards: local_outdir+f"/{wildcards.uid}/demultiplex",
        R1=lambda wildcards: local_outdir+f"/{wildcards.uid}/demultiplex/{'{{name}}'}-R1.fq.gz",
        R2=lambda wildcards: local_outdir+f"/{wildcards.uid}/demultiplex/{'{{name}}'}-R2.fq.gz"

    run:
        # print(params.R1,params.R2)
        shell(f"mkdir -p {params.outdir}")
        shell(f"cutadapt -Z -e 0.01 --no-indels -g file:{params.random_index_fa} -o  {params.R1} -p {params.R2} {input.R1} {input.R2} > {output.stats_out}")
         # for the reads startswith random index present in random_index_fa, will be taken and write into 1 fastq (1 cell),
         # cut the left 8 bp sequence and add the random index name (A2, P24) into the cell fastq name.
         # one uid will be broken down into 384 cells (if the number of multiplex group = 1: V2-single).

        # remote temporary fastq files
        os.remove(input.R1)
        os.remove(input.R2)

        #parse demultiplex.stats.txt for cell_qc
        random_index_fasta_path = os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') if barcode_version == 'V1' else os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa')
        index_seq_dict = _parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        single_df = _read_cutadapt_result(output.stats_out) #path="bican/salk010/UWA7648_CX182024_Idg_1_P1-1-K15/demultiplex.stats.txt"
        single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
        removed_index_names=single_df.loc[(single_df['Trimmed'] < int(config['total_read_pairs_min'])) | (single_df['Trimmed'] > int(config['total_read_pairs_max']))].index_name.unique().tolist()

         # rename & upload to GCP
        uids=wildcards.uid.split('-')
        input_fqs=glob.glob(local_outdir+f"/{wildcards.uid}/demultiplex/*-R*.fq.gz")
        for input_fq in input_fqs:
            index_name=os.path.basename(input_fq).split('-')[0]
            if index_name in removed_index_names: #cell qc
                print(f"{input_fq} was removed because Trimmed reads too large or too small")
                continue
            real_multiplex_group=index_name2multiplex_group(index_name)
            if real_multiplex_group=='NA':
                continue
            new_uid='-'.join([uids[0],str(real_multiplex_group),uids[-1]])
            read_type = os.path.basename(input_fq).rstrip('.fq.gz').split('-')[-1]
            new_path=local_outdir+f"/{new_uid}/fastq/{new_uid}-{index_name}-{read_type}.fq.gz"
            # print(input_fq,new_path)
            if run_on_gcp: #upload to GCP
                os.system(f"gsutil cp -n {input_fq} gs://{new_path}")
                os.remove(input_fq)
            else: #local, rename
                if not os.path.exists(os.path.dirname(new_path)):
                    os.makedirs(os.path.dirname(new_path),exist_ok=True)
                os.rename(input_fq,new_path)