import os,sys
import pandas as pd
import pathlib
import re
import cemba_data
PACKAGE_DIR=cemba_data.__path__[0]
from cemba_data.gcp import *
from cemba_data.demultiplex import _parse_index_fasta,_read_cutadapt_result

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
fq_ext=config["fq_ext"] if 'fq_ext' in config else 'fastq'
outdir=config["outdir"]) if 'outdir' in config else 'mapping'
barcode_version = config["barcode_version"] if 'barcode_version' in config else "V2"
env_name="yap" if 'env_name' not in config else config['env_name']

print(outdir)

df=get_fastq_info(fq_dir,outdir,fq_ext,run_on_gcp)

if barcode_version == 'V2' and df['multiplex_group'].nunique() == 1:
    print('Detect only single multiplex group in each plate, will use V2-single mode.')
    barcode_version = 'V2-single'

rule write_fastq_info:
    input:
        os.path.join(outdir,"stats/demultiplex.stats.csv")
    output:
        tsv=os.path.join(outdir,"stats/fastq_info.tsv")
    run:
        if os.path.exists("fastq_info.txt"):
            os.rename("fastq_info.txt",output.tsv)
        else:
            df.to_csv(output.tsv,sep='\t',index=False)


# rule demultiplex:
#     input:
#         os.path.join(outdir,"stats/fastq_info.tsv")

rule run_demultiplex: #{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}.fastq.gz
    input: #uid = {plate}-{multiplex_group}-{primer_name} # primer_name is pcr index?
        R1 = lambda wildcards: df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R1.iloc[0] if not run_on_gcp \
            else GS.remote(df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R1.iloc[0]),
        R2 = lambda wildcards: df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R2.iloc[0] if not run_on_gcp \
            else GS.remote(df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R2.iloc[0])

    output: #uid, lane, index_name, read_type; dynamic: https://stackoverflow.com/questions/52598637/unknown-output-in-snakemake
        stats_out ="{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt",
#         stats_out =dynamic("{dir}/{uid}/lanes/{uid}-{lane}-{name}.demultiplex.stats.txt"),
#         R1=dynamic("{dir}/{uid}/lanes/{uid}-{lane}-{name}-R1.fq.gz"),
#         R2=dynamic("{dir}/{uid}/lanes/{uid}-{lane}-{name}-R2.fq.gz"),
    conda:
        env_name

    params:
        random_index_fa=lambda wildcards: \
                        os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') \
                        if barcode_version=="V1" else \
                        os.path.join(PACKAGE_DIR,'files','random_index_v2',\
                        'random_index_v2.multiplex_group_'+wildcards.uid.split('-')[-2]+'.fa') \
                        if  barcode_version=="V2" else \
                        os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa'),
        outdir=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes", # will be copy to GCP in merge_lanes.smk
        outdir2=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes" if not run_on_gcp else \
                                                workflow.default_remote_prefix+f"/{wildcards.dir}/{wildcards.uid}/lanes",
        R1=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes/{wildcards.uid}-{wildcards.lane}-{{name}}-R1.fq.gz",
        R2=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes/{wildcards.uid}-{wildcards.lane}-{{name}}-R2.fq.gz"

    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.outdir2}
#         echo {output.stats_out}
        cutadapt -Z -e 0.01 --no-indels -g file:{params.random_index_fa} \
        -o  {params.R1} -p {params.R2} {input.R1} {input.R2} > {output.stats_out}
        """
         # for the reads startswith random index present in random_index_fa, will be taken and write into 1 fastq (1 cell),
         # cut the left 8 bp sequence and add the random index name (A2, P24) into the cell fastq name.
         # one uid will be broken down into 384 cells (if the number of multiplex group = 1: V2-single).
         # If run on GCP, the output R1 and R2 in the params are not really generated on GCP, cause it contains {name}, and
         # not included in output, so it will not be uploaded to GCP automatically, will be uploaded in merge_lanes.smk

rule summary_demultiplex:
    input:
        df['stats_out'].tolist() #{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt", one uid contains multiple lanes and index_names
    output:
        csv="{dir}/stats/demultiplex.stats.csv"
    params:
        stat_dir=lambda wildcards:os.path.join(wildcards.dir,"stats") if not run_on_gcp else \
                     os.path.join(workflow.default_remote_prefix,wildcards.dir,"stats")
    run:
#         print(params.stat_dir)
        shell(f"mkdir -p {params.stat_dir}")
        # pathlib.Path(params.stat_dir).mkdir(exist_ok=True)
        random_index_fasta_path=os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') if barcode_version=='V1' else \
                                os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa')
        index_seq_dict = _parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        stat_list = []
        for path in input:
            *uid, suffix = os.path.basename(path).split('-')
            lane = suffix.split('.')[0]
            uid = '-'.join(uid)
            single_df=_read_cutadapt_result(path)
            single_df['uid'] = uid
            single_df['lane'] = lane
            single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
            assert single_df['index_name'].isna().sum() == 0
            stat_list.append(single_df)
        df_stats = pd.concat(stat_list)
        df_stats['multiplex_group']=df_stats['uid'].apply(lambda x:x.split('-')[1])
        if barcode_version == 'V2' and df_stats['multiplex_group'].nunique() == 1:
            df_stats['real_multiplex_group']=df_stats.index_name.apply(lambda x:((int(x[1:])-1) % 12) // 2 + 1 if 'unknow' not in x.lower() else 'NA')
        else:
            df_stats['real_multiplex_group']=df_stats.multiplex_group.tolist()

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