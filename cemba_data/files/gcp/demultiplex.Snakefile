import os.path

import pandas as pd
import pathlib
import re
import cemba_data
# from snakemake.io import glob_wildcards
PACKAGE_DIR=cemba_data.__path__[0]

if 'gcp' in config and config["gcp"]:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
    fq_dir=config["fq_dir"]
else:
    fq_dir=pathlib.Path(config["fq_dir"]).absolute()
fq_ext=config["fq_ext"] if 'fq_ext' in config else 'fastq'
outdir=config["outdir"] if 'outdir' in config else 'mapping'
barcode_version = config["barcode_version"] if 'barcode_version' in config else "V2"

print(outdir)

def get_fastq_info(fq_dir,outdir):
    fq_df_path="fastq_info.txt"
    if not os.path.exists(fq_df_path):
        #For example: 220517-AMB-mm-na-snm3C_seq-NovaSeq-pe-150-WT-AMB_220510_8wk_12D_13B_2_P3-1-A11_S7_L001_R1_001.fastq.gz
        indirs,prefixes,plates,multiple_groups,primer_names,pns,lanes,\
        read_types,suffixes=glob_wildcards(os.path.join(str(fq_dir),\
                "{indir}/{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")) \
                if not config["gcp"] else \
                GS.glob_wildcards(fq_dir+"/{indir}/{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")
        indirs=[fq_dir + '/' + d for d in indirs]

        if len(indirs)==0: # depth=1
            prefixes,plates,multiple_groups,primer_names,pns,lanes,\
            read_types,suffixes=glob_wildcards(os.path.join(str(fq_dir),\
                    "{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")) \
                    if not config["gcp"] else \
                    GS.glob_wildcards(fq_dir+"/{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")
            indirs=[str(fq_dir)] * len(prefixes)

        df=pd.DataFrame.from_dict(
                                    {'indir':indirs,
                                    'prefix':prefixes,
                                    'plate':plates,
                                    'multiplex_group':multiple_groups,
                                    'primer_name':primer_names,
                                    'pns':pns,
                                    'lane':lanes,
                                    'read_type':read_types,
                                    'suffix':suffixes}
                                    )
        df['fastq_path']=df.apply(lambda row:os.path.join(row.indir,\
                                '-'.join(row.loc[['prefix','plate','multiplex_group','primer_name']].map(str).tolist())+\
                                "_"+"_".join(row.loc[['pns','lane','read_type','suffix']].map(str).tolist())+f".{fq_ext}.gz")
                                ,axis=1)
        df['uid']=df.plate.map(str)+'-'+df.multiplex_group.map(str)+'-'+df.primer_name.map(str) #f'{plate}-{multiplex_group}-{primer_name}')
        assert df.groupby(['lane','read_type'])['uid'].nunique().nunique()==1
        df=df.loc[df.read_type=='R1']
        df.rename(columns={'fastq_path':'R1'},inplace=True)
        df['R2']=df.R1.apply(lambda x:x.replace('_R1_','_R2_'))
        df['stats_out']=df.apply(lambda row: os.path.join(outdir, f"{row.uid}/lanes/{row.uid}-{row.lane}.demultiplex.stats.txt"),
                                            axis=1) #"{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt"
        df.to_csv(fq_df_path,sep='\t',index=False)

    df=pd.read_csv(fq_df_path,sep='\t')
    return df

def read_cutadapt_result(stat_path):
    """
    Parser of cutadapt output
    """
    with open(stat_path) as f:
        p = re.compile(
            r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times")
        series = []
        total_pairs = -1
        for line in f:
            if line.startswith('Total read pairs processed'):
                total_pairs = line.split(' ')[-1]
                total_pairs = int(total_pairs.replace(',','').replace(' ',''))

            m = p.search(line)
            if m is not None:
                result_dict = {}
                for i in m.group().split('; '):
                    k, v = i.split(': ')
                    result_dict[k] = v
                result_series = pd.Series(result_dict)
                series.append(result_series)
        total_df = pd.DataFrame(series)
        total_df['Trimmed'] = total_df['Trimmed'].apply(
            lambda c: c.split(' ')[0]).astype(int) # how many times the index_name appear in the beginning of read (meaning the number of reads in one cell)
        total_df['TotalPair'] = total_pairs
        total_df['Ratio'] = total_df['Trimmed'] / total_pairs
    return total_df

def parse_index_fasta(fasta_path):
    records = {}
    with open(fasta_path) as f:
        key_line = True
        for line in f:
            if key_line:
                key = line.lstrip('>').rstrip('\n')
                key_line = False
            else:
                value = line.lstrip('^').rstrip('\n')
                records[key] = value
                key_line = True
    return records

df=get_fastq_info(fq_dir,outdir)

if barcode_version == 'V2' and df['multiplex_group'].nunique() == 1:
    print('Detect only single multiplex group in each plate, will use V2-single mode.')
    barcode_version = 'V2-single'

rule write_fastq_info:
    input:
        os.path.join(outdir,"stats/demultiplex.stats.csv")
    output:
        tsv=os.path.join(outdir,"stats/fastq_info.tsv")
    run:
        df.to_csv(output.tsv,sep='\t',index=False)
        if os.path.exists("fastq_info.txt"):
            os.remove("fastq_info.txt")

# rule demultiplex:
#     input:
#         os.path.join(outdir,"stats/fastq_info.tsv")

rule run_demultiplex:
    input: #uid = {plate}-{multiplex_group}-{primer_name} # primer_name is pcr index?
        R1 = lambda wildcards: df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R1.iloc[0] if not config["gcp"] \
            else GS.remote(df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R1.iloc[0]),
        R2 = lambda wildcards: df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R2.iloc[0] if not config["gcp"] \
            else GS.remote(df.loc[(df.uid==wildcards.uid) & (df.lane==wildcards.lane)].R2.iloc[0])

    output: #uid, lane, index_name, read_type
        stats_out ="{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt",
    conda:
        "yap"

    params:
        random_index_fa=lambda wildcards: \
                        os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') \
                        if barcode_version=="V1" else \
                        os.path.join(PACKAGE_DIR,'files','random_index_v2',\
                        'random_index_v2.multiplex_group_'+wildcards.uid.split('-')[-2]+'.fa') \
                        if  barcode_version=="V2" else \
                        os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa'),
        outdir=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes", # will be copy to GCP in merge_lanes.smk
        outdir2=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes" if not config['gcp'] else \
                                                workflow.default_remote_prefix+f"/{wildcards.dir}/{wildcards.uid}/lanes",
        R1=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes/{wildcards.uid}-{wildcards.lane}-{{name}}-R1.fq.gz",
        R2=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/lanes/{wildcards.uid}-{wildcards.lane}-{{name}}-R2.fq.gz"

    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.outdir2}
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
        stat_dir=lambda wildcards:os.path.join(wildcards.dir,"stats") if not config['gcp'] else \
                     os.path.join(workflow.default_remote_prefix,wildcards.dir,"stats")
    run:
        print(params.stat_dir)
        shell(f"mkdir -p {params.stat_dir}")
        # pathlib.Path(params.stat_dir).mkdir(exist_ok=True)
        random_index_fasta_path=os.path.join(PACKAGE_DIR,'files','random_index_v1.fa') if barcode_version=='V1' else \
                                os.path.join(PACKAGE_DIR,'files','random_index_v2','random_index_v2.fa')
        index_seq_dict = parse_index_fasta(random_index_fasta_path)
        index_name_dict = {v: k for k, v in index_seq_dict.items()}
        stat_list = []
        for path in input:
            *uid, suffix = os.path.basename(path).split('-')
            lane = suffix.split('.')[0]
            uid = '-'.join(uid)
            single_df=read_cutadapt_result(path)
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