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
else:
outdir=config["outdir"] if 'outdir' in config else 'mapping'
barcode_version = config["barcode_version"] if 'barcode_version' in config else "V2"

def get_lanes_info(outdir,barcode_version):
    #  uid={plate}-{multiplex_group}-{primer_name}
    outname="lane_info.txt"
    if not os.path.exists(outname):
        uids,plates,multiple_groups,primer_names,lanes,index_names,read_types=\
                glob_wildcards(os.path.join(outdir,\
                "{uid}/lanes/{plate}-{multiplex_group}-{primer_name}-{lane}-{index_name}-{read_type}.fq.gz"))
        if len(uids)==0:
            print("Run demultiplex.smk first, then run merge_lanes.smk !")
            return None
        df=pd.DataFrame.from_dict(
                                    {'uid':uids,
                                    'plate':plates,
                                    'multiplex_group':multiple_groups,
                                    'primer_name':primer_names,
                                    'lane':lanes,
                                    'index_name':index_names,
                                    'read_type':read_types}
                                    ).drop_duplicates()
        df['fastq_path']=df.apply(lambda row:os.path.join(outdir,\
                                row.uid,"lanes",\
                                '-'.join(row.loc[['uid','lane','index_name','read_type']].map(str).tolist())+\
                                ".fq.gz"),axis=1)
        if barcode_version == 'V2' and df['multiplex_group'].nunique() == 1:
            df['real_multiplex_group']=df.index_name.apply(lambda x:((int(x[1:])-1) % 12) // 2 + 1 if 'unknow' not in x.lower() else 'NA')
        else:
            df['real_multiplex_group']=df.multiplex_group.tolist()
        df=df.loc[df.real_multiplex_group !='NA']
        # new uid (real uid)
        df['uid']= df.plate.map(str)+'-'+df.real_multiplex_group.map(str)+'-'+df.primer_name.map(str) #{plate}-{multiplex_group}-{primer_name}
        df.to_csv(outname,sep='\t',index=False)
    df=pd.read_csv(outname,sep='\t')
    return df

df=get_lanes_info(outdir,barcode_version)
if df is None:
    print("Merging is already done.")
    os._exit(1) #sys.exit()
df1=df.loc[:,['uid','index_name','read_type','fastq_path']].groupby(['uid','index_name','read_type'],as_index=False).agg(lambda x:x.tolist())
df1['fastq_out']=df1.apply(lambda row:os.path.join(outdir,row.uid,"fastq",\
                        '-'.join(row.loc[['uid','index_name','read_type']].map(str).tolist())+\
                        ".fq.gz"),axis=1)

rule target_merge_lanes:
    input:
        os.path.join(outdir,"stats/lane_info.tsv")

rule run_merge_lanes: #merge the lanes from the same cell_id and read_type, generating cell fastq
    input: # cell_id = uid-index_name
        fqs=lambda wildcards: [local(p) for p in df1.loc[(df1.uid==wildcards.uid) & \
        (df1.index_name==wildcards.index_name) & \
        (df1.read_type==wildcards.read_type)].fastq_path.iloc[0]] #local disk

    output: # plate-multiplex_group-primer_name-index_name-read_type; AMB_220510_8wk_12D_13B_2_P4-1-I15-A12-R1.fq.gz; index_name is random index?
        fq="{dir}/{uid}/fastq/{uid}-{index_name}-{read_type}.fq.gz" # uid = plate - multiplex_group - pcr_index(primer_name); 64 cells (128 fastq) under each uid.

    #params:
    #    dirname=lambda wildcards: f"{wildcards.dir}/{wildcards.uid}/fastq" if not config['gcp'] else \
    #                                            workflow.default_remote_prefix+f"/{wildcards.dir}/{wildcards.uid}/fastq"

    run:
        outdir=pathlib.Path(os.path.dirname(output.fq)).absolute()
        outdir.mkdir(exist_ok=True, parents=True)
        if len(input.fqs) > 1:
            print("More than 1 input were detected, running merge..")
            shell("gzip -cd {input.fqs} | gzip -5 > {output.fq} && rm -f {input.fqs}")
        else:
            os.rename(input.fqs[0],output.fq)

rule cell_qc:
    input:
        fqs=df1['fastq_out'].tolist(), #output of merge_lanes
        csv=os.path.join(outdir,"stats/demultiplex.stats.csv")
    output:
        csv=os.path.join(outdir , 'stats/UIDTotalCellInputReadPairs.csv')
    run:
        demultiplex_df = pd.read_csv(input.csv,index_col=0)
        total_read_pairs_min = int(config['total_read_pairs_min'])
        total_read_pairs_max = int(config['total_read_pairs_max'])

        too_large = demultiplex_df['CellInputReadPairs'] > total_read_pairs_max
        too_small = demultiplex_df['CellInputReadPairs'] < total_read_pairs_min
        judge = too_small | too_large
        unmapped_cells = demultiplex_df[judge]
        print(f'Skip {too_small.sum()} cells due to too less input read pairs (< {total_read_pairs_min})')
        print(f'Skip {too_large.sum()} cells due to too large input read pairs (> {total_read_pairs_max})')
        real_outdir=outdir if not config['gcp'] else workflow.default_remote_prefix+'/'+outdir

        for cell_id, row in unmapped_cells.iterrows():
            uid = row['UID']
            skipped_dir = pathlib.Path(real_outdir).absolute()  /  uid / 'fastq/skipped'
            skipped_dir.mkdir(exist_ok=True, parents=True)

            # move both R1 R2 to skipped files, it will not be included in Snakefile
            for read_type in ['R1', 'R2']:
                fastq_path = pathlib.Path(real_outdir).absolute() / uid / f'fastq/{cell_id}-{read_type}.fq.gz'
                new_path = skipped_dir / f'{cell_id}-{read_type}.fq.gz'
                # if CellInputReadPairs = 0, the FASTQ file do not actually exist, but it does have a row in metadata.
                if fastq_path.exists():
                    os.rename(str(fastq_path),str(new_path))

        # save UID total input reads, for command order
        uid_order = demultiplex_df[~judge].groupby('UID')['CellInputReadPairs'].sum().sort_values(ascending=False)
        uid_order.to_csv(output.csv,header=False)

        doc="""
            # clean
            lane_dirs=pathlib.Path(outdir).absolute() / "*" / "lanes"
            print(f"Remove dirs: {lane_dirs}")
            shell(f"rm -rf {lane_dirs}")
        """

rule write_lane_info:
    input:
        os.path.join(outdir,'stats/UIDTotalCellInputReadPairs.csv')
    output:
        tsv=os.path.join(outdir,"stats/lane_info.tsv")
    run:
        df1.to_csv(output.tsv,sep='\t',index=False)
        if os.path.exists("lane_info.txt"):
            os.remove("lane_info.txt")