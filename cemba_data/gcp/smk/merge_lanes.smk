import os,sys
import pandas as pd
import pathlib
import re
from glob import glob
from cemba_data.gcp import *

# config
default_config={
    'total_read_pairs_min':1,
    'total_read_pairs_max':6000000
}
for key in default_config:
    if key not in config:
        config[key]=default_config[key]

if 'gcp' in config and config["gcp"]:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
    run_on_gcp=True
else:
    run_on_gcp=False
outdir=config["outdir"] if 'outdir' in config else 'mapping'
local_outdir=outdir if not run_on_gcp else workflow.default_remote_prefix+"/"+outdir
barcode_version = config["barcode_version"] if 'barcode_version' in config else "V2"

df_lane=get_lanes_info(local_outdir,barcode_version)
df_lane['fastq_out'] = df_lane.apply(lambda row:os.path.join(outdir, row.uid, "fastq", '-'.join(row.loc[['uid', 'index_name', 'read_type']].map(str).tolist()) + ".fq.gz"), axis=1)
if df_lane is None:
    print("Merging is already done.")
    os._exit(1) #sys.exit()

rule write_lane_info:
    input:
        os.path.join(outdir,'stats/UIDTotalCellInputReadPairs.csv')
    output:
        tsv=os.path.join(outdir,"stats/lane_info.tsv")
    run:
        if os.path.exists("lane_info.txt"):
            os.rename("lane_info.txt",output.tsv)
        else:
            df_lane.to_csv(output.tsv,sep='\t',index=False)

rule merge_lanes: #merge the lanes from the same cell_id and read_type, generating cell fastq
    input: # cell_id = uid-index_name
        fqs=lambda wildcards: [local(p) for p in df_lane.loc[(df_lane.uid==wildcards.uid) & \
        (df_lane.index_name==wildcards.index_name) & \
        (df_lane.read_type==wildcards.read_type)].fastq_path.iloc[0]] #local disk

    output: # plate-multiplex_group-primer_name-index_name-read_type; AMB_220510_8wk_12D_13B_2_P4-1-I15-A12-R1.fq.gz; index_name is random index?
        fq="{dir}/{uid}/fastq/{uid}-{index_name}-{read_type}.fq.gz" # uid = plate - multiplex_group - pcr_index(primer_name); 64 cells (128 fastq) under each uid.

    run:
        outdir=pathlib.Path(os.path.dirname(output.fq)).absolute()
        outdir.mkdir(exist_ok=True, parents=True)
        if len(input.fqs) > 1:
            # print("More than 1 input were detected, running merge..")
            # shell("gzip -cd {input.fqs} | gzip -5 > {output.fq} && rm -f {input.fqs}")
            shell("cat {input.fqs} > {output.fq} && rm -f {input.fqs}")
        else:
            os.rename(input.fqs[0],output.fq)

rule cell_qc:
    input:
        fqs=df_lane['fastq_out'].tolist(), #output of merge_lanes
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
        real_outdir=outdir if not run_on_gcp else workflow.default_remote_prefix+'/'+outdir

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