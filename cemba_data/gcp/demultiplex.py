import os,sys
import pandas as pd
import cemba_data
import glob
from snakemake.io import glob_wildcards
PACKAGE_DIR=cemba_data.__path__[0]
from cemba_data.demultiplex.fastq_dataframe import _parse_v2_fastq_path

def make_v2_fastq_df(fq_dir,fq_ext,run_on_gcp):
	# For example: UWA7648_CX05_A10_2_P8-1-O4_22F25JLT3_S15_L001_I1_001.fastq.gz
	# depth = 2
	if run_on_gcp:
		from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
		GS = GSRemoteProvider()
		os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser('~/.config/gcloud/application_default_credentials.json')

		bucket_name = fq_dir.replace('gs://', '').split('/')[0]
		indir = '/'.join(fq_dir.replace('gs://', '').split('/')[1:])
		files = GS.client.list_blobs(bucket_name, prefix=indir)
		input_files = ["gs://"+bucket_name+"/"+file.name for file in files if file.name.endswith(f"{fq_ext}.gz")]
	else:
		input_files=glob.glob(fq_dir) # * should be included in fq_dir, is fastq_pattern
	R=[]
	for file in input_files:
		R.append(_parse_v2_fastq_path(file))
	df = pd.DataFrame(R)
	df.fastq_path=df.fastq_path.apply(lambda x:str(x))
	return df

def get_fastq_info(fq_dir,outdir,fq_ext,run_on_gcp):
	if os.path.exists("fastq_info.txt"):
		df=pd.read_csv("fastq_info.txt",sep='\t')
		# need to write to file, otherwise, snakemake will call this function multiple times.
		return df
	df=make_v2_fastq_df(fq_dir,fq_ext,run_on_gcp)
	# df['fastq_path']=df.apply(lambda row:os.path.join(row.indir,'-'.join(row.loc[['plate','multiplex_group','primer_name']].map(str).tolist())+"_"+"_".join(row.loc[['ID','pns','lane','read_type','suffix']].map(str).tolist())+f".{fq_ext}.gz"),axis=1)
	# df['uid']=df.plate.map(str)+'-'+df.multiplex_group.map(str)+'-'+df.primer_name.map(str) #f'{plate}-{multiplex_group}-{primer_name}')
	df=df.loc[df.read_type.isin(['R1','R2'])]
	assert df.groupby(['lane','read_type'])['uid'].nunique().nunique()==1
	df=df.loc[df.read_type=='R1']
	df.rename(columns={'fastq_path':'R1'},inplace=True)
	df['R2']=df.R1.apply(lambda x:x.replace('_R1_','_R2_'))
	df['stats_out']=df.apply(lambda row: os.path.join(outdir, f"{row.uid}/lanes/{row.uid}-{row.lane}.demultiplex.stats.txt"),axis=1) #"{dir}/{uid}/lanes/{uid}-{lane}.demultiplex.stats.txt"
	df.to_csv("fastq_info.txt",sep='\t',index=False)
	return df

def get_lanes_info(outdir,barcode_version):
	#  uid={plate}-{multiplex_group}-{primer_name}
	if os.path.exists("lane_info.txt"):
		df1=pd.read_csv("lane_info.txt",sep='\t')
		df1.fastq_path = df1.fastq_path.apply(lambda x: eval(x))
		return df1
	uids,plates,multiple_groups,primer_names,lanes,index_names,read_types=glob_wildcards(os.path.join(outdir,"{uid}/lanes/{plate}-{multiplex_group}-{primer_name}-{lane}-{index_name}-{read_type}.fq.gz"))
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
	df['fastq_path']=df.apply(
		lambda row:os.path.join(outdir,row.uid,"lanes",'-'.join(
			row.loc[['uid','lane','index_name','read_type']].map(str).tolist())+".fq.gz"),axis=1)
	if barcode_version == 'V2' and df['multiplex_group'].nunique() == 1:
		df['real_multiplex_group']=df.index_name.apply(
			lambda x:((int(x[1:])-1) % 12) // 2 + 1 if 'unknow' not in x.lower() else 'NA'
		)
	else:
		df['real_multiplex_group']=df.multiplex_group.tolist()
	df=df.loc[df.real_multiplex_group !='NA']
	# new uid (real uid)
	df['uid']= df.plate.map(str)+'-'+df.real_multiplex_group.map(str)+'-'+df.primer_name.map(str) #{plate}-{multiplex_group}-{primer_name}

	# Put multiple lanes fastq into one list
	df1 = df.loc[:, ['uid', 'index_name', 'read_type',
					 'fastq_path']].groupby(
		['uid', 'index_name', 'read_type'],as_index=False).agg(lambda x: x.tolist())
	df1['fastq_out'] = df1.apply(lambda row:os.path.join(outdir, row.uid, "fastq", '-'.join(row.loc[['uid', 'index_name', 'read_type']].map(str).tolist()) + ".fq.gz"), axis=1)
	df1.to_csv("lane_info.txt",sep='\t',index=False)
	return df1

def get_demultiplex_skypilot_yaml():
	skypilot_template = os.path.join(PACKAGE_DIR, "gcp", 'yaml', "demultiplex.yaml")
	with open(skypilot_template) as f:
		template = f.read()
	print(template)

def prepare_demultiplex(fq_dir="fastq",remote_prefix="mapping",outdir="test",
						fq_ext="fastq",barcode_version="V2",env_name='base',
						region='us-west1',keep_remote=False,gcp=True,
						skypilot_template=None,n_jobs=96,job_name="demultiplex",
						workdir="./",output=None):
	workdir = os.path.abspath(os.path.expanduser(workdir))
	CMD=f"yap-gcp run_demultiplex --fq_dir {fq_dir} --remote_prefix {remote_prefix} --outdir {outdir} \
					--fq_ext {fq_ext} --barcode_version {barcode_version} \
					--gcp {gcp} --region {region} --keep_remote {keep_remote} --n_jobs {n_jobs}"
	if not env_name is None:
		CMD=f"conda activate {env_name} \n  "+ CMD
	if skypilot_template is None:
		skypilot_template=os.path.join(PACKAGE_DIR,"gcp",'yaml',"demultiplex.yaml")
	with open(skypilot_template) as f:
		template = f.read()
	if output is None:
		print(template.format(job_name=job_name, workdir=workdir,
								CMD=CMD))
	else:
		with open(os.path.abspath(os.path.expanduser(output)), 'w') as f:
			f.write(template.format(job_name=job_name, workdir=workdir,
								CMD=CMD,env_name=env_name))

	# print(f"To run this job: sky spot launch -y -n {job_name} -y {output} [spot] \n")
	print(f"To run: sky launch -y -n {job_name} {output}")

def run_demultiplex(fq_dir="fastq",remote_prefix="mapping",outdir="test",
						fq_ext="fastq",barcode_version="V2",
						gcp=True,region='us-west1',keep_remote=False,
						n_jobs=96):
	smk1=os.path.join(PACKAGE_DIR,"gcp",'smk',"demultiplex.Snakefile")
	smk2 = os.path.join(PACKAGE_DIR, "gcp", 'smk', "merge_lanes.Snakefile")

	# Demultiplex
	config_str=f'--config gcp={gcp} fq_dir="{fq_dir}" outdir="{outdir}" fq_ext="{fq_ext}" barcode_version="{barcode_version}" '
	common_str=f"--default-remote-prefix {remote_prefix} --default-remote-provider GS --google-lifesciences-region {region} "
	if keep_remote:
		common_str+="--keep-remote "
	CMD1 = f"snakemake -s {smk1} {config_str} {common_str} -j {n_jobs} \n  "

	# Merge lanes
	CMD2 = f"snakemake -s {smk2} {config_str} {common_str} -j {n_jobs} \n  "

	for cmd in [CMD1, CMD2]:
		print(f"CMD: {cmd}")
		os.system(cmd)