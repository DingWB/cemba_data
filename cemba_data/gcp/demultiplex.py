import os,sys
import pandas as pd
import cemba_data
import glob
from snakemake.io import glob_wildcards
PACKAGE_DIR=cemba_data.__path__[0]
from cemba_data.demultiplex.fastq_dataframe import _parse_v2_fastq_path
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

def make_v2_fastq_df(fq_dir,run_on_gcp):
	# For example: UWA7648_CX05_A10_2_P8-1-O4_22F25JLT3_S15_L001_I1_001.fastq.gz
	# depth = 2
	if run_on_gcp:
		GS = GSRemoteProvider()
		os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser('~/.config/gcloud/application_default_credentials.json')

		bucket_name = fq_dir.replace('gs://', '').split('/')[0]
		indir = '/'.join(fq_dir.replace('gs://', '').split('/')[1:])
		files = GS.client.list_blobs(bucket_name, prefix=indir, match_glob='**{.fq.gz,.fastq.gz}')
		input_files = ["gs://"+bucket_name+"/"+file.name for file in files]
	else:
		input_files=glob.glob(fq_dir) # * should be included in fq_dir, is fastq_pattern
	R=[]
	for file in input_files:
		R.append(_parse_v2_fastq_path(file))
	df = pd.DataFrame(R)
	df.fastq_path=df.fastq_path.apply(lambda x:str(x))
	return df

def get_fastq_info(fq_dir,outdir,run_on_gcp):
	if os.path.exists("fastq_info.txt"):
		df=pd.read_csv("fastq_info.txt",sep='\t')
		# need to write to file, otherwise, snakemake will call this function multiple times.
		return df
	df=make_v2_fastq_df(fq_dir,run_on_gcp)
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

def get_fastq_dirs(remote_prefix=None):
	GS = GSRemoteProvider()
	os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
		'~/.config/gcloud/application_default_credentials.json')
	bucket_name = remote_prefix.split('/')[0]
	indir = '/'.join(remote_prefix.replace('gs://', '').split('/')[1:])
	if indir == '':
		prefix=None
	else:
		prefix=indir
	bucket = GS.client.bucket(bucket_name)
	files = bucket.list_blobs(prefix=prefix, match_glob='**{-R1.fq.gz,-R1.fastq.gz}')
	fastq_dirs=[]
	for file in files:
		if 'fastq/' not in file.name:
			continue
		path=file.name.split('/')[1]
		if path not in fastq_dirs:
			fastq_dirs.append(path)
	return fastq_dirs

def get_demultiplex_skypilot_yaml():
	skypilot_template = os.path.join(PACKAGE_DIR, "gcp", 'yaml', "skypilot.yaml")
	with open(skypilot_template) as f:
		template = f.read()
	print(template)

def prepare_demultiplex(fq_dir="fastq",remote_prefix="mapping",outdir="test",
						barcode_version="V2",env_name='base',
						region='us-west1',keep_remote=False,gcp=True,
						skypilot_template=None,n_jobs=96,job_name="demultiplex",
						workdir="./",output=None):
	workdir = os.path.abspath(os.path.expanduser(workdir))
	CMD=f"yap-gcp run_demultiplex --fq_dir {fq_dir} --remote_prefix {remote_prefix} --outdir {outdir} \
					--barcode_version {barcode_version} \
					--gcp {gcp} --region {region} --keep_remote {keep_remote} --n_jobs {n_jobs}"
	if not env_name is None:
		CMD=f"conda activate {env_name} \n  "+ CMD
	if skypilot_template is None:
		skypilot_template=os.path.join(PACKAGE_DIR,"gcp",'yaml',"skypilot.yaml")
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
						barcode_version="V2",
						gcp=True,region='us-west1',keep_remote=False,
						n_jobs=96):
	smk1=os.path.join(PACKAGE_DIR,"gcp",'smk',"demultiplex.Snakefile")
	smk2 = os.path.join(PACKAGE_DIR, "gcp", 'smk', "merge_lanes.Snakefile")

	# Demultiplex
	config_str=f'--config gcp={gcp} fq_dir="{fq_dir}" outdir="{outdir}" barcode_version="{barcode_version}" '
	common_str=f"--default-remote-prefix {remote_prefix} --default-remote-provider GS --google-lifesciences-region {region} "
	if keep_remote:
		common_str+="--keep-remote "
	CMD1 = f"snakemake -s {smk1} {config_str} {common_str} -j {n_jobs} \n  "

	# Merge lanes
	CMD2 = f"snakemake -s {smk2} {config_str} {common_str} -j {n_jobs} \n  "

	for cmd in [CMD1, CMD2]:
		print(f"CMD: {cmd}")
		os.system(cmd)

def make_gcp_snakefile(fastq_prefix,subdir):
	GS = GSRemoteProvider()
	os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
		'~/.config/gcloud/application_default_credentials.json')

	config = get_configuration(os.path.join(fastq_prefix,'mapping_config.ini'))
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	if mode == 'mc':
		config_str = mc_config_str(config)
	elif mode == 'mct':
		config_str = mct_config_str(config)
	elif mode == 'm3c':
		config_str = m3c_config_str(config)
	elif mode == '4m':
		config_str = _4m_config_str(config)
	else:
		raise ValueError(f'Unknown mode {mode}')
	print('Making Snakefile based on mapping config INI file. The parameters are:')
	print(config_str)

	with open(PACKAGE_DIR / f'mapping/Snakefile_template/{mode}.Snakefile') as f:
		snake_template = f.read()

	sub_folder=os.path.join(fastq_prefix,subdir)
	if not os.path.exists(sub_folder):
		os.makedirs(sub_folder,exist_ok=True)
	cell_ids = GS.glob_wildcards(os.path.join(sub_folder,"fastq/{cell_id}-R1.fq.gz"))[0]
	if len(cell_ids) == 0: # length should be 64
		raise ValueError(f"No cell fastq were identified under {sub_folder}/fastq")
	cell_id_str = f'CELL_IDS = {cell_ids}\n'

	total_snakefile = config_str + cell_id_str + snake_template
	with open(os.path.join(sub_folder,'Snakefile'), 'w') as f:
		f.write(total_snakefile)
	return


def prepare_mapping(fastq_prefix="gs://mapping_example/test_gcp",
					config_path="config.ini",aligner='hisat-3n',
					tmp_dir="mapping_gcp",chunk_size=2,
					region='us-west1',keep_remote=False,gcp=True,
					skypilot_template=None,job_name='mapping',
					env_name='base',n_jobs=96):
	outdir=os.path.abspath(os.path.expanduser(tmp_dir))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	fastq_dirs=get_fastq_dirs(fastq_prefix)
	if len(fastq_dirs)==0:
		raise ValueError(f"Please check gs://{remote_prefix}/{outdir} and make sure this is correct, cause no fastq dirs were detected")
	with open(os.path.join(outdir,"fastq_dirs.txt"),'w') as f:
		for d in fastq_dirs:
			f.write(d+'\n')

	# split fastq_dirs into multiple files, running with different skypilot nodes
	j=0
	i=0
	while i +chunk_size < len(fastq_dirs):
		with open(os.path.join(outdir, f"fastq_dirs_{j}"), 'w') as f:
			for d in fastq_dirs[i:i+chunk_size]:
				f.write(d + '\n')
		i+=chunk_size
		j+=1

	if skypilot_template is None:
		skypilot_template=os.path.join(PACKAGE_DIR,"gcp",'yaml',"skypilot.yaml")
	with open(skypilot_template) as f:
		template = f.read()

	CMD = f'yap-gcp run_mapping --fastq_prefix {fastq_prefix} \
						--config_path {config_path} --aligner {aligner} \
						--gcp {gcp} --region {region} \
						--keep_remote {keep_remote} --n_jobs {n_jobs} \
						--node_rank "$SKYPILOT_NODE_RANK"'
	if not env_name is None:
		CMD=f"conda activate {env_name} \n  "+ CMD
	with open(os.path.abspath(os.path.expanduser(output)), 'w') as f:
		f.write(template.format(job_name=job_name, workdir=outdir,
								CMD=CMD, env_name=env_name))

def run_mapping(fastq_prefix="gs://mapping_example/test_gcp",
				gcp=True,region='us-west1',keep_remote=False,
				config_path="config.ini",aligner='hisat-3n',
				n_jobs=96,node_rank=0):
	if not os.path.exists(fastq_prefix):
		os.makedirs(fastq_prefix,exist_ok=True) #on loal GCP VM machine
	os.system(f"cp {config_path} {fastq_prefix}/mapping_config.ini")

	input_fastq_dir=f"fastq_dirs_{node_rank}"
	with open(input_fastq_dir,'r') as f:
		subdirs=f.read().strip().split('\n')

	common_str = f'--config gcp={gcp} local_fastq=False -j {n_jobs} --default-remote-provider GS --google-lifesciences-region {region} '
	if keep_remote:
		config_str+="--keep-remote "
	cmds=[]
	for subdir in subdirs:
		make_gcp_snakefile(fastq_prefix,subdir) #
		# mapping_config.ini need to be under local_output_dir
		cmd_str=f"--default-remote-prefix {fastq_prefix}/{subdir} "
		# there should be fastq dir under default-remote-prefix
		cmd=f"snakemake -s {fastq_prefix}/{subdir}/Snakefile {common_str} {cmd_str}"
		cmds.append(cmd)

	for cmd in cmds:
		print(f"CMD: {cmd}")
		os.system(cmd)
