import os,sys
import pandas as pd
import cemba_data
import glob
import numpy as np
from cemba_data.mapping.pipelines import make_gcp_snakefile
from snakemake.io import glob_wildcards
PACKAGE_DIR=cemba_data.__path__[0]
from cemba_data.demultiplex.fastq_dataframe import _parse_v2_fastq_path
from cemba_data.demultiplex import _parse_index_fasta
os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
		'~/.config/gcloud/application_default_credentials.json')
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import json
with open(os.environ['GOOGLE_APPLICATION_CREDENTIALS'] ,'r') as f:
	D=json.load(f)
gcp_project=D['quota_project_id']

def make_v2_fastq_df(fq_dir,run_on_gcp):
	# For example: UWA7648_CX05_A10_2_P8-1-O4_22F25JLT3_S15_L001_I1_001.fastq.gz
	# depth = 2
	if run_on_gcp:
		GS = GSRemoteProvider(project=gcp_project)
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

def get_fastq_info(fq_dir,run_on_gcp):
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
	df.to_csv("fastq_info.txt",sep='\t',index=False)
	return df

def index_name2multiplex_group(x):
	if 'unknow' not in x.lower():
		return ((int(x[1:]) - 1) % 12) // 2 + 1
	else:
		return 'NA'

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
	df1.to_csv("lane_info.txt",sep='\t',index=False)
	return df1

def get_random_index(UIDs, barcode_version):
	if os.path.exists("random_index.txt"):
		df_index=pd.read_csv("random_index.txt",sep='\t')
		return df_index
	R=[]
	for uid in UIDs:
		random_index_fa = os.path.join(PACKAGE_DIR, 'files',
														 'random_index_v1.fa') if barcode_version == "V1" \
			else os.path.join(PACKAGE_DIR, 'files', 'random_index_v2',
							  'random_index_v2.multiplex_group_' + uid.split('-')[
								  -2] + '.fa') if barcode_version == "V2" \
			else os.path.join(PACKAGE_DIR, 'files', 'random_index_v2', 'random_index_v2.fa')
		index_seq_dict = _parse_index_fasta(random_index_fa)
		index_names=list(index_seq_dict.keys()) #384 random index names (A-P, 1-24)
		for index_name in index_names:
			for read_type in ['R1','R2']:
				R.append([uid,read_type,index_name])
	df_index=pd.DataFrame(R,columns=['uid','read_type','index_name'])
	df_index.rename(columns={'uid':'old_uid'},inplace=True)
	df_index['real_multiplex_group']=df_index.index_name.apply(index_name2multiplex_group)
	df_index['uid']=df_index.loc[:,['old_uid','real_multiplex_group']].apply(
		lambda x:'-'.join([x.old_uid.split('-')[0],str(x.real_multiplex_group),x.old_uid.split('-')[-1]]),axis=1
	)
	df_index.to_csv("random_index.txt", sep='\t', index=False)
	return df_index

def get_fastq_dirs(remote_prefix=None):
	GS = GSRemoteProvider(project=gcp_project)
	bucket_name = remote_prefix.replace('gs://', '').split('/')[0]
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
		path=file.name.split('/')[1] #uid
		if bucket.blob(f"{prefix}/{path}/MappingSummary.csv.gz").exists():
			continue # existed, skip
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
						tmp_dir="demultiplex_gcp_tmp",disk_size=3072,
						region='us-west1',keep_remote=False,gcp=True,
						skypilot_template=None,n_jobs=16,job_name="demultiplex",
						image="bican",output='run_demultiplex.yaml'):
	"""
		Prepare the skypilot yaml file to run demultiplex on GCP.

	Parameters
	----------
	fq_dir : path
		address of fastq files on GCP, for example: gs://mapping_example/fastq/test_fastq
	remote_prefix : path
		remote_prefix passed to snakemake --default-remote-prefix, used to write
		the output of demultiplex and mapping to this remote prefix/outdir, for
		example, --remote_prefix=="mapping_example", remote_prefix could be
		different from fq_dir.
	outdir : str
		Used together with remote_prefix, the resulting demultiplex and mapping
		will be upload to remote_prefix/outdir on GCP.
	barcode_version : str
		default is V2
	env_name : str
		conda env used to run the snakemake pipeline, default is base
	region : str
		default is us-west1
	keep_remote : bool
		Whether to keep remote fastq files, default is False.
	gcp : bool
		Whether to run on GCP, default is True, to run on local, please normal
		YAP pipelien.
	skypilot_template : path
		skypilot template used to generate output, default is None, will used
		PACKAGE_DIR/gcp/yaml/skypilot.yaml, see https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/gcp/yaml/skypilot.yaml,
		users could also provide their own template, but keep wildcard
		{job_name}, {workdir},{CMD},and {env_name} in this template.
	n_jobs : int
	job_name : str
	workdir : str
		default workdir is current dir.
	output : path
		output of the generated skypilot yaml file according to the template.

	Returns
	-------

	"""
	workdir = os.path.abspath(os.path.expanduser(tmp_dir))
	if not os.path.exists(workdir):
		os.makedirs(workdir)
	CMD=f"yap-gcp run_demultiplex --fq_dir {fq_dir} --remote_prefix {remote_prefix} --outdir {outdir} \
--barcode_version {barcode_version} \
--gcp {gcp} --region {region} --keep_remote {keep_remote} --n_jobs {n_jobs}"
	if not env_name is None:
		CMD=f"conda activate {env_name} \n  "+ CMD
	if skypilot_template is None:
		skypilot_template=os.path.join(PACKAGE_DIR,"gcp",'yaml',"skypilot.yaml")
	else:
		skypilot_template=os.path.expanduser(skypilot_template)
	with open(skypilot_template) as f:
		template = f.read()
	if output is None:
		print(template.format(job_name=job_name, workdir=workdir,
							  CMD=CMD,env_name=env_name,
							  n_node=1,image=image,disk_size=disk_size))
	else:
		with open(os.path.abspath(os.path.expanduser(output)), 'w') as f:
			f.write(template.format(job_name=job_name, workdir=workdir,
									CMD=CMD,env_name=env_name,
									n_node=1,image=image,disk_size=disk_size))

	# print(f"To run this job: sky spot launch -y -n {job_name} {output} [spot] \n")
	print(f"To run: sky launch -y -i 10 -n {job_name} {output}")

def run_demultiplex(fq_dir="fastq",remote_prefix="mapping",outdir="test",
						barcode_version="V2",
						gcp=True,region='us-west1',keep_remote=False,
						n_jobs=16,print_only=False):
	"""
		This function need to be executed on the GCP VM machine. Please see
		prepare_demultiplex for parameters.
	Parameters
	----------
	fq_dir :
	remote_prefix :
	outdir :
	barcode_version :
	gcp :
	region :
	keep_remote :
	n_jobs :

	Returns
	-------

	"""
	smk1=os.path.join(PACKAGE_DIR,"gcp",'smk',"demultiplex.smk")

	# Demultiplex
	config_str=f'--scheduler greedy --printshellcmds --rerun-incomplete --config gcp={gcp} fq_dir="{fq_dir}" outdir="{outdir}" barcode_version="{barcode_version}" '
	common_str=f"--default-remote-prefix {remote_prefix} --default-remote-provider GS --google-lifesciences-region {region} "
	if keep_remote:
		common_str+="--keep-remote "
	cmd = f"snakemake -s {smk1} {config_str} {common_str} -j {n_jobs} \n  "

	print(f"CMD: {cmd}")
	if not print_only:
		os.system(cmd)

def prepare_mapping(fastq_prefix="gs://mapping_example/test_gcp",
					config_path="config.ini",aligner='hisat-3n',
					tmp_dir="mapping_gcp_tmp",disk_size=500,
					chunk_size=None,n_node=12,image="bican",
					region='us-west1',keep_remote=False,gcp=True,
					skypilot_template=None,job_name='mapping',
					env_name='base',n_jobs=64,output="run_mapping.yaml"):
	"""
		Prepare the skypilot yaml file to run demultiplex on GCP.

	Parameters
	----------
	fastq_prefix : path
		The fastq files prefix on GCP bucket, for example, after running prepare_demultiplex
		and submited onto the GCP with parameter: --fq_dir gs://mapping_example/fastq/test_fastq \
              --remote_prefix mapping_example --outdir test_gcp_hisat3n;
        Then, the fastq_prefix should be gs://mapping_example/test_gcp_hisat3n,
        under this fastq_prefix, there are many subfolder (each subfolder is a uid),
        under each uid, there are 64 cell level fastq files.
	config_path :path
		config.ini, generated by yap default-mapping-config, for example:
		yap default-mapping-config --mode m3c --barcode_version V2 \
			--bismark_ref "~/Ref/hg38/hg38_ucsc_with_chrL.bismark1" \
			--genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" \
		-	-chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes"  > config.ini
	aligner : str
		default is hisat-3n
	tmp_dir : path
		default is mapping_gcp_tmp, the fastq dir names will be writen into this
		folder and path of this folder would be used as workdir (synced onto
		GCP VM machine sky_workdir)
	chunk_size : int
		default is 2, run two uid (64 cells for each uid) on each GCP node.
	region : str
		GCP VM machine region, default is us-west1
	keep_remote : bool
		default is False
	gcp : bool
		[True]
	skypilot_template : path
		skypilot template used to generate output, default is None, will used
		PACKAGE_DIR/gcp/yaml/skypilot.yaml, see https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/gcp/yaml/skypilot.yaml,
		users could also provide their own template, but keep wildcard
		{job_name}, {workdir},{CMD},and {env_name} in this template. Please Note:
		the instance_type in template should match n_jobs, for example, if n_jobs=96,
		type should be n2-standard-96, other options are n2-standard-64, n2-standard-8 and so on.
	job_name : str
	env_name : str
		conda env name to run snakemake pipeline
	n_jobs : int
		No of CPUs
	output : path
		output of the generated skypilot yaml file according to the template.

	Returns
	-------

	"""
	outdir=os.path.abspath(os.path.expanduser(tmp_dir))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	os.system(f"cp {config_path} {outdir}/mapping_config.ini")
	fastq_dirs=get_fastq_dirs(fastq_prefix)
	if len(fastq_dirs)==0:
		raise ValueError(f"Please check {fastq_prefix} and make sure this is correct, cause no fastq dirs were detected")
	with open(os.path.join(outdir,"fastq_dirs.txt"),'w') as f:
		for d in fastq_dirs:
			f.write(d+'\n')

	# split fastq_dirs into multiple files, running with different skypilot nodes
	if not n_node is None:
		chunk_size=int(np.ceil(len(fastq_dirs)/n_node))
		print(f"n_node:{n_node}; chunk_size: {chunk_size}")
	else:
		print(f"chunk_size: {chunk_size}")
	j=0
	i=0
	while i < len(fastq_dirs):
		with open(os.path.join(outdir, f"fastq_dirs_{j}"), 'w') as f:
			for d in fastq_dirs[i:i+chunk_size]:
				f.write(d + '\n')
		i+=chunk_size
		j+=1
	n_node=j

	if skypilot_template is None:
		skypilot_template=os.path.join(PACKAGE_DIR,"gcp",'yaml',"skypilot.yaml")
	else:
		skypilot_template=os.path.expanduser(skypilot_template)
	with open(skypilot_template) as f:
		template = f.read()

	CMD = f'yap-gcp run_mapping --fastq_prefix {fastq_prefix} \
--config_path "mapping_config.ini" --aligner {aligner} \
--gcp {gcp} --region {region} \
--keep_remote {keep_remote} --n_jobs {n_jobs} \
--node_rank "$SKYPILOT_NODE_RANK"'
	if not env_name is None:
		CMD=f"conda activate {env_name} \n  "+ CMD
	with open(os.path.abspath(os.path.expanduser(output)), 'w') as f:
		f.write(template.format(job_name=job_name, workdir=outdir,
								CMD=CMD, env_name=env_name,
								n_node=n_node,image=image,disk_size=disk_size))
	print(f"To run this job: \nsky spot launch -y {output} [spot] \n")
	print(f"Or: \nsky launch -y -n {job_name} {output}")

def run_mapping(fastq_prefix="gs://mapping_example/test_gcp",
				gcp=True,region='us-west1',keep_remote=False,
				config_path="mapping_config.ini",aligner='hisat-3n',
				n_jobs=64,node_rank=0,print_only=False):
	output_dir=fastq_prefix.replace("gs://","")
	if not os.path.exists(output_dir):
		os.makedirs(output_dir,exist_ok=True) #on loal GCP VM machine
	os.system(f"cp {config_path} {output_dir}/mapping_config.ini")

	if node_rank < 0:
		input_fastq_dir="fastq_dirs.txt"
	elif os.path.exists(f"fastq_dirs_{node_rank}"):
		input_fastq_dir=f"fastq_dirs_{node_rank}"
	else:
		input_fastq_dir = "fastq_dirs.txt"
	info=f"Node Rank: {node_rank}; input: {input_fastq_dir}'"
	print(info)
	log_path=os.path.join(output_dir,f"logs_{node_rank}.txt")
	with open(log_path,'w') as f:
		f.write(info+'\n')
	with open(input_fastq_dir,'r') as f:
		subdirs=f.read().strip().split('\n')

	common_str = f'--default-resources mem_mb=100 --resources mem_mb=50000 --printshellcmds --scheduler greedy --rerun-incomplete --config gcp={gcp} local_fastq=False -j {n_jobs} --default-remote-provider GS --google-lifesciences-region {region} '
	if keep_remote:
		common_str+="--keep-remote "
	cmds=[]
	for subdir in subdirs:
		make_gcp_snakefile(output_dir,subdir,aligner=aligner) #
		# mapping_config.ini need to be under local_output_dir
		cmd_str=f"--default-remote-prefix {output_dir}/{subdir}"
		# there should be fastq dir under default-remote-prefix
		cmd=f"snakemake -s {output_dir}/{subdir}/Snakefile {common_str} {cmd_str}"
		# workdir should be current, not {output_dir}/{subdir}
		cmds.append(cmd)

	for cmd in cmds:
		print(f"CMD: {cmd}")
		if print_only:
			continue
		with open(log_path, 'a') as f:
			f.write(cmd + '\n')
		os.system(cmd)

def yap_pipeline(
	fq_dir="gs://mapping_example/fastq/salk10_test",
	remote_prefix='bican',outdir='salk010_test',
	barcode_version="V2", env_name='base',
	region='us-west1', keep_remote=False, gcp=True,
	n_jobs1=16,n_jobs2=64,
	image="bican",demultiplex_template=None,
	mapping_template=None, genome="~/Ref/hg38_Broad/hg38.fa",
	hisat3n_dna_ref="~/Ref/hg38_Broad/hg38",
	mode='m3c',bismark_ref='~/Ref/hg38/hg38_ucsc_with_chrL.bismark1',
	chrom_size_path='~/Ref/hg38_Broad/hg38.chrom.sizes',
	aligner='hisat-3n',n_node=12,sky_env='sky',disk_size1=3072,
	disk_size2=500):
	if not demultiplex_template is None:
		demultiplex_template=os.path.expanduser(demultiplex_template)
	if not mapping_template is None:
		mapping_template=os.path.expanduser(mapping_template)

	cmd=f'conda activate {env_name} && yap-gcp prepare_demultiplex --fq_dir {fq_dir} --remote_prefix {remote_prefix} \
--outdir {outdir} --barcode_version {barcode_version} --env_name {env_name} \
--region {region} --keep_remote {keep_remote} --gcp {gcp} \
--skypilot_template {demultiplex_template} --n_jobs {n_jobs1} \
--job_name demultiplex --image {image} --disk_size {disk_size1} \
--output run_demultiplex.yaml'
	print(cmd)
	print(f"conda activate {sky_env} && sky launch -y -i 5 -n demultiplex run_demultiplex.yaml")
	fastq_prefix=f"gs://{remote_prefix}/{outdir}"
	os.system(f'yap default-mapping-config --mode {mode} \
--barcode_version {barcode_version} \
--bismark_ref "{bismark_ref}" --genome "{genome}" \
--chrom_size_path "{chrom_size_path}" \
--hisat3n_dna_ref  "{hisat3n_dna_ref}" > config.ini')
	cmd=f'conda activate {env_name} && yap-gcp prepare_mapping --fastq_prefix {fastq_prefix} \
--config_path config.ini --aligner {aligner} \
--tmp_dir mapping_gcp_tmp --n_node {n_node} --image {image} \
--region {region} --keep_remote {keep_remote} --gcp {gcp} \
--skypilot_template {mapping_template} --job_name mapping \
--env_name {env_name} --n_jobs {n_jobs2} --disk_size {disk_size2} \
--output run_mapping.yaml'
	print(cmd)
	print(f"conda activate {sky_env} && sky spot launch -y run_mapping.yaml")

def check_demultiplex(fastq_prefix="gs://mapping_example/novaseq_mapping"):
	GS = GSRemoteProvider(project=gcp_project)
	os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
		'~/.config/gcloud/application_default_credentials.json')
	bucket_name = fastq_prefix.replace('gs://', '').split('/')[0]
	indir = '/'.join(fastq_prefix.replace('gs://', '').split('/')[1:])
	files = GS.client.list_blobs(bucket_name, prefix=indir, match_glob='**{.fq.gz,.fastq.gz}')
	R=[]
	for file in files:
		values=file.name.split('/')
		R.append([values[-3],values[-1]])
	df=pd.DataFrame(R,columns=['uid','fq'])
	print(f"uids: {df.uid.nunique()}")
	print(df.groupby('uid').fq.count())
	a=df.groupby('uid').fq.count()
	print(a.min())
	return df.groupby('uid').fq.count().min()

def cell_qc(fastq_prefix="gs://bican/salk010",
			total_read_pairs_max=6000000,total_read_pairs_min=1,
			sky_env='sky'):
	gsutil = os.path.join(os.path.abspath(os.path.join(os.path.dirname(sys.executable),'../../')),
						  f"{sky_env}/bin/gsutil")
	GS = GSRemoteProvider(project=gcp_project)
	os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
		'~/.config/gcloud/application_default_credentials.json')
	if not os.path.exists("demultiplex.stats.csv"):
		os.system(f"{gsutil} cp -n {fastq_prefix}/stats/demultiplex.stats.csv ./")
	df_cell=pd.read_csv("demultiplex.stats.csv",index_col=0)
	df_cell['col384'] = df_cell.IndexName.apply(lambda x:int(x[1:]) - 1)
	df_cell['row384'] = df_cell.IndexName.apply(lambda x: ord(x[0]) - 65)
	too_large = df_cell['CellInputReadPairs'] > total_read_pairs_max
	too_small = df_cell['CellInputReadPairs'] < total_read_pairs_min
	judge = too_small | too_large
	unmapped_cells = df_cell[judge]
	unmapped_cells['Position']=unmapped_cells.loc[:,['row384','col384']].apply(
		lambda x:'row: '+str(x[0])+'; col: '+str(x[1]),axis=1
	)
	print("Total No. of cells to be removed: ",unmapped_cells.shape[0])
	print(unmapped_cells.Position.value_counts())
	f=open("cell_qc.sh",'w')
	for cell_id,uid in unmapped_cells.UID.items():
		for read_type in ['R1','R2']:
			cmd=f"{gsutil} rm -f {fastq_prefix}/{uid}/fastq/{cell_id}-{read_type}.fq.gz"
			os.system(f"{gsutil} ls -lh {fastq_prefix}/{uid}/fastq/{cell_id}-{read_type}.fq.gz")
			f.write(cmd+'\n')
	f.close()

