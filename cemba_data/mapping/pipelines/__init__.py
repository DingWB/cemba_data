import glob
import pathlib
import subprocess
import os
import pandas as pd
import cemba_data
from .m3c import m3c_config_str
from .mc import mc_config_str
from .mct import mct_config_str
from ._4m import _4m_config_str
from ...utilities import get_configuration
from .mhap import bam2mhap
# from cemba_data.utilities import get_configuration
# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])
def prepare_uid_snakefile(uid_dir, config_str, snake_template):
	cell_ids = [path.name.split('.')[0][:-3] for path in (uid_dir / 'fastq').glob('*R1.fq.gz')]
	cell_id_str = f'CELL_IDS = {cell_ids}\n'

	# no file in this UID, do not make snakefile
	if len(cell_ids) == 0:
		print(f'There is no cell_id parsed from FASTQ files, '
			  f'check the {uid_dir}, make sure things are intact.')
		return

	total_snakefile = config_str + cell_id_str + snake_template
	with open(uid_dir / 'Snakefile', 'w') as f:
		f.write(total_snakefile)
	return


def validate_mapping_config(output_dir):
	output_dir = pathlib.Path(output_dir).absolute()
	config = get_configuration(output_dir / 'mapping_config.ini')
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	if mode.split('-')[0] == 'mc':
		config_str = mc_config_str(config)
	elif mode.split('-')[0] == 'mct':
		config_str = mct_config_str(config)
	elif mode.split('-')[0] == 'm3c':
		config_str = m3c_config_str(config)
	elif mode.split('-')[0] == '4m':
		config_str = _4m_config_str(config)
	else:
		raise ValueError(f'Unknown mode {mode}')

	print(f'Mapping config file looks good. Here is what will be used in generating Snakefile:\n{config_str}')
	return

def make_snakefile(output_dir,aligner="bismark"):
	output_dir = pathlib.Path(output_dir).absolute()
	mapping_config_name = list(output_dir.glob('mapping_config.*'))[0].name
	config = get_configuration(output_dir / mapping_config_name)
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	if mode.split('-')[0] == 'mc':
		config_str = mc_config_str(config)
	elif mode.split('-')[0] == 'mct':
		config_str = mct_config_str(config)
	elif mode.split('-')[0] == 'm3c':
		config_str = m3c_config_str(config)
	elif mode.split('-')[0] == '4m':
		config_str = _4m_config_str(config)
	else:
		raise ValueError(f'Unknown mode {mode}')
	# print('Making Snakefile based on mapping config INI file. The parameters are:')
	# print(config_str)

	if aligner.lower()=="bismark":
		snakefile_path=os.path.join(PACKAGE_DIR, f'files/smk/bismark/{mode.lower()}.smk')
	elif aligner.lower() in ['hisat3n', 'hisat-3n', 'hisat_3n', 'hisat']:
		snakefile_path = os.path.join(PACKAGE_DIR, f'files/smk/hisat3n/{mode.lower()}.smk')
	else:
		raise ValueError(f"Unknown aligner: {aligner}")
	with open(snakefile_path) as f:
		snake_template = f.read()

	for sub_dir in output_dir.iterdir():
		if sub_dir.is_dir():
			if sub_dir.name not in ['stats', 'snakemake']:
				prepare_uid_snakefile(uid_dir=sub_dir,
									  config_str=config_str,
									  snake_template=snake_template)
	return

def make_all_snakefile(output_dir, subdir=None, aligner="hisat-3n", gcp=True,
					   snakemake_template=None, pattern="fastq/{cell_id}-R1.fq.gz"):
	"""

	Parameters
	----------
	output_dir :
	subdir :
	aligner :
	gcp :
	snakemake_template :
	pattern : str
		used to get cell_ids
		fastq/{cell_id}-R1.fq.gz, or "CELL_IDS.tsv"

	Returns
	-------

	"""
	if gcp:
		from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
		import json
		os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser(
			'~/.config/gcloud/application_default_credentials.json')
		with open(os.environ['GOOGLE_APPLICATION_CREDENTIALS'], 'r') as f:
			D = json.load(f)
		gcp_project = D['quota_project_id']
		GS = GSRemoteProvider(project=gcp_project)
	else:
		from snakemake.io import glob_wildcards
	# assert os.path.exists(os.path.join(output_dir,'mapping_config.ini'))
	try:
		mapping_config_name = [file for file in os.listdir(output_dir) if file.startswith('mapping_config.')][0]
	except:
		raise ValueError(f"Could not find mapping_config.* under {output_dir}")
	config = get_configuration(os.path.join(output_dir,mapping_config_name))
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	if mode.split('-')[0] == 'mc':
		config_str = mc_config_str(config)
	elif mode.split('-')[0] == 'mct':
		config_str = mct_config_str(config)
	elif mode.split('-')[0] == 'm3c':
		config_str = m3c_config_str(config)
	elif mode.split('-')[0] =='4m':
		config_str = _4m_config_str(config)
	else:
		print(mode)
		raise ValueError(f'Unknown mode {mode}')
	# print('Making Snakefile based on mapping config INI file. The parameters are:')
	# print(config_str)

	if not snakemake_template is None:
		snakefile_path=os.path.expanduser(snakemake_template)
	else:
		if aligner.lower() == "bismark":
			snakefile_path = os.path.join(PACKAGE_DIR, f'files/smk/bismark/{mode.lower()}.smk')
		elif aligner.lower() in ['hisat3n', 'hisat-3n', 'hisat_3n', 'hisat']:
			snakefile_path = os.path.join(PACKAGE_DIR, f'files/smk/hisat3n/{mode.lower()}.smk')
		else:
			raise ValueError(f"Unknown aligner: {aligner}")

	with open(snakefile_path) as f:
		snake_template = f.read()

	if not subdir is None:
		sub_folder=os.path.join(output_dir,subdir)
		if not os.path.exists(sub_folder):
			os.makedirs(sub_folder,exist_ok=True)
	else:
		sub_folder=output_dir
	if pattern=='CELL_IDS':
		cell_ids=pd.read_csv(os.path.join(sub_folder,pattern),sep='\t',index_col=0).index.tolist()
	elif gcp:
		cell_ids = GS.glob_wildcards(os.path.join(sub_folder,pattern))[0]
		#sub_folder can startwith gs://, if gs:// not present at the beginning, it is also OK
	else:
		cell_ids = glob_wildcards(os.path.join(sub_folder, pattern))[0]

	if len(cell_ids) == 0: # length should be 64
		raise ValueError(f"No cell fastq were identified under {sub_folder}/fastq")
	cell_id_str = f'CELL_IDS = {cell_ids}\n'

	if aligner=="bismark":
		total_snakefile = config_str + cell_id_str + snake_template
	else: # hisat-3n
		total_snakefile = cell_id_str + snake_template
		if not os.path.exists(os.path.join(output_dir,'snakemake')):
			os.makedirs(os.path.join(output_dir,'snakemake'),exist_ok=True)
		subprocess.run(['touch', os.path.join(output_dir,'snakemake/hisat3n')], check=True)
	with open(os.path.join(sub_folder,'Snakefile'), 'w') as f:
	# with open(f'{subdir}.smk', 'w') as f:
		f.write(total_snakefile)
	return

def make_snakefile_hisat3n(output_dir,aligner='hisat-3n'):
	output_dir = pathlib.Path(output_dir)

	mapping_config_name = list(output_dir.glob('mapping_config.*'))[0].name

	config = get_configuration(output_dir / mapping_config_name)
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')

	skip_dirs = ['stats', 'snakemake', 'scool']
	mapping_job_dirs = [p for p in output_dir.glob('*')
						if p.is_dir() and (p.name not in skip_dirs)]

	snakemake_dir = output_dir / 'snakemake'
	snakemake_dir.mkdir(exist_ok=True)
	stats_dir = output_dir / 'stats'
	stats_dir.mkdir(exist_ok=True)

	package_dir = cemba_data.__path__[0]
	if aligner.lower()=="bismark":
		snakefile_path=os.path.join(PACKAGE_DIR, f'files/smk/bismark/{mode.lower()}.smk')
	elif aligner.lower() in ['hisat3n', 'hisat-3n', 'hisat_3n', 'hisat']:
		snakefile_path = os.path.join(PACKAGE_DIR, f'files/smk/hisat3n/{mode.lower()}.smk')
	else:
		raise ValueError(f"Unknown aligner: {aligner}")
	if not pathlib.Path(snakefile_path).exists():
		print('Possible snakefile templates:')
		for p in pathlib.Path(f'{package_dir}/hisat3n/snakefile/').glob('*.smk'):
			print(p)
		raise ValueError(f'Mode {mode} not supported, '
						 f'because Snakefile {snakefile_path} not found.')

	for p in mapping_job_dirs:
		subprocess.run(['cp', f'{output_dir}/{mapping_config_name}',
						f'{p}/{mapping_config_name}'], check=True)
		subprocess.run(['cp', snakefile_path, f'{p}/Snakefile'], check=True)

	# leave a flag to indicate using hisat-3n pipeline
	subprocess.run(['touch', f'{output_dir}/snakemake/hisat3n'], check=True)
	return

def write_qsub_commands(output_dir, cores_per_job, total_memory_gb=None,
						script_dir=None,fastq_server='local'):
	if total_memory_gb is None:
		total_memory_gb=2*cores_per_job
	if fastq_server!='local':
		config_par=f"--config fastq_server='{fastq_server}' "
	else:
		config_par=''
	cmds = {}
	snake_files = list(output_dir.glob('*/Snakefile'))
	for snake_file in snake_files:
		uid = snake_file.parent.name
		cmd = f"""snakemake -d {snake_file.parent} --snakefile {snake_file} {config_par} -j {cores_per_job} --rerun-incomplete --scheduler greedy --default-resources mem_mb=100 \
--resources mem_mb={int(1024 * total_memory_gb)} && rm -rf {snake_file.parent}/.snakemake"""
		cmds[uid] = cmd #--resources mem_mb is the limitation.
	script_path = script_dir / 'snakemake_cmd.txt'
	with open(script_path, 'w') as f:
		try:
			uid_order = pd.read_csv(
				output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0,header=None
			).squeeze().sort_values(ascending=False)
			for uid in uid_order.index:
				if uid in cmds:
					f.write(cmds.pop(uid) + '\n')
			try:
				assert len(cmds) == 0
			except AssertionError as e:
				print(cmds)
				print(uid_order)
				raise e
		except FileNotFoundError:
			# uid_order file do not exist (when starting from cell FASTQs)
			for cmd in cmds.values():
				f.write(cmd + '\n')
	return script_path

def write_gcp_skypolit_yaml(output_dir, template_path):
	output_dir=pathlib.Path(output_dir).absolute()
	config = get_configuration(output_dir / 'mapping_config.ini')
	try:
		mode = config['mode']
	except KeyError:
		raise KeyError('mode not found in the config file.')
	if template_path is None:
		print("Using template: "+str(PACKAGE_DIR)+ f'/files/skypilot_template.yaml')
		with open(PACKAGE_DIR / f'files/skypilot_template.yaml') as f:
			template = f.read()
	else:
		with open(template_path) as f:
			template = f.read()
	sky_dir=output_dir/"snakemake/gcp"
	sky_dir.mkdir(exist_ok=True, parents=True)
	snake_files = list(output_dir.glob('*/Snakefile'))
	f_cmd=open(sky_dir / "sky_spot.sh",'w')
	for snake_file in snake_files:
		uid = snake_file.parent.name
		yaml_path = sky_dir / f"{uid}.yaml"
		outdir=output_dir.name
		workdir=str(output_dir)+f"/{uid}"
		print(yaml_path)
		name=uid.lower().replace('_','-')
		with open(yaml_path,'w') as f:
			f.write(template.format(name=name,uid=uid,workdir=workdir,outdir=outdir))
		f_cmd.write(f"sky spot launch -n {name} -y "+str(yaml_path)+"\n")
	f_cmd.close()

def write_sbatch_commands(output_dir, cores_per_job, script_dir, total_mem_mb, queue):
	output_dir_name = output_dir.name
	outdir=str(output_dir.absolute())
	cmds = {}
	snake_files = list(output_dir.glob('*/Snakefile'))
	for snake_file in snake_files:
		uid = snake_file.parent.name
		cmd = f'snakemake ' \
			  f'-d {outdir}/{snake_file.parent.name} ' \
			  f'--snakefile {outdir}/{snake_file.parent.name}/Snakefile ' \
			  f'-j {cores_per_job} ' \
			  f'--default-resources mem_mb=100 --rerun-incomplete ' \
			  f'--resources mem_mb={total_mem_mb} ' \
			  f'--rerun-incomplete ' \
			  f'&& test -f "{outdir}/{snake_file.parent.name}/MappingSummary.csv.gz && && rm -rf {outdir}/{snake_file.parent.name}/.snakemake"'
		cmds[uid] = cmd
	script_path = script_dir / f'snakemake_{queue}_cmd.txt'
	with open(script_path, 'w') as f:
		try:
			uid_order = pd.read_csv(
				output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0,header=None
			).squeeze().sort_values(ascending=False)
			for uid in uid_order.index:
				if uid in cmds:
					f.write(cmds.pop(uid) + '\n')
			try:
				assert len(cmds) == 0
			except AssertionError as e:
				print(cmds)
				print(uid_order)
				raise e
		except FileNotFoundError:
			# uid_order file do not exist (when starting from cell FASTQs)
			for cmd in cmds.values():
				f.write(cmd + '\n')
	return f'{outdir}/snakemake/sbatch/snakemake_{queue}_cmd.txt'

def prepare_qsub(name, snakemake_dir, total_jobs, cores_per_job, total_memory_gb,fastq_server):
	memory_gb_per_core = int(total_memory_gb / cores_per_job) if not total_memory_gb is None else 2
	output_dir = snakemake_dir.parent
	qsub_dir = snakemake_dir / 'qsub'
	qsub_dir.mkdir(exist_ok=True)
	script_path = write_qsub_commands(output_dir, cores_per_job, total_memory_gb,
									  script_dir=qsub_dir,fastq_server=fastq_server)
	qsub_str = f"""
#!/bin/bash
#$ -N yap{name}
#$ -V
#$ -l h_rt=99:99:99
#$ -l s_rt=99:99:99
#$ -wd {qsub_dir}
#$ -e {qsub_dir}/qsub.error.log
#$ -o {qsub_dir}/qsub.output.log
#$ -pe smp 1
#$ -l h_vmem=3G

yap qsub \
--command_file_path {script_path} \
--working_dir {qsub_dir} \
--project_name y{name} \
--total_cpu {int(cores_per_job * total_jobs)} \
--qsub_global_parms "-pe smp={cores_per_job};-l h_vmem={memory_gb_per_core}G"
"""
	qsub_total_path = qsub_dir / 'qsub.sh'
	with open(qsub_total_path, 'w') as f:
		f.write(qsub_str)
	print('#' * 60)
	print(f"IF YOU USE QSUB ON GALE: ")
	print(f"All snakemake commands need to be executed "
		  f"were included in {qsub_total_path}")
	print(f"You just need to qsub this script to "
		  f"map the whole library in {output_dir}")
	print(f"You can also change the per job parameters in {script_path} "
		  f"or change the global parameters in {qsub_total_path}")
	print(f"Read 'yap qsub -h' if you want to have more options about sbatch. "
		  f"Alternatively, you can sbatch the commands in {script_path} by yourself, "
		  f"as long as they all get successfully executed.")
	print('#' * 60 + '\n')
	return

def prepare_sbatch(name, snakemake_dir, queue,total_memory_gb=None):
	input_total_mem_mb = total_memory_gb * 1024 if not total_memory_gb is None else None
	output_dir = snakemake_dir.parent
	output_dir_name = output_dir.name
	outdir=str(output_dir.absolute())
	mode = get_configuration(output_dir / 'mapping_config.ini')['mode']

	if queue == 'skx-normal':
		sbatch_cores_per_job = 96
		if mode.split('-')[0] == 'm3c':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 204800
		elif mode.split('-')[0] == '4m':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 204800
		elif mode.split('-')[0] == 'mc':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 204800
		elif mode.split('-')[0] == 'mct':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 204800
		else:
			raise KeyError(f'Unknown mode {mode}')
	elif queue == 'normal':
		sbatch_cores_per_job = 64
		if mode.split('-')[0] == 'm3c':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == '4m':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == 'mc':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == 'mct':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		else:
			raise KeyError(f'Unknown mode {mode}')
	else: # queue == 'shared':
		sbatch_cores_per_job = 64
		if mode.split('-')[0] == 'm3c':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == '4m':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == 'mc':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		elif mode.split('-')[0] == 'mct':
			time_str = "48:00:00"
			total_mem_mb = input_total_mem_mb if not input_total_mem_mb is None else 64*2*1024
		else:
			raise KeyError(f'Unknown mode {mode}')
	# else:
	#     raise ValueError(f'Unknown queue {queue}')
	sbatch_dir = snakemake_dir / 'sbatch'
	sbatch_dir.mkdir(exist_ok=True)

	script_path = write_sbatch_commands(output_dir,
										cores_per_job=sbatch_cores_per_job,
										script_dir=sbatch_dir,
										total_mem_mb=total_mem_mb,
										queue=queue)
	# the path here is using stampede path
	sbatch_cmd = f'yap sbatch ' \
				 f'--project_name {name} ' \
				 f'--command_file_path {script_path} ' \
				 f'--working_dir {outdir}/snakemake/sbatch ' \
				 f'--time_str {time_str} ' \
				 f'--queue {queue}'
	sbatch_total_path = sbatch_dir / f'sbatch-{queue}-queue.sh'
	with open(sbatch_total_path, 'w') as f:
		f.write(sbatch_cmd)

	print('#' * 40)
	print(f'For running jobs on the STAMPEDE2 {queue} queue:')
	print(f"All snakemake commands need to be executed "
		  f"were included in {sbatch_total_path}")
	print(f"You just need to run this script to "
		  f"map the whole library in {output_dir}. "
		  f"Note that this script will not return until all the mapping job finished. "
		  f"So you better run this script with nohup or in a screen.")
	print(f"You can also change "
		  f"the per job parameters in {script_path} "
		  f"or change the global parameters in {sbatch_total_path}")
	print(f"Read 'yap sbatch -h' if you want to have more options about sbatch. "
		  f"Alternatively, you can sbatch the commands in "
		  f"{outdir}/snakemake/sbatch/sbatch.sh by yourself, "
		  f"as long as they all get successfully executed.")
	print('#' * 40 + '\n')
	return

def prepare_run(output_dir, total_jobs=12, cores_per_job=10, total_memory_gb=None,
				name=None,fastq_server='local'):
	config = get_configuration(output_dir / 'mapping_config.ini')
	mode = config['mode']
	if mode.split('-')[0] in ['mc', 'm3c'] and cores_per_job < 4:
		raise ValueError(f'cores must >= 4 to run this pipeline.')
	elif mode.split('-')[0] in ['mct', '4m'] and cores_per_job < 10:
		raise ValueError(f'cores must >= 10 to run this pipeline.')

	output_dir = pathlib.Path(output_dir).absolute()
	if name is None:
		name = output_dir.name
	snakemake_dir = output_dir / 'snakemake'
	snakemake_dir.mkdir(exist_ok=True)

	# this is only some automatic code for ecker lab...
	# so conditioned by the host name
	try:
		host_name = os.environ['HOSTNAME']
	except KeyError:
		host_name = subprocess.run('hostname', stdout=subprocess.PIPE, encoding='utf8').stdout
		if not isinstance(host_name, str):
			host_name = 'unknown'
	prepare_qsub(name=name,
					snakemake_dir=snakemake_dir,
					total_jobs=total_jobs,
					cores_per_job=cores_per_job,
					total_memory_gb=total_memory_gb,
				 fastq_server=fastq_server)
	prepare_sbatch(name=name, snakemake_dir=snakemake_dir, queue='normal',total_memory_gb=total_memory_gb)
	prepare_sbatch(name=name, snakemake_dir=snakemake_dir, queue='skx-normal',total_memory_gb=total_memory_gb)
	prepare_sbatch(name=name, snakemake_dir=snakemake_dir, queue='shared',total_memory_gb=total_memory_gb)
	# else:
	#     script_path = write_qsub_commands(output_dir, cores_per_job, memory_gb_per_core, script_dir=snakemake_dir)
	#     print(f"All snakemake commands need to be executed were summarized in {script_path}")
	#     print(f"You need to execute them based on the computational environment you have "
	#           f"(e.g., use a job scheduler or run locally).")

	print(f"Once all commands are executed successfully, use 'yap summary' to generate final mapping summary.")
	return

def start_from_cell_fastq(output_dir, fastq_pattern, config_path,aligner='bismark',n_group=64,
						  n_jobs=64,total_memory_gb=None):
	output_dir = pathlib.Path(output_dir).absolute()
	if output_dir.exists():
		raise FileExistsError(f'Output dir {output_dir} already exist, please delete it or use another path.')
	output_dir.mkdir()
	subprocess.run(['cp', config_path, f'{output_dir}/mapping_config.ini'], check=True)
	stats_dir = output_dir / 'stats'
	stats_dir.mkdir(exist_ok=True)

	# parse fastq patterns
	fastq_paths = [pathlib.Path(p).absolute() for p in glob.glob(fastq_pattern)]
	r1_records = {}
	r2_records = {}
	for path in fastq_paths:
		*cell_id, suffix = path.name.split('-')
		cell_id = '-'.join(cell_id)
		if suffix == 'R1.fq.gz':
			if cell_id in r1_records:
				raise ValueError(f'Found duplicated cell ID: {cell_id}, '
								 f'File caused this error: {path}')
			r1_records[cell_id] = path
		elif suffix == 'R2.fq.gz':
			if cell_id in r2_records:
				raise ValueError(f'Found duplicated cell ID: {cell_id}, '
								 f'File caused this error: {path}')
			r2_records[cell_id] = path
		else:
			raise ValueError(
				f'Unable to parse read type. Expect file name ends with "-R1.fq.gz" or "-R2.fq.gz", '
				f'File caused this error: {path}'
			)
	fastq_df = pd.DataFrame({'R1Path': r1_records, 'R2Path': r2_records})

	# make symlink of fastq files, using dir structure of demultiplex
	# the cells are randomly grouped though, max group is 64
	groups = min(n_group, fastq_df.shape[0])
	for i, (cell_id, (r1_path, r2_path)) in enumerate(fastq_df.sample(fastq_df.shape[0]).iterrows()):
		group_id = i % groups
		fastq_dir = output_dir / f'Group{group_id}/fastq'
		fastq_dir.mkdir(exist_ok=True, parents=True)

		# make symlinks
		new_r1_path = fastq_dir / r1_path.name
		new_r1_path.symlink_to(r1_path)
		new_r2_path = fastq_dir / r2_path.name
		new_r2_path.symlink_to(r2_path)

	# prepare scripts
	if aligner.lower() == 'bismark':
		make_snakefile(output_dir)
	else:
		make_snakefile_hisat3n(output_dir)
	if total_memory_gb is None:
		total_memory_gb = 2 * n_jobs
	prepare_run(output_dir,cores_per_job=n_jobs,total_memory_gb=total_memory_gb)
	return
