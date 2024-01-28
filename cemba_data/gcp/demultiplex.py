import os,sys
import pandas as pd
from snakemake.io import glob_wildcards

def get_fastq_info(fq_dir,outdir,fq_ext,run_on_gcp):
	if os.path.exists("fastq_info.txt"):
		df=pd.read_csv("fastq_info.txt",sep='\t')
		return df
	#For example: 220517-AMB-mm-na-snm3C_seq-NovaSeq-pe-150-WT-AMB_220510_8wk_12D_13B_2_P3-1-A11_S7_L001_R1_001.fastq.gz
	indirs,prefixes,plates,multiple_groups,primer_names,pns,lanes,\
	read_types,suffixes=glob_wildcards(os.path.join(str(fq_dir),\
			"{indir}/{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")) \
			if not run_on_gcp else \
			GS.glob_wildcards(fq_dir+"/{indir}/{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")
	indirs=[fq_dir + '/' + d for d in indirs]

	if len(indirs)==0: # depth=1
		prefixes,plates,multiple_groups,primer_names,pns,lanes,\
		read_types,suffixes=glob_wildcards(os.path.join(str(fq_dir),\
				"{prefixes}-{plates}-{multiplex_groups}-{primer_names}_{pns}_{lanes}_{read_types}_{suffixes}."+f"{fq_ext}.gz")) \
				if not run_on_gcp else \
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
	df.to_csv("fastq_info.txt",sep='\t',index=False)
	return df

def get_lanes_info(outdir,barcode_version):
	#  uid={plate}-{multiplex_group}-{primer_name}
	if os.path.exists("lane_info.txt"):
		df1=pd.read_csv("lane_info.txt",sep='\t')
		df1.fastq_path = df1.fastq_path.apply(lambda x: eval(x))
		return df1
	uids,plates,multiple_groups,primer_names,lanes,index_names,read_types=\
			glob_wildcards(os.path.join(outdir,
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
	df['fastq_path']=df.apply(
		lambda row:os.path.join(outdir,row.uid,"lanes",'-'.join(
			row.loc[['uid','lane','index_name','read_type']].map(str).tolist())+\
							".fq.gz"),axis=1)
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
	df1['fastq_out'] = df1.apply(lambda row:
								 os.path.join(
									 outdir, row.uid, "fastq", '-'.join(
										 row.loc[['uid', 'index_name', 'read_type']].map(
								  str).tolist()) + ".fq.gz"), axis=1)
	df1.to_csv("lane_info.txt",sep='\t',index=False)
	return df1