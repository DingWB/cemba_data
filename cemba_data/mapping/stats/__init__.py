import pathlib
import subprocess
import pysam
from pysam import AlignmentFile
import pandas as pd
from papermill import execute_notebook, PapermillExecutionError
from ..pipelines import PACKAGE_DIR
from .utilities import parse_trim_fastq_stats, parse_trim_fastq_stats_mct, \
	parse_bismark_report, parse_deduplicate_stat, \
	generate_allc_stats

# mc
def mc_mapping_stats(output_dir=None,fastq_dir=None,mode="m3c",mc_stat_feature='CHN CGN CCC',
					mc_stat_alias='mCH mCG mCCC',num_upstr_bases=0):
	"""this may apply to single UID dir, so config is provided as parameter"""
	output_dir = pathlib.Path(output_dir).absolute()
	bam_dir = output_dir / 'bam'
	allc_dir = output_dir / 'allc'
	cell_stats = []
	cell_ids = [path.name.split('.')[0]
				for path in allc_dir.glob(f'*.allc.tsv.gz')]

	for cell_id in cell_ids:
		print(f'Parsing stats of {cell_id}.')
		total_stats = []
		for read_type in ['R1', 'R2']:
			if mode in ['4m', 'mct']:
				total_stats.append(
					parse_trim_fastq_stats_mct(
						fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.txt'))
			else:
				total_stats.append(
					parse_trim_fastq_stats(
						fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.tsv'))
			total_stats.append(
				parse_bismark_report(
					bam_dir / f'{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt'))
			total_stats.append(
				parse_deduplicate_stat(
					bam_dir / f'{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt'
				))
		cell_stats.append(pd.concat(total_stats))
	mapping_df = pd.DataFrame(cell_stats)
	mapping_df.index.name = 'cell_id'

	# add allc stats
	allc_df = generate_allc_stats(output_dir, mc_stat_feature,mc_stat_alias,num_upstr_bases)
	final_df = pd.concat([mapping_df, allc_df], sort=True, axis=1)
	return final_df


def mc_additional_cols(final_df):
	"""Additional columns for mC mapping summary"""
	final_df = final_df.copy()
	final_df['CellInputReadPairs'] = final_df['R1InputReads'].astype(int)  # == final_df['R2InputReads']
	if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
		cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
										for _, i in final_df.groupby('PCRIndex')])
		final_df['CellBarcodeRatio'] = cell_barcode_ratio

	final_df['FinalmCReads'] = final_df['R1FinalBismarkReads'] + final_df['R2FinalBismarkReads']
	return final_df

# m3c
def m3c_bam_unique_read_counts(bam_path, read_type_int):
	unique_reads = set()
	with AlignmentFile(bam_path) as bam:
		for read in bam:
			unique_reads.add(read.query_name.split(f'_{read_type_int}:N:0:')[0])
	return len(unique_reads)


def m3c_count_bams(bam_dir, cell_id, read_type):
	bam_path_dict = {
		f'{read_type}UniqueMappedReads': bam_dir / f'{cell_id}-{read_type}.two_mapping.filter.bam',
		f'{read_type}DeduppedReads': bam_dir / f'{cell_id}-{read_type}.two_mapping.deduped.bam',
	}
	read_counts = {name: m3c_bam_unique_read_counts(path, 1 if read_type == 'R1' else 2)
				   for name, path in bam_path_dict.items()}
	return pd.Series(read_counts, name=cell_id)


def m3c_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases):
	"""this may apply to single UID dir, so config is provided as parameter"""
	output_dir = pathlib.Path(output_dir).absolute()
	#fastq_dir = output_dir / 'fastq'
	bam_dir = output_dir / 'bam'
	hic_dir = output_dir / 'hic'
	cell_stats = []
	cell_ids = [path.name.split('.')[0]
				for path in bam_dir.glob('*.3C.sorted.bam')]

	for cell_id in cell_ids:
		total_stats = []  # list of series
		for read_type in ['R1', 'R2']:
			# fastq reads
			if mode in ['4m', 'mct']:
				total_stats.append(
					parse_trim_fastq_stats_mct(
						fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.txt'))
			else:
				total_stats.append(
					parse_trim_fastq_stats(
						fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.tsv'))
			# bam reads
			total_stats.append(
				m3c_count_bams(bam_dir, cell_id, read_type)
			)
		# contacts
		contact_counts = pd.read_csv(hic_dir / f'{cell_id}.3C.contact.tsv.counts.txt',
									 header=None, index_col=0).squeeze()
		contact_counts.name = cell_id
		total_stats.append(contact_counts)

		cell_stats.append(pd.concat(total_stats))
	total_df = pd.DataFrame(cell_stats)

	# add allc stats
	allc_df = generate_allc_stats(output_dir, mc_stat_feature,mc_stat_alias,num_upstr_bases)
	final_df = pd.concat([total_df, allc_df], sort=True, axis=1)
	return final_df


def m3c_additional_cols(final_df):
	final_df['FinalmCReads'] = final_df['R1DeduppedReads'] + final_df['R2DeduppedReads']
	final_df['CellInputReadPairs'] = final_df['R1InputReads']
	# use % to be consistent with others
	final_df['R1MappingRate'] = final_df['R1UniqueMappedReads'] / final_df['R1TrimmedReads'] * 100
	final_df['R2MappingRate'] = final_df['R2UniqueMappedReads'] / final_df['R2TrimmedReads'] * 100
	final_df['R1DuplicationRate'] = (1 - final_df['R1DeduppedReads'] / final_df['R1UniqueMappedReads']) * 100
	final_df['R2DuplicationRate'] = (1 - final_df['R2DeduppedReads'] / final_df['R2UniqueMappedReads']) * 100

	if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
		cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
										for _, i in final_df.groupby('PCRIndex')])
		final_df['CellBarcodeRatio'] = cell_barcode_ratio

	final_df['TotalContacts'] = final_df[
		['CisShortContact', 'CisLongContact', 'TransContact']].sum(axis=1)
	final_df['CisShortRatio'] = final_df['CisShortContact'] / final_df['TotalContacts']
	final_df['CisLongRatio'] = final_df['CisLongContact'] / final_df['TotalContacts']
	final_df['TransRatio'] = final_df['TransContact'] / final_df['TotalContacts']
	return final_df

# mct
def _count_reads_by_rg_in_star_bam(bam_path):
	try:
		bam = pysam.AlignmentFile(bam_path)
	except ValueError:
		# empty bam file
		return

	cell_read_counts = {cell['ID']: 0 for cell in bam.header['RG']}

	for read in bam:
		cell = read.get_tag('RG')
		cell_read_counts[cell] += 1
	read_count = pd.Series(cell_read_counts, name='Reads')
	read_count.index.name = 'cell_id'
	return read_count


def summary_rna_mapping(output_dir):
	output_dir = pathlib.Path(output_dir)

	# summarize read counts for each cell before filter by mC rate
	total_star_mapped_reads = _count_reads_by_rg_in_star_bam(output_dir / 'rna_bam/TotalRNAAligned.filtered.bam')

	# feature count summary
	total_counts = pd.read_csv(output_dir / 'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary',
							   sep='\t', index_col=0).T
	total_counts.index = total_counts.index.map(lambda i: i.split(':')[-1])
	feature_count_summary = total_counts[['Assigned']].copy()
	feature_count_summary['FinalRNAReads'] = total_counts.sum(axis=1)
	feature_count_summary.columns = ['FinalCountedReads', 'FinalRNAReads']

	total_rna_stat = feature_count_summary.copy()
	total_rna_stat['RNAUniqueMappedReads'] = total_star_mapped_reads
	total_rna_stat['SelectedRNAReadsRatio'] = total_rna_stat['FinalRNAReads'] / total_rna_stat['RNAUniqueMappedReads']
	total_rna_stat.index.name = 'cell_id'
	return total_rna_stat

def summarize_select_dna_reads(output_dir,mc_rate_max_threshold,dna_cov_min_threshold):
	bam_dir = pathlib.Path(output_dir) / 'bam'
	mc_rate_max_threshold = float(mc_rate_max_threshold) * 100
	cov_min_threshold = float(dna_cov_min_threshold)

	records = []
	select_dna_reads_stat_list = list(bam_dir.glob('*.reads_profile.csv'))
	for path in select_dna_reads_stat_list:
		try:
			_df = pd.read_csv(path)
		except pd.errors.EmptyDataError:
			# means the bam file is empty
			continue

		cell_id = path.name.split('.')[0]
		if cell_id.endswith('-R1') or cell_id.endswith('-R2'):
			# select DNA preformed in R1 R2 separately:
			cell_id = cell_id[:-3]
		_df['cell_id'] = cell_id
		_df['mc_rate_max_threshold'] = mc_rate_max_threshold
		_df['cov_min_threshold'] = cov_min_threshold
		records.append(_df)
	total_stats_df = pd.concat(records)

	selected_reads = total_stats_df[
		(total_stats_df['cov'] >= cov_min_threshold)
		& (total_stats_df['mc_frac'] < mc_rate_max_threshold)]

	selected_reads = selected_reads.groupby('cell_id')['count'].sum()
	selected_ratio = selected_reads / total_stats_df.groupby('cell_id')['count'].sum()
	final_stat = pd.DataFrame({'FinalDNAReads': selected_reads, 'SelectedDNAReadsRatio': selected_ratio})
	final_stat.index.name = 'cell_id'
	return final_stat

def mct_mapping_stats(output_dir, fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases,
					  mc_rate_max_threshold,dna_cov_min_threshold):
	"""this may apply to single UID dir, so config is provided as parameter"""
	mc_stats_df = mc_mapping_stats(output_dir, fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases)
	select_dna_stats_df = summarize_select_dna_reads(output_dir, mc_rate_max_threshold,dna_cov_min_threshold)
	rna_stats_df = summary_rna_mapping(output_dir)
	final_df = pd.concat([mc_stats_df, select_dna_stats_df, rna_stats_df], axis=1)
	return final_df


def aggregate_feature_counts(output_dir):
	output_dir = pathlib.Path(output_dir)
	cell_data = []

	count_paths = list(output_dir.glob('*/bam/*.feature_count.tsv'))
	if len(count_paths) == 0:
		return

	data = None
	for path in count_paths:
		data = pd.read_csv(path, sep='\t', index_col=0, comment='#')
		cell_data.append(data.iloc[:, 5:])
	cell_data = pd.concat(cell_data, axis=1, sort=True)
	cell_data.columns = cell_data.columns.str.split(':').str[1]

	# all count table should have the same info, so read the last one
	# chr, start, end, strand, length
	gene_info = data.iloc[:, :5]
	with pd.HDFStore(output_dir / 'TotalRNAData.h5', mode='w', complevel=5) as hdf:
		hdf['data'] = cell_data.T  # cell by gene
		hdf['gene'] = gene_info
		hdf['stats'] = pd.DataFrame({'GenesDetected': (cell_data > 0).sum()})
	return

def mct_additional_cols(final_df, output_dir):
	final_df = final_df.copy()
	input_read_col='InputReadPairs' if 'R1InputReads' not in final_df.columns.tolist() else 'R1InputReads'
	final_df['CellInputReadPairs'] = final_df[input_read_col].astype(int)  # == final_df['R2InputReads']
	if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
		cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
										for _, i in final_df.groupby('PCRIndex')])
		final_df['CellBarcodeRatio'] = cell_barcode_ratio

	stats = pd.read_hdf(output_dir / 'TotalRNAData.h5', key='stats')
	final_df['GenesDetected'] = stats['GenesDetected']

	# calculate some mCT specific ratios
	final_df['DNAReadsYield'] = final_df['FinalDNAReads'] / (
			final_df['CellInputReadPairs'] * 2)
	final_df['RNAReadsYield'] = final_df['FinalRNAReads'] / final_df[
		'CellInputReadPairs']
	if 'R1FinalBismarkReads' in final_df.columns:
		final_df['RNA/(DNA+RNA)'] = final_df['FinalRNAReads'].fillna(0) / (
				final_df['R1FinalBismarkReads'].fillna(0) + 1)
	else: #20240321 added, not sure whether it is correct
		final_df['RNA/(DNA+RNA)'] = final_df['FinalRNAReads'].fillna(0) / (
			final_df['FinalRNAReads'].fillna(0)+ final_df['FinalDNAReads'].fillna(0) + 1)
	return final_df

# 4m
def _4m_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases,
					mc_rate_max_threshold,dna_cov_min_threshold):
	"""this may apply to single UID dir, so config is provided as parameter"""
	m3c_stats_df = m3c_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases)
	select_dna_stats_df = summarize_select_dna_reads(output_dir, mc_rate_max_threshold,dna_cov_min_threshold)
	rna_stats_df = summary_rna_mapping(output_dir)
	final_df = pd.concat([m3c_stats_df, select_dna_stats_df, rna_stats_df], axis=1)
	return final_df


def _4m_additional_cols(final_df, output_dir):
	final_df = final_df.copy()
	final_df['CellInputReadPairs'] = final_df['R1InputReads'].astype(int)
	if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
		cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
										for _, i in final_df.groupby('PCRIndex')])
		final_df['CellBarcodeRatio'] = cell_barcode_ratio

	# snm3C part
	final_df['FinalmCReads'] = final_df['R1DeduppedReads'] + final_df['R2DeduppedReads']
	# use % to be consistent with others
	final_df['R1MappingRate'] = final_df['R1UniqueMappedReads'] / final_df['R1TrimmedReads'] * 100
	final_df['R2MappingRate'] = final_df['R2UniqueMappedReads'] / final_df['R2TrimmedReads'] * 100
	final_df['R1DuplicationRate'] = (1 - final_df['R1DeduppedReads'] / final_df['R1UniqueMappedReads']) * 100
	final_df['R2DuplicationRate'] = (1 - final_df['R2DeduppedReads'] / final_df['R2UniqueMappedReads']) * 100
	final_df['TotalContacts'] = final_df[
		['CisShortContact', 'CisLongContact', 'TransContact']].sum(axis=1)
	final_df['CisShortRatio'] = final_df['CisShortContact'] / final_df['TotalContacts']
	final_df['CisLongRatio'] = final_df['CisLongContact'] / final_df['TotalContacts']
	final_df['TransRatio'] = final_df['TransContact'] / final_df['TotalContacts']

	# snmCT part
	stats = pd.read_hdf(output_dir / 'TotalRNAData.h5', key='stats')
	final_df['GenesDetected'] = stats['GenesDetected']
	# calculate some mCT specific ratios
	final_df['DNAReadsYield'] = final_df['FinalDNAReads'] / (
			final_df['CellInputReadPairs'] * 2)
	final_df['RNAReadsYield'] = final_df['FinalRNAReads'] / final_df[
		'CellInputReadPairs']
	final_df['RNA/(DNA+RNA)'] = final_df['FinalRNAReads'].fillna(0) / (
			final_df['R1DeduppedReads'].fillna(0) + 1)
	return final_df

# plate info
def _parse_cell_id_v1(cell_id):
	plate1, plate2, pcr_index, random_index = cell_id.split('-')
	if random_index.upper() in {'AD001', 'AD002', 'AD004', 'AD006'}:
		plate = plate1
	else:
		plate = plate2
	# 96 pos
	col96 = int(pcr_index[1:]) - 1
	row96 = ord(pcr_index[0]) - 65  # convert A-H to 0-8
	# 384 pos
	ad_index_384_dict = {
		'AD001': (0, 0),
		'AD002': (0, 1),
		'AD004': (1, 0),
		'AD006': (1, 1),
		'AD007': (0, 0),
		'AD008': (0, 1),
		'AD010': (1, 0),
		'AD012': (1, 1)
	}
	col384 = 2 * col96 + ad_index_384_dict[random_index][0]
	row384 = 2 * row96 + ad_index_384_dict[random_index][1]
	record = pd.Series({
		'Plate': plate,
		'PCRIndex': pcr_index,
		'RandomIndex': random_index,
		'Col384': col384,
		'Row384': row384
	})
	return record


def _parse_cell_id_v2(cell_id):
	plate, multiplex_group, pcr_index, random_index = cell_id.split('-')
	# 384 pos
	col384 = int(random_index[1:]) - 1
	row384 = ord(random_index[0]) - 65  # convert A-P to 0-23
	record = pd.Series({
		'Plate': plate,
		'PCRIndex': pcr_index,
		'MultiplexGroup': multiplex_group,
		'RandomIndex': random_index,
		'Col384': col384,
		'Row384': row384
	})
	return record


def get_plate_info(cell_ids, barcode_version):
	if barcode_version == 'V1':
		func = _parse_cell_id_v1
	else:
		func = _parse_cell_id_v2
	try:
		plate_info = pd.DataFrame([func(cell_id) for cell_id in cell_ids],
								  index=cell_ids)
	except Exception:
		print('Errors occur during parsing the plate info, this happens '
			  'when the input FASTQ file name is not generated by yap. '
			  'The `yap summary` also can not generate html report due to missing the plate info. '
			  'In this case, you need to add the plateinfo by yourself in order to make the plate view plots. '
			  'These information is not necessary for following analysis though.')
		plate_info = pd.DataFrame([], index=cell_ids)
	return plate_info

# main function
def mapping_stats(output_dir=None,fastq_dir=None,mode='m3c',barcode_version='V2',
				mc_stat_feature='CHN CGN CCC',mc_stat_alias='mCH mCG mCCC',num_upstr_bases=0,
				mc_rate_max_threshold=0.5,dna_cov_min_threshold=3):
	"""This is UID level mapping summary, the config file is in parent dir"""
	output_dir = pathlib.Path(output_dir).absolute()
	if fastq_dir is None:
		fastq_dir = output_dir / 'fastq'
	else:
		fastq_dir=pathlib.Path(fastq_dir).absolute()
	if mode == 'mc':
		final_df = mc_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases)
	elif mode == 'mct':
		final_df = mct_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases,
					  mc_rate_max_threshold,dna_cov_min_threshold)
	elif mode == 'm3c':
		final_df = m3c_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases)
	elif mode == '4m':
		final_df = _4m_mapping_stats(output_dir,fastq_dir,mode,mc_stat_feature,mc_stat_alias,num_upstr_bases,
					mc_rate_max_threshold,dna_cov_min_threshold)
	else:
		raise ValueError

	# plate info, which is tech independent.
	_plate_info = get_plate_info(final_df.index, barcode_version=barcode_version)
	final_df = pd.concat([_plate_info, final_df], axis=1)

	# save
	final_df.to_csv(output_dir / 'MappingSummary.csv.gz')
	return


def final_summary(output_dir, cleanup=True, notebook=None,mode='m3c'):
	output_dir = pathlib.Path(output_dir).absolute()
	path_to_remove = []

	# Before running summary,
	# first make sure all the UID dir having Snakefile also has mapping summary (means successful)
	snakefile_list = list(output_dir.glob('*/Snakefile'))
	summary_paths = []
	missing_summary_dirs = []
	for path in snakefile_list:
		uid_dir = path.parent
		summary_path = uid_dir / 'MappingSummary.csv.gz'
		if summary_path.exists():
			summary_paths.append(summary_path)
		else:
			missing_summary_dirs.append(uid_dir)

	if len(missing_summary_dirs) != 0:
		print('These sub dir missing MappingSummary files:')
		for p in missing_summary_dirs:
			print(p)
		raise FileNotFoundError(f'Note that all sub dir should be successfully mapped '
								f'before generating final summary. \n'
								f'The MappingSummary.csv.gz is the final target file of snakefile in {path}. \n'
								f'Run the corresponding snakemake command again to retry mapping.\n'
								f'The snakemake commands can be found in output_dir/snakemake/*/snakemake_cmd.txt')

	# aggregate mapping summaries
	total_mapping_summary = pd.concat([pd.read_csv(path, index_col=0)
									   for path in summary_paths])
	total_mapping_summary_path = output_dir / 'stats/MappingSummary.csv.gz'

	# if this is mct, aggregate all the gene counts
	if mode in ['mct', '4m']:
		aggregate_feature_counts(output_dir)

	# add additional columns based on some calculation
	if mode == 'mc':
		total_mapping_summary = mc_additional_cols(total_mapping_summary)
	elif mode == 'mct':
		total_mapping_summary = mct_additional_cols(total_mapping_summary, output_dir=output_dir)
	elif mode == 'm3c':
		total_mapping_summary = m3c_additional_cols(total_mapping_summary)
	elif mode == '4m':
		total_mapping_summary = _4m_additional_cols(total_mapping_summary, output_dir=output_dir)
	else:
		raise

	# save total mapping summary
	total_mapping_summary.to_csv(total_mapping_summary_path)

	# add .snakemake files to deletion
	snakemake_hiding_dirs = list(output_dir.glob('*/.snakemake'))
	path_to_remove += snakemake_hiding_dirs

	# add temp dir in the bam dirs to deletion
	mapping_temp_dirs = list(output_dir.glob('*/bam/temp'))
	path_to_remove += mapping_temp_dirs

	# write a ALLC path file for generating MCDS
	allc_paths = pd.Series({path.name.split('.')[0]: str(path)
							for path in output_dir.glob('*/allc/*tsv.gz')})
	allc_paths.to_csv(output_dir / 'stats/AllcPaths.tsv', sep='\t', header=False)

	if 'Plate' in total_mapping_summary.columns:  # only run notebook when plate info exist
		# run summary notebook
		nb_path = output_dir / 'stats/MappingSummary.ipynb'
		try:
			if notebook is None:
				template_notebook = PACKAGE_DIR / f'files/mapping_summary_template/{mode}_template.ipynb'
			else:
				template_notebook = str(notebook)
			print(f'Using notebook template from {template_notebook}')
			print('Executing summary plotting notebook...')
			execute_notebook(
				input_path=str(template_notebook),
				output_path=str(nb_path),
				parameters=dict(output_dir=str(output_dir))
			)
			print('Summary notebook successfully executed. Exporting HTML...')
			subprocess.run(['jupyter', 'nbconvert', '--to', 'html', str(nb_path)])
			print(f'See the summary plots here: {str(nb_path)[:-5]}html')
			print(f'Or customize the summary plots here: {nb_path}')
		except PapermillExecutionError:
			print(f'Ops, summary plotting notebook got some error, check the information in {nb_path}')
			cleanup = False

	# delete
	if cleanup:
		print('Clean up snakemake log (might take several minutes) ...')
		for path in path_to_remove:
			subprocess.run(['rm', '-rf', str(path)], check=True)
	return


