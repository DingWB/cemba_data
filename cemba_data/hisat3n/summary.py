import os

from .stats_parser import *


def snmc_summary(outname="MappingSummary.csv.gz",indir="."):
	"""
	Generate snmC pipeline MappingSummary.csv.gz and save into cwd

	Returns
	-------
	pd.DataFrame
	"""
	all_stats = []

	# fastq trimming stats
	df = parse_single_stats_set(path_pattern=indir+'/fastq/*.trimmed.stats.txt',
								parser=cell_parser_cutadapt_trim_stats,indir=indir)
	all_stats.append(df)

	# hisat-3n mapping
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna_summary.txt',
								parser=cell_parser_hisat_summary,indir=indir)
	all_stats.append(df)

	# uniquely mapped reads dedup
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.unique_align.deduped.matrix.txt',
								parser=cell_parser_picard_dedup_stat,
								prefix='UniqueAlign',indir=indir)
	all_stats.append(df)

	# multi mapped reads dedup
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.multi_align.deduped.matrix.txt',
								parser=cell_parser_picard_dedup_stat,
								prefix='MultiAlign',indir=indir)
	all_stats.append(df)

	# allc count
	df = parse_single_stats_set(path_pattern=indir+'/allc/*.allc.tsv.gz.count.csv',
								parser=cell_parser_allc_count,indir=indir)
	all_stats.append(df)

	# concatenate all stats
	all_stats = pd.concat(all_stats, axis=1)
	all_stats.index.name = 'cell'
	if all_stats.shape[0] > 0:
		all_stats.to_csv(outname)
	else:
		print(f'Nothing in {outname}')
	return all_stats


def snmct_summary(outname="MappingSummary.csv.gz",indir="."):
	"""
	Generate snmCT pipeline MappingSummary.csv.gz and save into cwd

	Returns
	-------
	pd.DataFrame
	"""
	all_stats = []

	# fastq trimming stats
	df = parse_single_stats_set(path_pattern=indir+'/fastq/*.trimmed.stats.txt',
								parser=cell_parser_cutadapt_trim_stats,indir=indir)
	all_stats.append(df)

	# hisat-3n DNA mapping
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna_summary.txt',
								parser=cell_parser_hisat_summary, prefix='DNA',indir=indir)
	all_stats.append(df)

	# hisat-3n RNA mapping
	df = parse_single_stats_set(path_pattern=indir+'/rna_bam/*.hisat3n_rna_summary.txt',
								parser=cell_parser_hisat_summary, prefix='RNA',indir=indir)
	all_stats.append(df)

	# uniquely mapped reads dedup
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.unique_align.deduped.matrix.txt',
								parser=cell_parser_picard_dedup_stat,
								prefix='DNAUniqueAlign',indir=indir)
	all_stats.append(df)

	# multi mapped reads dedup
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.multi_align.deduped.matrix.txt',
								parser=cell_parser_picard_dedup_stat,
								prefix='DNAMultiAlign',indir=indir)
	all_stats.append(df)

	# uniquely mapped dna reads selection
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv',
								parser=cell_parser_reads_mc_frac_profile,
								prefix='UniqueAlign',indir=indir)
	all_stats.append(df)

	# multi mapped dna reads selection
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna.multi_align.deduped.dna_reads.reads_mch_frac.csv',
								parser=cell_parser_reads_mc_frac_profile,
								prefix='MultiAlign',indir=indir)
	all_stats.append(df)

	# uniquely mapped rna reads selection
	df = parse_single_stats_set(path_pattern=indir+'/rna_bam/*.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv',
								parser=cell_parser_reads_mc_frac_profile,indir=indir)
	all_stats.append(df)

	# allc count
	df = parse_single_stats_set(path_pattern=indir+'/allc/*.allc.tsv.gz.count.csv',
								parser=cell_parser_allc_count,indir=indir)
	all_stats.append(df)

	# feature count
	df = parse_single_stats_set(path_pattern=indir+'/rna_bam/*.feature_count.tsv.summary',
								parser=cell_parser_feature_count_summary,indir=indir)
	all_stats.append(df)

	# concatenate all stats
	all_stats = pd.concat(all_stats, axis=1)
	all_stats.index.name = 'cell'
	if all_stats.shape[0] > 0:
		all_stats.to_csv(outname)
	else:
		print(f'Nothing in {outname}')
	return all_stats


def snm3c_summary(outname="MappingSummary.csv.gz",indir="."):
	"""
	Generate snm3C pipeline MappingSummary.csv.gz and save into cwd

	Returns
	-------
	pd.DataFrame
	"""
	print(f"CWD: {os.getcwd()}")
	print(f"indir: {indir}, outname: {outname}")
	all_stats = []

	# fastq trimming stats
	df = parse_single_stats_set(path_pattern=indir+'/fastq/*.trimmed.stats.txt',
								parser=cell_parser_cutadapt_trim_stats,indir=indir)
	all_stats.append(df)

	# hisat-3n mapping PE
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna_summary.txt',
								parser=cell_parser_hisat_summary,indir=indir)
	all_stats.append(df)

	# hisat-3n mapping split-reads SE
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna_split_reads_summary.R1.txt',
								parser=cell_parser_hisat_se_summary, prefix='R1SplitReads',
								indir=indir) #single end summary
	all_stats.append(df)

	df = parse_single_stats_set(path_pattern=indir+'/bam/*.hisat3n_dna_split_reads_summary.R2.txt',
								parser=cell_parser_hisat_se_summary, prefix='R2SplitReads',
								indir=indir)
	all_stats.append(df)

	# uniquely mapped reads dedup
	df = parse_single_stats_set(path_pattern=indir+'/bam/*.all_reads.deduped.matrix.txt',
								parser=cell_parser_picard_dedup_stat, prefix='UniqueAlign',
								indir=indir)
	all_stats.append(df)

	# call chromatin contacts
	df = parse_single_stats_set(path_pattern=indir+'/hic/*.all_reads.contact_stats.csv',
								parser=cell_parser_call_chromatin_contacts,
								indir=indir)
	all_stats.append(df)

	# allc count
	df = parse_single_stats_set(path_pattern=indir+'/allc/*.allc.tsv.gz.count.csv',
								parser=cell_parser_allc_count,indir=indir)
	all_stats.append(df)
	# concatenate all stats
	all_stats = pd.concat(all_stats, axis=1)
	all_stats.index.name = 'cell'
	if all_stats.shape[0] > 0:
		all_stats.to_csv(outname)
	else:
		print(f'Nothing in {outname}')
	return all_stats
