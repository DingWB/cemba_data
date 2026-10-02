import pathlib
import pandas as pd

from ...utilities import parse_mc_pattern


def parse_trim_fastq_stats(stat_path):
    # example trim fastq stats
    """
status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
0	OK	1490	213724	0	0	0	1490	4	0	213712
1	status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
2	OK	1490	213712	0	0	0	1482	0	1300	182546
"""
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    trim_stats = pd.read_csv(stat_path, sep='\t')
    trim_stats = trim_stats.iloc[[0, 2], :].reset_index()  # skip the duplicated title row

    data = pd.Series({
        f'{read_type}InputReads': trim_stats['in_reads'][0],
        f'{read_type}InputReadsBP': trim_stats['in_bp'][0],
        f'{read_type}WithAdapters': trim_stats['w/adapters'][0],
        f'{read_type}QualTrimBP': trim_stats['qualtrim_bp'][1],
        f'{read_type}TrimmedReads': trim_stats['out_reads'][1],
        f'{read_type}TrimmedReadsBP': trim_stats['out_bp'][1],
        f'{read_type}TrimmedReadsRate': int(trim_stats['out_reads'][1]) / int(trim_stats['in_reads'][0])
    }, name=cell_id)
    return data


def parse_trim_fastq_stats_mct(stat_path):
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    with open(stat_path) as f:
        lines = f.readlines()
        adapter_str = ''.join(lines[:-2])
        trim_lines = lines[-2:]

    # adapter counts
    total_dict = {}
    for line in adapter_str.replace('===\n\n', '; ').replace('=== Adapter ', 'Adapter: ').split('\n'):
        if line.startswith('=== Summary'):
            total_dict[f'{read_type}InputReads'] = int(line.strip('\n').split(' ')[-1].replace(',', ''))
        if line.startswith('Adapter: '):
            line_list = line.split('; ')
            line_dict = {}
            for l in line_list:
                k, v = l.split(': ')
                line_dict[k] = v
            name = line_dict.pop('Adapter').strip(' ')
            total_dict[f'{read_type}With{name}'] = int(line_dict['Trimmed'][:-5])
    data = pd.Series(total_dict, name=cell_id)

    # add trimmed counts, the last two rows in tsv format, the same as normal mc
    trim_data = pd.DataFrame([
        line.strip('\n').split('\t') for line in trim_lines
    ]).T.set_index(0)[1]
    data[f'{read_type}QualTrimBP'] = int(trim_data['qualtrim_bp'])
    data[f'{read_type}TrimmedReads'] = int(trim_data['out_reads'])
    data[f'{read_type}TrimmedReadsBP'] = int(trim_data['out_bp'])
    data[f'{read_type}TrimmedReadsRate'] = data[f'{read_type}TrimmedReads'] / data[f'{read_type}InputReads']
    return data


def parse_bismark_report(stat_path):
    """
    parse bismark output
    """
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    term_dict = {
        'Number of alignments with a unique best hit from the different alignments': f'{read_type}UniqueMappedReads',
        'Mapping efficiency': f'{read_type}MappingRate',
        'Sequences with no alignments under any condition': f'{read_type}UnmappedReads',
        'Sequences did not map uniquely': f'{read_type}UnuniqueMappedReads',
        'CT/CT': f'{read_type}OT',
        'CT/GA': f'{read_type}OB',
        'GA/CT': f'{read_type}CTOT',
        'GA/GA': f'{read_type}CTOB',
        "Total number of C's analysed": f'{read_type}TotalC',
        'C methylated in CpG context': f'{read_type}TotalmCGRate',
        'C methylated in CHG context': f'{read_type}TotalmCHGRate',
        'C methylated in CHH context': f'{read_type}TotalmCHHRate'}

    with open(stat_path) as rep:
        report_dict = {}
        for line in rep:
            try:
                start, rest = line.split(':')
            except ValueError:
                continue  # more or less than 2 after split
            try:
                report_dict[term_dict[start]] = rest.strip().split('\t')[0].strip('%')
            except KeyError:
                pass
    return pd.Series(report_dict, name=cell_id)


def parse_deduplicate_stat(stat_path):
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    try:
        dedup_result_series = pd.read_csv(stat_path, comment='#', sep='\t').T[0]
        rename_dict = {
            'UNPAIRED_READS_EXAMINED': f'{read_type}MAPQFilteredReads',
            'UNPAIRED_READ_DUPLICATES': f'{read_type}DuplicatedReads',
            'PERCENT_DUPLICATION': f'{read_type}DuplicationRate'
        }
        dedup_result_series = dedup_result_series.loc[rename_dict.keys()].rename(rename_dict)

        dedup_result_series[f'{read_type}FinalBismarkReads'] = dedup_result_series[f'{read_type}MAPQFilteredReads'] - \
                                                               dedup_result_series[f'{read_type}DuplicatedReads']
        dedup_result_series.name = cell_id
    except pd.errors.EmptyDataError:
        # if a BAM file is empty, picard matrix is also empty
        dedup_result_series = pd.Series({f'{read_type}MAPQFilteredReads': 0,
                                         f'{read_type}DuplicatedReads': 0,
                                         f'{read_type}FinalBismarkReads': 0,
                                         f'{read_type}DuplicationRate': 0}, name=cell_id)
    return dedup_result_series


def generate_allc_stats(output_dir, mc_stat_feature,mc_stat_alias,num_upstr_bases, mc_format='allc'):
    output_dir = pathlib.Path(output_dir).absolute()
    # Select the methylation count files by the chosen output format:
    #   mc_format='cz' -> cytozip bam_to_cz counts (cz/<cell>.cz.count.csv)
    #   otherwise      -> ALLCools bam-to-allc counts (allc/<cell>.*count.csv)
    # Both share the same context-indexed mc/cov/genome_cov layout. Reading
    # only the chosen format avoids picking up stale count files from a
    # previous run with a different mc_format.
    if str(mc_format).lower() == 'cz':
        allc_list = []  # lambda stats come from the lambda_* columns of cz count.csv
        allc_stats_dict = {p.name.split('.')[0]: p for p in output_dir.glob('cz/*count.csv')}
    else:
        allc_list = list(output_dir.glob('allc/*tsv.gz'))
        allc_stats_dict = {p.name.split('.')[0]: p for p in output_dir.glob('allc/*count.csv')}

    # no methylation count files at all (e.g. run failed) -> nothing to add
    if len(allc_stats_dict) == 0:
        import warnings
        warnings.warn(f'no {mc_format} count tables found under {output_dir}; '
                      f'MappingSummary will lack mC / Lambda columns')
        return pd.DataFrame()

    # patterns = config['mc_stat_feature'].split(' ')
    # patterns_alias = config['mc_stat_alias'].split(' ')
    patterns = mc_stat_feature.split(' ')
    patterns_alias = mc_stat_alias.split(' ')
    pattern_translate = {k: v for k, v in zip(patterns, patterns_alias)}

    # real all cell stats
    total_stats = []
    for cell_id, path in allc_stats_dict.items():
        allc_stat = pd.read_csv(path, index_col=0)
        allc_stat['cell_id'] = cell_id
        total_stats.append(allc_stat)
    total_stats = pd.concat(total_stats)
    cell_genome_cov = pd.Series(total_stats.set_index('cell_id')['genome_cov'].to_dict())
    # aggregate into patterns
    cell_records = []
    for pattern in pattern_translate.keys():
        contexts = parse_mc_pattern(pattern)
        pattern_stats = total_stats[total_stats.index.isin(contexts)]
        cell_level_data = pattern_stats.groupby('cell_id')[['mc', 'cov']].sum()
        cell_level_data['frac'] = cell_level_data['mc'] / cell_level_data['cov']

        # prettify col name
        _pattern = pattern_translate[pattern]
        cell_level_data = cell_level_data.rename(
            columns={'frac': f'{_pattern}Frac',
                     'mc': f'{_pattern}mC',
                     'cov': f'{_pattern}Cov'})
        cell_records.append(cell_level_data)
    final_df = pd.concat(cell_records, axis=1, sort=True).reindex(allc_stats_dict.keys())
    final_df['GenomeCov'] = cell_genome_cov

    # add lambda DNA stats (Lambda{CA,CC,CT,CH,CY,CG}{mC,Cov,Frac})
    if str(mc_format).lower() == 'cz':
        lambda_frac = get_cz_lambda_frac(list(allc_stats_dict.values()), num_upstr_bases)
    else:
        lambda_frac = get_allc_lambda_frac(allc_list, num_upstr_bases)
    for col, data in lambda_frac.items():
        final_df[col] = data

    final_df.index.name = 'cell_id'
    return final_df


def _lambda_table(paths, parser):
    from ...hisat3n.stats_parser import _lambda_records
    rows = {}
    for path in paths:
        cell = pathlib.Path(path).name.split('.')[0]
        record = parser(path)
        if record.empty:  # no chrL in this cell
            zeros = dict.fromkeys(['CA', 'CC', 'CT', 'CG'], 0)
            record = _lambda_records(zeros, dict(zeros), cell)
        rows[cell] = record
    return pd.DataFrame(rows).T


def get_allc_lambda_frac(allc_list, num_upstr_bases=None):
    """Lambda (chrL) stats of tabix-indexed ALLC files; the C position is inferred from each count table."""
    from ...hisat3n.stats_parser import cell_parser_allc_lambda
    return _lambda_table(allc_list, cell_parser_allc_lambda)


def get_cz_lambda_frac(count_list, num_upstr_bases=None):
    """Same as get_allc_lambda_frac, from the lambda_mc/lambda_cov columns of cytozip <cell>.cz.count.csv."""
    from ...hisat3n.stats_parser import cell_parser_cz_lambda
    return _lambda_table(count_list, cell_parser_cz_lambda)


def mc_file_summary(input=None, output='mc_summary.csv.gz', output_dir=None,
                    mc_format='auto', config_path=None, reference_cz=None,
                    mc_stat_feature='CHN CGN CCC', mc_stat_alias='mCH mCG mCCC',
                    num_upstr_bases=0, lambda_chrom='chrL', use_count_csv=True, cpu=1):
    """
    MappingSummary-style mC table for a set of single-cell ALLC / .cz files.

    One row per cell with ``{alias}mC/Cov/Frac`` (default mCH, mCG, mCCC),
    ``GenomeCov`` and the lambda spike-in stats
    ``Lambda{CA,CC,CT,CH,CY,CG}{mC,Cov,Frac}`` (incl. ``LambdaCYFrac``,
    ``LambdaCYCov``). Counts come from ``<file>.count.csv`` when present,
    otherwise they are recomputed from the file (via ``cytozip.mc_summary``).

    Parameters
    ----------
    input : list or str, optional
        ``*.allc.tsv.gz`` / ``*.cz`` paths: a list, directory, glob,
        comma-separated string, or a text file listing one path per line.
    output : str or None
        Output csv path (``.gz`` compresses); None to only return the table.
    output_dir : str, optional
        A yap mapping output directory; used when ``input`` is None to pick
        ``allc/*.allc.tsv.gz`` or ``cz/*.cz`` according to ``mc_format``.
    mc_format : {'auto', 'allc', 'cz'}
        File type to collect from ``output_dir``; 'auto' prefers cz/ if present.
    config_path : str, optional
        Mapping config .ini; when given, ``mc_stat_feature``,
        ``mc_stat_alias``, ``num_upstr_bases`` and ``reference_cz`` are read
        from it (explicit ``reference_cz`` still wins).
    reference_cz : str, optional
        Reference .cz, needed for .cz files without ``.count.csv``.
    cpu : int
        Number of parallel processes.

    Returns
    -------
    pandas.DataFrame indexed by cell_id.
    """
    from cytozip import mc_summary
    from ...utilities import get_configuration

    if config_path is not None:
        config = get_configuration(config_path)
        mc_stat_feature = config.get('mc_stat_feature', mc_stat_feature)
        mc_stat_alias = config.get('mc_stat_alias', mc_stat_alias)
        num_upstr_bases = int(config.get('num_upstr_bases', num_upstr_bases))
        if reference_cz is None and config.get('reference_cz', '') not in ('', 'None'):
            reference_cz = config['reference_cz']

    if input is None:
        if output_dir is None:
            raise ValueError('provide either input or output_dir')
        output_dir = pathlib.Path(output_dir).expanduser().absolute()
        cz_files = sorted(output_dir.glob('cz/*.cz'))
        allc_files = sorted(output_dir.glob('allc/*.allc.tsv.gz'))
        fmt = str(mc_format).lower()
        if fmt == 'auto':
            fmt = 'cz' if cz_files else 'allc'
        input = [str(p) for p in (cz_files if fmt == 'cz' else allc_files)]
        if not input:
            raise FileNotFoundError(f'no {fmt} files found under {output_dir}')

    return mc_summary(input=input, reference=reference_cz, output=output,
                      mc_stat_feature=mc_stat_feature, mc_stat_alias=mc_stat_alias,
                      num_upstr_bases=num_upstr_bases, lambda_chrom=lambda_chrom,
                      use_count_csv=use_count_csv, jobs=cpu)
