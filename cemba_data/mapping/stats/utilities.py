import os
import pathlib
import numpy as np
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

    count_tables = {cell_id: pd.read_csv(path, index_col=0)
                    for cell_id, path in allc_stats_dict.items()}
    final_df = _mc_pattern_table(count_tables, mc_stat_feature, mc_stat_alias)

    # add lambda DNA stats (Lambda{CA,CC,CT,CH,CY,CG}{mC,Cov,Frac})
    if str(mc_format).lower() == 'cz':
        lambda_frac = get_cz_lambda_frac(list(allc_stats_dict.values()), num_upstr_bases)
    else:
        lambda_frac = get_allc_lambda_frac(allc_list, num_upstr_bases)
    for col, data in lambda_frac.items():
        final_df[col] = data

    final_df.index.name = 'cell_id'
    return final_df


def _mc_pattern_table(count_tables, mc_stat_feature, mc_stat_alias):
    """{alias}mC/Cov/Frac and GenomeCov per cell from {cell_id: count.csv-layout table}."""
    pattern_translate = dict(zip(mc_stat_feature.split(), mc_stat_alias.split()))
    total_stats = []
    for cell_id, table in count_tables.items():
        table = table.copy()
        table['cell_id'] = cell_id
        total_stats.append(table)
    total_stats = pd.concat(total_stats)
    cell_genome_cov = pd.Series(total_stats.set_index('cell_id')['genome_cov'].to_dict())
    cell_records = []
    for pattern, alias in pattern_translate.items():
        contexts = parse_mc_pattern(pattern)
        pattern_stats = total_stats[total_stats.index.isin(contexts)]
        cell_level_data = pattern_stats.groupby('cell_id')[['mc', 'cov']].sum()
        cell_level_data['frac'] = cell_level_data['mc'] / cell_level_data['cov']
        cell_level_data = cell_level_data.rename(
            columns={'frac': f'{alias}Frac', 'mc': f'{alias}mC', 'cov': f'{alias}Cov'})
        cell_records.append(cell_level_data)
    final_df = pd.concat(cell_records, axis=1, sort=True).reindex(list(count_tables))
    final_df['GenomeCov'] = cell_genome_cov
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


def _resolve_allc_inputs(input):
    """Expand an ALLC path / list / directory / glob / comma-separated string / path-list file into ALLC paths."""
    import glob
    if isinstance(input, (list, tuple)):
        paths = [str(p) for p in input]
    else:
        s = os.path.expanduser(str(input))
        if os.path.isdir(s):
            paths = sorted(glob.glob(os.path.join(s, '*.allc.tsv.gz')))
        elif any(ch in s for ch in '*?['):
            paths = sorted(glob.glob(s))
        elif ',' in s:
            paths = [p for p in s.split(',') if p]
        elif s.endswith(('.gz', '.bgz')):
            paths = [s]
        else:
            # text file listing one path per line (first tab-separated column)
            with open(s) as f:
                paths = [line.split('\t')[0].strip() for line in f if line.strip()]
    if not paths:
        raise ValueError(f'no ALLC files found from {input!r}')
    paths = [os.path.abspath(os.path.expanduser(p)) for p in paths]
    not_allc = [p for p in paths if not p.endswith(('.gz', '.bgz'))]
    if not_allc:
        raise ValueError(f'mc_file_summary only accepts bgzipped ALLC files, got {not_allc[:3]}')
    return paths


def _allc_cell_stats(args):
    """(cell_id, count table, lambda record) of one ALLC; writes ``<allc>.count.csv`` when missing."""
    from ...hisat3n.stats_parser import cell_parser_allc_lambda, _lambda_records
    path, overwrite, num_upstr_bases = args
    cell_id = pathlib.Path(path).name.split('.')[0]
    count_path = f'{path}.count.csv'
    if not overwrite and os.path.exists(count_path):
        table = pd.read_csv(count_path, index_col=0)
        if 'genome_cov' not in table.columns:
            table['genome_cov'] = np.nan
    else:
        parts = []
        for chunk in pd.read_csv(path, sep='\t', header=None, usecols=[3, 4, 5],
                                 dtype={3: str}, chunksize=5_000_000):
            chunk.columns = ['context', 'mc', 'cov']
            parts.append(chunk.groupby('context')[['mc', 'cov']].sum())
        table = (pd.concat(parts).groupby(level=0).sum().astype('int64') if parts
                 else pd.DataFrame(columns=['mc', 'cov'], dtype='int64'))
        table.index.name = None
        table['mc_rate'] = table['mc'] / table['cov']
        # genome_cov needs every covered position (not in an ALLC): keep the old value if any
        genome_cov = np.nan
        if os.path.exists(count_path):
            old = pd.read_csv(count_path, index_col=0)
            if 'genome_cov' in old.columns and len(old):
                genome_cov = old['genome_cov'].iloc[0]
        table['genome_cov'] = genome_cov
        try:
            table.to_csv(count_path)
        except OSError as e:
            import warnings
            warnings.warn(f'could not write {count_path}: {e}')
    lam = cell_parser_allc_lambda(path, c_pos=num_upstr_bases)
    if lam.empty:  # no chrL in this cell
        zeros = dict.fromkeys(['CA', 'CC', 'CT', 'CG'], 0)
        lam = _lambda_records(zeros, dict(zeros), cell_id)
    return cell_id, table, lam


def mc_file_summary(input=None, output='mc_summary.csv.gz', output_dir=None, config_path=None,
                    mc_stat_feature='CHN CGN CCC', mc_stat_alias='mCH mCG mCCC',
                    num_upstr_bases=None, overwrite=False, cpu=1):
    """
    MappingSummary-style mC table for a set of single-cell ALLC files.

    One row per cell with ``{alias}mC/Cov/Frac`` (default mCH, mCG, mCCC),
    ``GenomeCov`` and the lambda (chrL) spike-in stats
    ``Lambda{CA,CC,CT,CH,CY,CG}{mC,Cov,Frac}`` (incl. ``LambdaCYFrac``,
    ``LambdaCYCov``), computed exactly as in the pipeline MappingSummary.
    For cytozip ``.cz`` files use ``cytozip.mc_summary`` instead.

    Per-cell mC counts come from ``{cell_id}.allc.tsv.gz.count.csv`` next to
    each ALLC (written by ``allcools bam-to-allc --save_count_df`` in the
    pipeline). When it is missing, or ``overwrite=True``, it is generated by
    summing the ALLC per context and saved in the same layout (``mc``,
    ``cov``, ``mc_rate``, ``genome_cov``). Lambda stats are always read from
    chrL via the ``.tbi`` index.

    Parameters
    ----------
    input : list or str, optional
        tabix-indexed ``*.allc.tsv.gz`` files: a single path, a list,
        a comma-separated string, a glob, a directory (all ``*.allc.tsv.gz``
        in it), or a text file listing one path per line.
    output : str or None
        Output csv path (``.gz`` compresses); None to only return the table.
    output_dir : str, optional
        A yap mapping output directory; used when ``input`` is None to take
        ``allc/*.allc.tsv.gz``.
    config_path : str, optional
        Mapping config .ini; ``mc_stat_feature``, ``mc_stat_alias`` and
        ``num_upstr_bases`` are read from it.
    num_upstr_bases : int, optional
        Bases upstream of the C in the context (1 for NOMe), used to pick the
        lambda dinucleotide. None infers it from each ``.count.csv``.
    overwrite : bool
        Regenerate every ``.count.csv`` from its ALLC even if it exists. The
        old ``genome_cov`` is kept, since it cannot be derived from an ALLC;
        generated tables otherwise have ``genome_cov`` (and ``GenomeCov``) NaN.
    cpu : int
        Number of parallel processes (one file each).

    Returns
    -------
    pandas.DataFrame indexed by cell_id.
    """
    from concurrent.futures import ProcessPoolExecutor
    from ...utilities import get_configuration

    if config_path is not None:
        config = get_configuration(config_path)
        mc_stat_feature = config.get('mc_stat_feature', mc_stat_feature)
        mc_stat_alias = config.get('mc_stat_alias', mc_stat_alias)
        if 'num_upstr_bases' in config:
            num_upstr_bases = int(config['num_upstr_bases'])
    if len(mc_stat_feature.split()) != len(mc_stat_alias.split()):
        raise ValueError('mc_stat_feature and mc_stat_alias must have the same length')

    if input is None:
        if output_dir is None:
            raise ValueError('provide either input or output_dir')
        input = str(pathlib.Path(output_dir).expanduser().absolute() / 'allc' / '*.allc.tsv.gz')
    paths = _resolve_allc_inputs(input)
    cell_ids = [pathlib.Path(p).name.split('.')[0] for p in paths]
    from collections import Counter
    dup = sorted(c for c, n in Counter(cell_ids).items() if n > 1)
    if dup:
        raise ValueError(f'duplicated cell ids in input: {dup[:5]}')

    tasks = [(p, overwrite, num_upstr_bases) for p in paths]
    if cpu > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(min(cpu, len(tasks))) as executor:
            results = list(executor.map(_allc_cell_stats, tasks))
    else:
        results = [_allc_cell_stats(t) for t in tasks]

    final_df = _mc_pattern_table({cell: table for cell, table, _ in results},
                                 mc_stat_feature, mc_stat_alias)
    lambda_df = pd.DataFrame({cell: lam for cell, _, lam in results}).T
    for col, data in lambda_df.items():
        final_df[col] = data
    final_df.index.name = 'cell_id'
    if output is not None:
        final_df.to_csv(os.path.abspath(os.path.expanduser(output)))
    return final_df
