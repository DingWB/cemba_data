from subprocess import run

from collections import defaultdict
from typing import Union, Tuple, Callable
from ._open import open_allc
from .utilities import tabix_allc
from ..utilities import parse_mc_pattern


def _merge_cg_strand(in_path, out_path):
    """
    Merge Strand used after extract ALLC step in extract_allc, so no need to check context.
    """
    prev_line = None
    cur_chrom = None

    with open_allc(in_path) as allc, \
            open_allc(out_path, 'w') as out_allc:
        for line in allc:
            cur_line = line.strip('\n').split('\t')
            if cur_line[0] != cur_chrom:
                if prev_line is not None:
                    out_allc.write('\t'.join(prev_line) + '\n')
                prev_line = cur_line
                cur_chrom = cur_line[0]
                continue
            if prev_line is None:
                prev_line = cur_line
                continue
            else:
                # pos should be continuous, strand should be reverse
                if int(prev_line[1]) + 1 == int(cur_line[1]) and prev_line[2] != cur_line[2]:
                    new_line = prev_line[:4] + [str(int(prev_line[4]) + int(cur_line[4])),
                                                str(int(prev_line[5]) + int(cur_line[5])), '1']
                    out_allc.write('\t'.join(new_line) + '\n')
                    prev_line = None
                # otherwise, only write and update prev_line
                else:
                    out_allc.write('\t'.join(prev_line) + '\n')
                    prev_line = cur_line
    return None


def _check_strandness_parameter(strandness) -> str:
    strandness = str(strandness).lower()
    if strandness in {'both', 'b'}:
        return 'Both'
    elif strandness in {'merge', 'm'}:
        # first getting both, deal with strand merge later
        return 'MergeTmp'
    elif strandness in {'split', 's'}:
        return 'Split'
    else:
        raise ValueError(f'Unknown value for strandness: {strandness}')


def _check_out_format_parameter(out_format) -> Tuple[str, Callable[[list], str]]:
    def _extract_allc_format(allc_line_list):
        # keep allc format
        return '\t'.join(allc_line_list)

    def _extract_bedgraph_cov_format(allc_line_list):
        # only chrom, pos, pos, cov
        allc_line_list = [allc_line_list[i] for i in [0, 2, 2, 5]]
        return '\t'.join(allc_line_list)

    def _extract_bedgraph_rate_format(allc_line_list):
        # only chrom, pos, pos, mc/cov
        allc_line_list = [allc_line_list[i] for i in [0, 2, 2]] + \
                         [f'{(int(allc_line_list[4]) / int(allc_line_list[5])):.3f}']
        return '\t'.join(allc_line_list)

    def _extract_bed5_format(allc_line_list):
        # only chrom, pos, pos, mc, cov
        allc_line_list = [allc_line_list[i] for i in [0, 2, 2, 4, 5]]
        return '\t'.join(allc_line_list)

    out_format = str(out_format).lower()
    if out_format == 'allc':
        return 'allc.tsv.gz', _extract_allc_format
    elif out_format == 'bed5':
        return 'bed5.bed.gz', _extract_bed5_format
    elif out_format == 'bg-cov':
        return 'cov.bg.gz', _extract_bedgraph_cov_format
    elif out_format == 'bg-rate':
        return 'rate.bg.gz', _extract_bedgraph_rate_format
    else:
        raise ValueError(f'Unknown value for out_format: {out_format}')


def extract_allc(allc_path: str,
                 out_prefix: str,
                 mc_contexts: Union[str, list] = 'CGN',
                 strandness: str = 'both',
                 output_format: str = 'allc',
                 region: str = None) -> None:
    """
    Extract ALLC file by context and/or strand. Save to several different format.

    Parameters
    ----------
    allc_path
        Path of the input ALLC file.
    out_prefix
        Path prefix of the output ALLC file.
    mc_contexts
        Space separated mC contexts to extract.
    strandness
        What to do with strand information, possible values are:
        1. both: save +/- strand together in one file;
        2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-);
        3. merge: This will only merge the count on adjacent CpG in +/- strands, only work for CpG like context.
                  For non-CG context, its the same as both.
    output_format
        Output format of extracted information, possible values are:
        1. allc: keep the allc format
        2. bed5: 5-column bed format, chrom, pos, pos, mc, cov
        3. bg-cov: bedgraph format, chrom, pos, pos, cov
        4. bg-rate: bedgraph format, chrom, pos, pos, mc/cov
    region
        Only extract records from certain genome region(s) via tabix, multiple region can be provided in tabix form.

    Returns
    -------

    """
    # TODO add parallel to this
    # prepare parms
    out_prefix = out_prefix.rstrip('.')
    if isinstance(mc_contexts, str):
        mc_contexts = mc_contexts.split(' ')
    mc_contexts = list(set(mc_contexts))
    strandness = _check_strandness_parameter(strandness)
    out_suffix, line_func = _check_out_format_parameter(output_format)

    # because mc_contexts can overlap (e.g. CHN, CAN)
    # each context may associate to multiple handle
    context_handle = defaultdict(list)
    handle_collect = []
    for mc_context in mc_contexts:
        parsed_context_set = parse_mc_pattern(mc_context)
        if strandness == 'Split':
            w_handle = open_allc(out_prefix + f'.{mc_context}-Watson.{out_suffix}', 'w')
            handle_collect.append(w_handle)
            c_handle = open_allc(out_prefix + f'.{mc_context}-Crick.{out_suffix}', 'w')
            handle_collect.append(c_handle)
            for mc_pattern in parsed_context_set:
                # handle for Watson/+ strand
                context_handle[(mc_pattern, '+')].append(w_handle)
                # handle for Crick/- strand
                context_handle[(mc_pattern, '-')].append(c_handle)
        else:
            # handle for both strand
            _handle = open_allc(out_prefix + f'.{mc_context}-{strandness}.{out_suffix}', 'w')
            handle_collect.append(_handle)
            for mc_pattern in parsed_context_set:
                context_handle[mc_pattern].append(_handle)

    # split file first
    # strandness function
    with open_allc(allc_path, region=region) as allc:
        if strandness == 'Split':
            for line in allc:
                cur_line = line.split('\t')
                try:
                    # key is (context, strand)
                    [h.write(line_func(cur_line)) for h in context_handle[(cur_line[3], cur_line[2])]]
                except KeyError:
                    continue
        else:
            for line in allc:
                cur_line = line.split('\t')
                try:
                    # key is context
                    [h.write(line_func(cur_line)) for h in context_handle[cur_line[3]]]
                except KeyError:
                    continue
    for handle in handle_collect:
        handle.close()

    if 'allc' in out_suffix:
        for mc_context in mc_contexts:
            # tabix ALLC file
            if strandness == 'Split':
                in_path = out_prefix + f'.{mc_context}-Watson.{out_suffix}'
                tabix_allc(in_path)
                in_path = out_prefix + f'.{mc_context}-Crick.{out_suffix}'
                tabix_allc(in_path)
            elif strandness == 'MergeTmp':
                in_path = out_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
                if 'CG' in mc_context:
                    out_path = out_prefix + f'.{mc_context}-Merge.{out_suffix}'
                    _merge_cg_strand(in_path, out_path)
                    run(['rm', '-f', in_path], check=True)
                else:
                    # for non-CG, there is no need to merge strand
                    out_path = out_prefix + f'.{mc_context}-Both.{out_suffix}'
                    run(['mv', in_path, out_path], check=True)
                tabix_allc(out_path)
            else:
                in_path = out_prefix + f'.{mc_context}-{strandness}.{out_suffix}'
                tabix_allc(in_path)
    return
