"""
Convert bisulfite alignments (hisat-3n, Bismark or other directional aligners)
into mHap files (chrom, start, end, haplotype, count, strand).

Conventions, the same as mHapTools / mHapSuite:

- start / end are the 1-based positions of the first / last site of a record,
  and the haplotype has one 0/1 per consecutive site of the site list.
- CpG outputs (every context starts with CG): sites are the top-strand C of each
  CpG for both strands; bottom-strand reads are called on the G (C + 1).
- Other contexts (e.g. CHN) are strand specific: a bottom-strand site is the
  bottom-strand C (the top-strand G), with its own site list.
"""
import os
import sys
from collections import Counter, OrderedDict

import numpy as np
import pysam

_BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate('ACGT'):
	_BASE_LUT[ord(_b)] = _i
	_BASE_LUT[ord(_b.lower())] = _i
_COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)  # A<->T, C<->G, N->N
_ALIGNED_OPS = {0, 7, 8}  # M, =, X
_QUERY_ONLY_OPS = {1, 4}  # I, S
_REF_ONLY_OPS = {2, 3}  # D, N
_MAX_OUTPUTS = 4  # two flag bits (top / bottom) per output in a uint8


def _context_code(context):
	return int(_BASE_LUT[ord(context[0])]) * 25 + int(_BASE_LUT[ord(context[1])]) * 5 + int(
		_BASE_LUT[ord(context[2])])


def _parse_contexts(pattern):
	"""IUPAC pattern(s) such as 'CGN', 'CHN', 'CGA,CGC' or a list -> set of 3-base contexts."""
	from cemba_data.utilities import parse_mc_pattern
	tokens = pattern.replace(',', ' ').split() if isinstance(pattern, str) else [str(p) for p in pattern]
	contexts = set()
	for token in tokens:
		token = token.strip().upper()
		if len(token) != 3 or token[0] != 'C':
			raise ValueError(f'mhap contexts must have 3 bases and start with C (e.g. CGN, CHN, CGA), got {token!r}')
		contexts |= {c for c in parse_mc_pattern(token) if set(c) <= set('ACGT')}
	if not contexts:
		raise ValueError(f'no valid context in pattern {pattern!r}')
	return frozenset(contexts)


class _Spec:
	"""One output: its contexts, flag bits and per-chromosome record counter."""

	def __init__(self, k, name, contexts, output):
		self.name = name
		self.contexts = contexts
		# a 3rd base that is N / IUPAC / outside the contig matches when every 3rd base of the
		# prefix is wanted (e.g. CGN, CHN), so a CpG before an N or at a contig end is kept
		self.lut_contexts = set(contexts) | {
			p + 'N' for p in {c[:2] for c in contexts} if all(p + b in contexts for b in 'ACGT')}
		self.cpg = all(c.startswith('CG') for c in contexts)
		self.top_bit = 1 << (2 * k)
		# CpG: both strands use the top-strand C list, bottom reads are read at C + 1
		self.bottom_bit = self.top_bit if self.cpg else 1 << (2 * k + 1)
		self.bottom_offset = 1 if self.cpg else 0
		# expected partner base on the read (G after a top C, C before a bottom G): CpG must keep it,
		# other contexts without CG must not show it; None when the contexts mix both
		second = {c[1] for c in contexts}
		self.partner = 'eq' if second == {'G'} else ('ne' if 'G' not in second else None)
		self.output = output
		self.plain = output[:-3]
		self.counter = Counter()

	def site(self, bottom):
		return (self.bottom_bit, self.bottom_offset) if bottom else (self.top_bit, 0)


def _context_luts(specs):
	top = np.zeros(125, dtype=np.uint8)
	bottom = np.zeros(125, dtype=np.uint8)
	for spec in specs:
		for context in spec.lut_contexts:
			code = _context_code(context)
			top[code] |= spec.top_bit
			if not spec.cpg:
				bottom[code] |= spec.bottom_bit
	return top, bottom


class _Fasta:
	"""Reference access; a plain FASTA is read directly through its .fai offsets (faster than pysam fetch)."""

	def __init__(self, path):
		self.path = os.path.expanduser(str(path))
		self._fa = pysam.FastaFile(self.path)  # also builds the .fai if missing
		self.references = set(self._fa.references)
		self._fai = None
		if not (self.path.endswith('.gz') or os.path.exists(self.path + '.gzi')):
			with open(self.path + '.fai') as f:
				self._fai = {cols[0]: tuple(map(int, cols[1:5])) for cols in (line.split('\t') for line in f)}

	def length(self, chrom):
		return self._fa.get_reference_length(chrom)

	def codes(self, chrom, start, end):
		"""Base codes (A0 C1 G2 T3 other 4) of [start, end)."""
		if self._fai is None:
			return _BASE_LUT[np.frombuffer(self._fa.fetch(chrom, start, end).encode('ascii'), dtype=np.uint8)]
		_, offset, line_bases, line_width = self._fai[chrom]
		b0 = offset + (start // line_bases) * line_width + start % line_bases
		b1 = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases + 1
		raw = np.fromfile(self.path, dtype=np.uint8, count=b1 - b0, offset=b0)
		raw = raw[raw > 32]  # drop line breaks (\n, \r)
		if len(raw) != end - start:
			raise ValueError(f'unexpected FASTA line layout of {chrom} in {self.path}')
		return _BASE_LUT[raw]

	def close(self):
		self._fa.close()


def _flags_from_fasta(fasta, chrom, luts, chunk_size=10_000_000):
	"""uint8 site flags of one chromosome; bit 2k: top-strand site of output k, bit 2k+1: bottom-strand site."""
	top_lut, bottom_lut = luts
	use_bottom = bool(bottom_lut.any())
	length = fasta.length(chrom)
	flags = np.zeros(length, dtype=np.uint8)
	for start in range(0, length, chunk_size):
		end = min(start + chunk_size, length)
		lo, hi = max(start - 2, 0), min(end + 2, length)
		codes = np.full(end - start + 4, 4, dtype=np.uint8)  # bases at start-2 .. end+1, N outside the chromosome
		offset = lo - (start - 2)
		codes[offset:offset + hi - lo] = fasta.codes(chrom, lo, hi)
		n = end - start
		mid = codes[2:2 + n]
		ic = np.flatnonzero(mid == 1)  # top-strand C: context C, +1, +2
		flags[start + ic] = top_lut[25 + codes[3 + ic] * 5 + codes[4 + ic]]
		if use_bottom:
			ig = np.flatnonzero(mid == 2)  # bottom-strand C (top G): context C, comp(-1), comp(-2)
			flags[start + ig] = bottom_lut[25 + _COMP[codes[1 + ig]] * 5 + _COMP[codes[ig]]]
	return flags


def _open_bam(bam_path):
	"""Open a BAM for streaming; no index is needed, so silence htslib's missing-index message."""
	verbosity = pysam.set_verbosity(0)
	try:
		return pysam.AlignmentFile(os.path.expanduser(bam_path), 'rb')
	finally:
		pysam.set_verbosity(verbosity)


def determine_method(bam_path, n_reads=1000):
	"""'hisat-3n' (YZ tag), 'bismark' (XG tag) or 'directional' (strand from the SAM flags)."""
	n_untagged = 0
	with _open_bam(bam_path) as bam:
		for read in bam.fetch(until_eof=True):
			if read.is_unmapped:
				continue
			if read.has_tag('YZ') and read.get_tag('YZ') in ('+', '-'):
				return 'hisat-3n'
			if read.has_tag('XG') and read.get_tag('XG') in ('CT', 'GA'):
				return 'bismark'
			n_untagged += 1
			if n_untagged >= n_reads:
				break
	return 'directional'


def _strand_func(method):
	"""read -> True for bottom-strand (G/A) reads, False for top-strand (C/T) reads, None if unknown."""
	method = method.lower()
	if method in ('hisat3n', 'hisat-3n', 'hisat3', 'hisat'):
		# YZ:A:+ / - is the conversion of the alignment; both mates (and m3C split parts) share it
		table = {'+': False, '-': True}
		tag = 'YZ'
	elif method == 'bismark':
		table = {'CT': False, 'GA': True}
		tag = 'XG'
	elif method == 'directional':
		# OT: R1 forward / R2 reverse, OB: R1 reverse / R2 forward
		return lambda read: read.is_reverse != read.is_read2 if read.is_paired else read.is_reverse
	else:
		raise ValueError(f"Unknown method {method!r}, choose from 'hisat-3n', 'bismark', 'directional'")

	def func(read):
		return table.get(read.get_tag(tag)) if read.has_tag(tag) else None

	return func


def _read_calls(read, bottom, flags, specs, min_bq, check_context):
	"""{site: 1/0} per output for one read, using the CIGAR to place every base on the reference."""
	seq = read.query_sequence
	qual = read.query_qualities
	blocks = []
	r, q = read.reference_start, 0
	for op, n in read.cigartuples:
		if op in _ALIGNED_OPS:
			blocks.append((r, q, n))
			r += n
			q += n
		elif op in _QUERY_ONLY_OPS:
			q += n
		elif op in _REF_ONLY_OPS:
			r += n
	if bottom:
		meth, unmeth, step, partner_base = 'G', 'A', -1, 'C'
		# a partner just before a block is aligned only as the last base of the previous block
		edge = {rs + n - 1: qs + n - 1 for rs, qs, n in blocks}
	else:
		meth, unmeth, step, partner_base = 'C', 'T', 1, 'G'
		edge = {rs: qs for rs, qs, _ in blocks}
	n_flags = len(flags)
	calls = []
	for spec in specs:
		bit, off = spec.site(bottom)
		rule = spec.partner if check_context else None
		spec_calls = {}
		for rs, qs, n in blocks:
			c0, c1 = max(rs - off, 0), min(rs + n - off, n_flags)
			if c1 <= c0:
				continue
			for h in np.flatnonzero(flags[c0:c1] & bit).tolist():
				c = c0 + h
				bp = c + off  # reference position of the called base
				qp = qs + bp - rs
				if qual is not None and qual[qp] < min_bq:
					continue
				base = seq[qp]
				if base == meth:
					value = 1
				elif base == unmeth:
					value = 0
				else:
					continue
				if rule is not None:
					pr = bp + step
					pq = qp + step if rs <= pr < rs + n else edge.get(pr)
					# partner outside the aligned read (read end, deletion): keep the call, as ALLC does
					if pq is not None and (seq[pq] == partner_base) != (rule == 'eq'):
						continue  # SNP / misalignment on the read breaks the context
				spec_calls[c] = value
		calls.append(spec_calls)
	return calls


class _Fragment:
	__slots__ = ('first_start', 'bottom', 'spans', 'calls')

	def __init__(self, first_start, bottom, n_specs):
		self.first_start = first_start
		self.bottom = bottom
		self.spans = []
		self.calls = [[] for _ in range(n_specs)]

	def add(self, start, end, calls):
		self.spans.append((start, end))
		for spec_calls, read_calls in zip(self.calls, calls):
			spec_calls.append(read_calls)


def _finalize_fragment(frag, flags, specs, stats):
	"""Merge the calls of all reads of a fragment and add its mHap records to the counters."""
	strand = '-' if frag.bottom else '+'
	for k, spec in enumerate(specs):
		merged, conflict = {}, set()
		for read_calls in frag.calls[k]:
			for c, v in read_calls.items():
				if merged.setdefault(c, v) != v:
					conflict.add(c)  # mates disagree: sequenced, but no call
		sites = sorted(c for c in merged if c not in conflict)
		if not sites:
			continue
		bit, off = spec.site(frag.bottom)
		lo_c = sites[0]
		all_sites = np.flatnonzero(flags[lo_c:sites[-1] + 1] & bit)
		cuts = []
		if len(all_sites) != len(sites):  # some sites between the first and last call have no call
			all_sites = (all_sites + lo_c).tolist()
			rank = {c: i for i, c in enumerate(all_sites)}
			drop = False
			for i in range(1, len(sites)):
				ra, rb = rank[sites[i - 1]], rank[sites[i]]
				if rb - ra == 1:
					continue
				# sequenced (inside a read, e.g. low BQ, deletion, mate conflict) vs. not sequenced (between mates)
				sequenced = any(s <= m + off < e for m in all_sites[ra + 1:rb] for s, e in frag.spans)
				if sequenced and not spec.split:
					drop = True
					break
				cuts.append(i)
			if drop:
				stats[f'dropped_fragments_{spec.name}'] += 1
				continue
		bounds = [0] + cuts + [len(sites)]
		for lo, hi in zip(bounds[:-1], bounds[1:]):
			haplotype = ''.join('1' if merged[c] else '0' for c in sites[lo:hi])
			spec.counter[(sites[lo] + 1, sites[hi - 1] + 1, strand, haplotype)] += 1


def _parse_chroms(chroms):
	if chroms is None:
		return None
	if isinstance(chroms, str):
		path = os.path.expanduser(chroms)
		if os.path.exists(path):
			with open(path) as f:
				return {line.split()[0] for line in f if line.strip() and not line.startswith('#')}
		if os.sep in chroms or chroms.endswith(('.sizes', '.txt', '.tsv')):
			raise FileNotFoundError(f'chrom size file {chroms} not found')
		return {c for c in chroms.replace(',', ' ').split() if c}
	return {str(c) for c in chroms}


def bam2mhap(bam_path=None, reference_fasta=None, output="test.mhap", method=None, pattern="CGN",
			 outputs=None, chroms=None, min_mapq=10, min_base_quality=20,
			 max_frag=1000, split=None, check_context=True):
	"""
	Convert a coordinate-sorted bisulfite BAM into bgzipped, tabix-indexed mHap file(s),
	similar to ``mhaptools convert`` (mHapTools >= 0.11).

	- reads: skip unmapped, secondary, supplementary, QC-fail and duplicate reads, and MAPQ < min_mapq;
	- bases: placed on the reference with the CIGAR (soft clips, insertions, deletions);
	  a site is called only from a C/T (top strand) or G/A (bottom strand) base with
	  base quality >= min_base_quality, and (check_context) when the read keeps the
	  context: the G of a CpG after a top-strand C, its C before a bottom-strand G
	  (no G / C there for CH); not checked when that base is outside the read;
	- fragments: alignments with the same read name and strand within max_frag bp
	  (mates, m3C split parts) are merged so overlapping bases are counted once;
	  sites where they disagree get no call;
	- records: a fragment is cut where sites were not sequenced (between mates);
	  a sequenced site without a call between two calls drops the fragment,
	  or cuts it with split (default: drop for CpG outputs as mHapTools, cut for other contexts).

	Parameters
	----------
	bam_path :
		coordinate-sorted BAM (hisat-3n, Bismark or a directional library).
	reference_fasta :
		genome fasta (with .fai, plain or bgzipped) the BAM was aligned to; the CpG / CH
		sites are taken from it.
	output :
		output path, ``.gz`` is appended if missing (``{output}.gz`` + ``.tbi``).
	method :
		'hisat-3n', 'bismark' or 'directional'; detected from the BAM tags when None.
	pattern :
		context pattern of ``output``: IUPAC like 'CGN' / 'CHN', comma separated or a list.
	outputs :
		{pattern: output} to write several outputs in one pass over the BAM,
		e.g. {'CGN': 'cell.CG.mhap.gz', 'CHN': 'cell.CH.mhap.gz'}; overrides pattern / output.
	chroms :
		chromosomes to keep: chrom size file, list or comma-separated string; None for all.
	min_mapq :
		minimum mapping quality.
	min_base_quality :
		minimum base quality of a call.
	max_frag :
		maximum distance (bp) between the starts of two alignments of the same fragment.
	split :
		True: cut instead of drop a fragment at a sequenced site without a call; False: drop;
		None: drop for CpG outputs, cut for other contexts (dense CH sites would lose many fragments).
	check_context :
		skip calls whose context is broken on the read (SNP / misalignment), as mHapTools does.

	Returns
	-------
	dict of counts (reads, filtered reads, fragments, records ...)
	"""
	if bam_path is None:
		raise ValueError('bam_path is required')
	if reference_fasta is None:
		raise ValueError('reference_fasta is required')
	if outputs is None:
		outputs = [(pattern, output)]
	elif isinstance(outputs, dict):
		outputs = list(outputs.items())
	if not 0 < len(outputs) <= _MAX_OUTPUTS:
		raise ValueError(f'1 to {_MAX_OUTPUTS} outputs are supported, got {len(outputs)}')
	specs = []
	for k, (pat, out) in enumerate(outputs):
		out = os.path.abspath(os.path.expanduser(str(out)))
		if not out.endswith('.gz'):
			out += '.gz'
		name = pat if isinstance(pat, str) else ','.join(map(str, pat))
		spec = _Spec(k, name, _parse_contexts(pat), out)
		spec.split = (not spec.cpg) if split is None else bool(split)
		specs.append(spec)
	if len({s.output for s in specs}) != len(specs):
		raise ValueError('outputs must be different files')
	for spec in specs:
		for path in [spec.plain, spec.output, f'{spec.output}.tbi', f'{spec.output}.csi']:
			if os.path.exists(path):
				os.remove(path)

	bam_path = os.path.expanduser(bam_path)
	if method is None:
		method = determine_method(bam_path)
	strand_of = _strand_func(method)
	keep_chroms = _parse_chroms(chroms)
	fasta = _Fasta(reference_fasta)
	luts = _context_luts(specs)
	stats = Counter()
	exclude = 0x4 | 0x100 | 0x200 | 0x400 | 0x800  # unmapped, secondary, QC fail, duplicate, supplementary

	def load_flags(chrom, length):
		if keep_chroms is not None and chrom not in keep_chroms:
			return None
		if chrom not in fasta.references:
			print(f'[bam2mhap] {chrom} not in the reference fasta, skipped', file=sys.stderr)
			return None
		if fasta.length(chrom) != length:
			raise ValueError(f'{chrom} has length {length} in the BAM but '
							 f'{fasta.length(chrom)} in {reference_fasta}')
		return _flags_from_fasta(fasta, chrom, luts)

	def finish_chrom(chrom, flags, pending, handles):
		while pending:
			_finalize_fragment(pending.popitem(last=False)[1], flags, specs, stats)
		for spec, handle in zip(specs, handles):
			for (start, end, strand, hap), count in sorted(spec.counter.items()):
				handle.write(f'{chrom}\t{start}\t{end}\t{hap}\t{count}\t{strand}\n')
				stats[f'records_{spec.name}'] += 1
			spec.counter.clear()

	handles = [open(spec.plain, 'w') for spec in specs]
	try:
		with _open_bam(bam_path) as bam:
			cur_tid, chrom, flags, last_start = None, None, None, -1
			done_tids = set()
			pending = OrderedDict()
			for read in bam.fetch(until_eof=True):
				stats['reads'] += 1
				tid = read.reference_id
				if tid < 0:  # unplaced unmapped reads are at the end of a sorted BAM
					break
				if tid != cur_tid:
					if flags is not None:
						finish_chrom(chrom, flags, pending, handles)
					if tid in done_tids:
						raise ValueError(f'{bam_path} is not coordinate-sorted')
					done_tids.add(tid)
					cur_tid, chrom, last_start = tid, read.reference_name, -1
					flags = load_flags(chrom, bam.get_reference_length(chrom))
				start = read.reference_start
				if start < last_start:
					raise ValueError(f'{bam_path} is not coordinate-sorted')
				last_start = start
				if flags is None:
					stats['reads_skipped_chrom'] += 1
					continue
				if (read.flag & exclude or read.mapping_quality < min_mapq or read.query_sequence is None
						or not read.cigartuples):
					stats['reads_filtered'] += 1
					continue
				bottom = strand_of(read)
				if bottom is None:
					stats['reads_no_strand'] += 1
					continue
				stats['reads_used'] += 1
				while pending:
					first = next(iter(pending.values()))
					if first.first_start >= start - max_frag:
						break
					_finalize_fragment(pending.popitem(last=False)[1], flags, specs, stats)
				key = (read.query_name, bottom)
				frag = pending.get(key)
				if frag is not None and start - frag.first_start > max_frag:
					_finalize_fragment(pending.pop(key), flags, specs, stats)
					frag = None
				if frag is None:
					frag = pending[key] = _Fragment(start, bottom, len(specs))
					stats['fragments'] += 1
				frag.add(start, read.reference_end,
						 _read_calls(read, bottom, flags, specs, min_base_quality, check_context))
			if flags is not None:
				finish_chrom(chrom, flags, pending, handles)
	finally:
		for handle in handles:
			handle.close()
		fasta.close()

	for spec in specs:
		pysam.tabix_index(spec.plain, seq_col=0, start_col=1, end_col=2, zerobased=False, force=True)
	stats = dict(stats)
	print(f'[bam2mhap] {bam_path} ({method}): {stats}', file=sys.stderr)
	if not stats.get('reads_used'):
		print(f'[bam2mhap] WARNING: no read of {bam_path} was used; check chroms, the reference and the strand tags',
			  file=sys.stderr)
	return stats


if __name__ == "__main__":
	import fire

	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(bam2mhap, serialize=lambda x: print(x) if x is not None else print(""))
