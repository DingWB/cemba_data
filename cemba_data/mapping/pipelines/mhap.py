import os, sys
import pandas as pd
import fire
import numpy as np
import gzip
import pysam
from collections import defaultdict

def determine_method(bam_path):
	fbam = pysam.AlignmentFile(os.path.expanduser(bam_path), 'rb')
	read = fbam.__next__()
	fbam.close()
	# Using the YZ tag of hisat-3n BAM file or XG tag of bismark BAM file.
	# If the YZ tag is "+" or XG tag is "CT", the read is C to T conversion, I change the flag to forward mapping
	# no matter R1 or R2 by read.is_forward = True
	# If the YZ tag is "-" or XG tag is "GA", the read is G to A conversion, I change the flag to reverse mapping
	# no matter R1 or R2 by read.is_forward = False
	# In this case, the read orientation is the same as bismark bam file, and the following base
	if read.has_tag('YZ') and read.get_tag('YZ') in ['+', '-']:  # hisat-3n
		return 'hisat-3n'
	elif read.has_tag('XG') and read.get_tag('XG') in ['CT', 'GA']:  # bismark
		return 'bismark'
	else:
		return 'directional'
def _is_hisat3n_ct_read(read):
	return read.get_tag('YZ') == '+'
def _is_bismark_ct_read(read):
	return read.get_tag('XG') == 'CT'
def _is_forward_read(read):
	return read.is_forward

def _bismark_is_read1(read):
	read1=True if read.get_tag('XR')=='GA' else False # R1: G/A converted, R2: C/T converted
	# for strand, if read.get_tag('XG')=='GA', then, it must be reverse strand, elsewise, it is forward strand (mapped to C/T reference)
	return read1

def _hisat3n_is_read1(read):
	return read.is_read1

def bam2mhap(bam_path=None,annotation="~/Ref/hg38/annotations/hg38_allc.gz",
		 output="test.mhap",method=None,pattern="CGN"):
	"""
	convert bam into .mhap.gz, similar to:
		prefix=""
		java -jar ${mHapSuite} convert --inputFile ${prefix}.bam --cpgPath ${cpgPath} --outPutFile ${prefix}.mhap.gz
		zcat ${prefix}.mhap.gz | sort -k1,1 -k2,2n | bgzip > ${prefix}.mhap.gz.tmp.gz
		rm -f ${prefix}.mhap.gz && mv ${prefix}.mhap.gz.tmp.gz ${prefix}.mhap.gz
		tabix -b 2 -e 3 ${prefix}.mhap.gz

	Parameters
	----------
	bam_path :
	annotation :
		path to *_allc annotation path.
	output :

	Returns
	-------

	"""
	def write_mhap(R,output):
		df = pd.DataFrame(R, columns=['chrom', 'start', 'end', 'vector', 'count', 'strand'])
		df = df.groupby(['chrom', 'start', 'end', 'vector', 'strand'], as_index=False).sum()
		df.sort_values(['chrom', 'start', 'end', 'count'], ascending=[True, True, True, False], inplace=True)
		df = df.loc[:, ['chrom', 'start', 'end', 'vector', 'count', 'strand']]
		if not os.path.exists(output):
			df.to_csv(output, sep='\t', header=False, index=False)
		else:
			df.to_csv(output, sep='\t', header=False, index=False, mode='a')

	def get_cytosine_positions(f,chrom,patterns):
		forward_cytosine_pos=[]
		reverse_cytosine_pos = []
		for line in f.fetch(chrom):
			values=line.strip().split('\t') #chrom,start,end,context, strand
			pattern=values[3]
			if pattern not in patterns:
				continue
			pos=int(values[1]) # start and end are 1-based
			strand=values[-1]
			if strand=='+':
				forward_cytosine_pos.append(pos)
			else:
				reverse_cytosine_pos.append(pos)
		return set(forward_cytosine_pos),set(reverse_cytosine_pos)

	output=os.path.expanduser(output)
	if isinstance(pattern,str):
		if pattern=='CGN':
			patterns=['CGA','CGC','CGT','CGG']
		elif pattern=='CHN':
			patterns = ['CAA', 'CCA', 'CTA'] + [
				'CAC', 'CCC', 'CTC'] + [
				'CAT', 'CCT', 'CTT'] + [
				'CAG', 'CCG', 'CTG']
		else:
			raise ValueError("Currently only support CGN and CHN")
	else:
		assert isinstance(pattern,(list,tuple))
		patterns=pattern
	# print(patterns)
	patterns=set(patterns)

	for file in [output,f"{output}.gz",f"{output}.gz.tbi"]:
		if os.path.exists(file):
			print(f"Deleting existed file: {file}")
			os.remove(file)
	if method is None:
		try:
			method=determine_method(bam_path)
		except:
			method="directional" # bam file is empty, will touch an empty output file finally
	if method.lower() in ['hisat3n','hisat-3n','hisat3']:
		ct_read_func=_is_hisat3n_ct_read
		is_read1_func=_hisat3n_is_read1
	elif method.lower()=='bismark':
		ct_read_func = _is_bismark_ct_read
		is_read1_func =_bismark_is_read1
	else:
		ct_read_func = _is_forward_read
		is_read1_func = _hisat3n_is_read1
	fbam = pysam.AlignmentFile(os.path.expanduser(bam_path), 'rb')
	f = pysam.TabixFile(os.path.expanduser(annotation)) # for allc.bed.gz, 1-based
	for chrom in f.contigs:
		forward_cytosine_pos,reverse_cytosine_pos=get_cytosine_positions(f,chrom,patterns) # 1-based
		if len(forward_cytosine_pos) == 0 and len(reverse_cytosine_pos)==0:  # cytosine_positions is 1-based, pos of C
			continue
		print(chrom)
		R = []
		pre_pos,pre_ct_read,pre_read1=None,None,None
		for read in fbam.fetch(chrom):
			if 'D' in read.cigarstring or 'I' in read.cigarstring or read.cigarstring.count('M') != 1:
				continue # skip the reads with Delete, Insertion or other indel or variation
			if sum([cigar[1] for cigar in read.cigar]) != len(read.seq):
				continue
			if sum([cigar[1] for cigar in read.cigar if cigar[0]==0]) != len(read.positions):
				continue # cigar: 0 means Match, if read.is_secondary: #not read.is_proper_pair:
			# read.pos is 0-based, read.aend = read.pos+1 + cigar length - 1 # 1-based
			ct_read = ct_read_func(read)
			is_read1=is_read1_func(read)
			if (read.pos,ct_read,is_read1) == (pre_pos,pre_ct_read,pre_read1):
				continue #duplicates
			seq=read.seq
			for cigar,length in read.cigar:
				if cigar!=0: # not equal to M, not matched
					seq=seq[length:] # seq and sum of cigar[1] are equal
					continue #start pos in bam file is the start pos of the first matched 'M' base, skipping unmatched base
				else:
					break #get the first matched cigar and length
			# assert cigar==0
			positions = read.positions # 0-based
			seq=seq[:length]
			# assert len(seq)==len(read.positions)
			if ct_read: #CT read (OT or CTOT, mapped to C/T ref), posision of C
				ovlp_cytosine_idx = [idx for idx, pos in enumerate(positions) if
									 pos + 1 in forward_cytosine_pos] #pos is 0-based, forward_cytosine_pos is 1-based
				# overlapped_cytosine_idx is the index for C
				# cytosine=''.join([seq[idx] for idx in overlapped_cytosine_idx])
				vector = ''.join(['1' if seq[idx]=='C' else '0' for idx in ovlp_cytosine_idx]) #1 for methylated, 0 for unmethylated
				# read.is_forward = True  # to be consistent with bismark: if read.get_tag('XG')=='GA', then, it must be reverse strand, elsewise, it is forward strand (mapped to C/T reference)
				strand='+'
			else: # G/A read (OB or CTOB, Mapped to G/A ref), postion of G (next of C)
				ovlp_cytosine_idx = [idx for idx, pos in enumerate(positions) if
									 pos + 1 in reverse_cytosine_pos] #index for the pos of G (+) or C (-)
				vector = ''.join(['1' if seq[idx]=='G' else '0' for idx in ovlp_cytosine_idx]) # G means the C is in the reverse strand
				# read.is_forward = False  # to be consistent with bismark
				strand = '-'
			if len(ovlp_cytosine_idx)>0:
				R.append([read.reference_name, positions[ovlp_cytosine_idx[0]]+1,
						  positions[ovlp_cytosine_idx[-1]]+1, vector, 1, strand])
			pre_pos,pre_ct_read,pre_read1=read.pos,ct_read,is_read1
		if len(R) > 0:
			write_mhap(R, output)
	if not os.path.exists(output):
		os.system(f"touch {output}")
	os.system(f"bgzip {output} && tabix -b 2 -e 3 {output}.gz")
	fbam.close()
	f.close()

if __name__ == "__main__":
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(serialize=lambda x:print(x) if not x is None else print(""))