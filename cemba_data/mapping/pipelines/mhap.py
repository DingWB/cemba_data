import os, sys
import pandas as pd
import fire
import pysam

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

def bam2mhap(bam_path=None,cpg_path="~/Ref/hg38/annotations/hg38_CpG.gz",
		 output="test.mhap",method=None):
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
	cpg_path :
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

	output=os.path.expanduser(output)
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
	cpg = pysam.TabixFile(os.path.expanduser(cpg_path)) # positions is 1-based
	for chrom in cpg.contigs:
		cytosine_positions = set(
			[int(line.split('\t')[1]) for line in cpg.fetch(chrom)])  # 1-based
		if len(cytosine_positions) == 0:  # cytosine_positions is 1-based, pos of C
			continue
		print(chrom)
		R = [] #positions=[int(line.split('\t')[1]) for line in cpg.fetch(chrom,read.pos+1,read.aend)]
		pre_pos,pre_ct_read,pre_read1=None,None,None
		for read in fbam.fetch(chrom):
			if 'D' in read.cigarstring or 'I' in read.cigarstring or read.cigarstring.count('M') != 1:
				continue
			if sum([cigar[1] for cigar in read.cigar]) != len(read.seq):
				continue
			if sum([cigar[1] for cigar in read.cigar if cigar[0]==0]) != len(read.positions):
				continue # if read.is_secondary: #not read.is_proper_pair:
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
			if ct_read: # C/T read, C could be the last base pair
				positions = read.positions # positions=read.positions[i:i+cigar[1]] #0-based
				read.is_forward=True # to be consistent with bismark: if read.get_tag('XG')=='GA', then, it must be reverse strand, elsewise, it is forward strand (mapped to C/T reference)
			else: # G/A read, C must be not in the last base pair
				# positions = read.positions[i:i + cigar[1]][:-1] #exclude the potential C in the last base pair
				positions = read.positions[:-1]
				read.is_forward = False #to be consistent with bismark
			overlapped_cytosine_idx = [idx for idx,pos in enumerate(positions) if pos+1 in cytosine_positions] #set(positions).intersection(set(cytosine_positions))
			# overlapped_cytosine_idx is the index for C
			if len(overlapped_cytosine_idx) > 0:
				seq=seq[:length]
				# assert len(seq)==len(read.positions)
				if ct_read: #CT read, posision of C
					# cytosine=''.join([seq[idx] for idx in overlapped_cytosine_idx])
					vector = ''.join(['1' if seq[idx]=='C' else '0' for idx in overlapped_cytosine_idx]) #1 for methylated, 0 for unmethylated
				else: # G/A read, postion of G (next of C)
					# cytosine=''.join([seq[idx+1] for idx in overlapped_cytosine_idx])
					vector = ''.join(['1' if seq[idx+1]=='G' else '0' for idx in overlapped_cytosine_idx])
				strand='+' if read.is_forward else '-'
				R.append([read.reference_name, positions[overlapped_cytosine_idx[0]]+1,
						  positions[overlapped_cytosine_idx[-1]]+1, vector, 1, strand])
			pre_pos,pre_ct_read,pre_read1=read.pos,ct_read,is_read1
		if len(R) > 0:
			write_mhap(R, output)
	if not os.path.exists(output):
		os.system(f"touch {output}")
	os.system(f"bgzip {output} && tabix -b 2 -e 3 {output}.gz")
	fbam.close()
	cpg.close()

if __name__ == "__main__":
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(serialize=lambda x:print(x) if not x is None else print(""))