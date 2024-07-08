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

def mhl(vectors,counts,k=10):
	"""
	Calculate MHL, (methylated haplotype load): partitioned each segment
	into methylation haplotype blocks (MHBs). MHBs were defined as the
	genomic region in which the r2 value of two adjacent CpG sites is no less
	than 0.5. MHL=normalized fraction of methylated haplotypes at different lengths.
	Ref: Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tumor tissue-of-origin mapping from plasma DNA
	for example:
		vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
		vectors2=['1100']*8+['0011']*8
		vectors3=['1100','0110','0011','1110','0111']
		vectors4=['1010','0101','1001','1011','1101']
		counts=[1]*16
		print(mhl(vectors1,counts))
		print(mhl(vectors2,counts))
		print(mhl(['0000']*16,counts))
		print(mhl(['1111']*16,counts))
		print(mhl(vectors3,[1]*5))
		print(mhl(vectors4,[1]*5))
	Parameters
	----------
	vectors :

	Returns
	-------

	"""
	mhl_dict={}
	max_l= max([len(vector) for vector in vectors])
	max_k=min(max_l,k)
	for i in range(1,max_k+1):
		unreachable_n=i-1
		s,total_count=0,0
		for vector, count in zip(vectors, counts):
			n=len(vector)-unreachable_n
			if n<=0:
				continue
			s+=n * count
			total_count+=count * sum([1 if vector[j:j+i] == '1'*i else 0 for j in range(0,len(vector)-unreachable_n)])
		mhl_dict[i]=total_count / s
	return sum([i*mhl_dict[i] for i in mhl_dict]) / sum([i for i in mhl_dict])

def methylation_entropy(vectors,counts): #wrong, deprecated
	"""
	Calculate methylatyion entropy, for example:
		vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
		vectors2=['1100']*8+['0011']*8
		counts=[1]*16
		print(methylation_entropy(vectors1,counts))
		print(methylation_entropy(vectors2,counts))
		print(methylation_entropy(['0000']*16,counts))
		print(methylation_entropy(['1111']*16,counts))
	Parameters
	----------
	vectors :
	counts :

	Returns
	-------

	"""
	df=pd.DataFrame.from_dict(dict(vectors=vectors,counts=counts)).groupby('vectors',as_index=False)['counts'].sum()
	df['P']=df.counts / df.counts.sum() #P for each haplotype
	return df.P.apply(lambda x:-x*np.log2(x)).sum() / df.vectors.apply(len).max()

def entropy(vectors,counts,k=4):
	"""
	Calculate MHL, (methylated haplotype load): partitioned each segment
	into methylation haplotype blocks (MHBs). MHBs were defined as the
	genomic region in which the r2 value of two adjacent CpG sites is no less
	than 0.5. MHL=normalized fraction of methylated haplotypes at different lengths.
	Ref: Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tumor tissue-of-origin mapping from plasma DNA
	for example:
		vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
		vectors2=['1100']*8+['0011']*8
		counts=[1]*16
		print(entropy(vectors1,counts))
		print(entropy(vectors2,counts))
		print(entropy(['0000']*16,counts))
		print(entropy(['1111']*16,counts))
		print(entropy(['1111']*8+['0000']*8,counts))
	Parameters
	----------
	vectors :

	Returns
	-------

	"""
	kmer_count=defaultdict(lambda : 0)
	for vector,count in zip(vectors,counts):
		if len(vector) < k:
			continue
		for i in range(len(vector) - k + 1):
			kmer_count[vector[i: i + k]]+=count
	values=[kmer_count[k] for k in kmer_count]
	t=sum(values)
	return sum([-1*(v / t)*np.log2(v/t) for v in values]) / k

def mbs(vectors,counts,k=4):
	"""
	MBS (methylation block score): N / nReads
	ref: Ultrasensitive detection of circulating tumour DNA via deep methylation sequencing aided by machine learning
	Parameters
	----------
	vectors :
	counts :

	Returns
	-------

	"""
	r=0
	s=0
	for vector,count in zip(vectors,counts): #for each read
		if len(vector) < k:
			continue
		n=len(vector)
		s+=count
		constant=count / (n*n)
		# find the longest successive methylated CpG
		for substr in vector.split('0'):
			length=len(substr)
			r+=length*length*constant
	try:
		return r / s
	except:
		return np.nan

def buildBinaryMatrix(mhap_path,cpg_path,chrom,start,end,shift=500):
	f_cpg = pysam.TabixFile(os.path.expanduser(cpg_path))
	posArray = [int(line.split('\t')[1]) for line in f_cpg.fetch(chrom, start - shift, end + shift)] #cpg start position array
	f_cpg.close()
	Nc = len(posArray)  # No. of CpGs

	Matrix=[]
	rows=[]
	f_mhap = pysam.TabixFile(os.path.expanduser(mhap_path))
	for line in f_mhap.fetch(chrom, start=start, end=end):
		values = line.split('\t')
		start1=int(values[1]) #start cpg pos of this vector (read)
		assert start1 in posArray
		end1=int(values[2])
		vector=values[3]
		rows.append(f"{start1}-{end1}-{int(values[4])}-{values[5]}") #start-end-count-strand
		matrix_row=[]
		i=0
		for pos in posArray:
			if pos < start1 or pos > end1:
				matrix_row.append(np.nan)
				continue
			matrix_row.append(int(vector[i]))
			i+=1
		Matrix.append(matrix_row)
	f_mhap.close()
	df=pd.DataFrame(Matrix,index=rows,columns=posArray)
	# df.loc[:,df.notna().sum()>0], np.sum, axis=1,sum by rows
	# df: each row is a record (read), each column is a CpG site, index are: start-end-count-strand
	return df

def merge_mhaps(mhap_paths=[],max_open=100):
	if len(mhap_paths) <= max_open:
		f_mhaps=[pysam.TabixFile(os.path.expanduser(mhap_path)) for mhap_path in mhap_paths]
	def get_records(f_mhaps, chrom, start, end):
		vectors, counts = [], []
		for f_mhap in f_mhaps:
			for line in f_mhap.fetch(chrom, start=start, end=end):
				values = line.split('\t')
				vectors.append(values[3])
				counts.append(int(values[4]))
		return vectors, counts

def stat_single_region(f_mhap,chrom,start,end):
	vectors,counts=[],[]
	for line in f_mhap.fetch(chrom,start=start,end=end):
		values=line.split('\t')
		vectors.append(values[3])
		counts.append(int(values[4]))
	if len(vectors)==0:
		return [0,0,0,0,0,0,0,np.nan,np.nan,
				np.nan,np.nan,np.nan,np.nan,np.nan]
	nReads=0 #sum(counts)
	mBase=0 #sum([vector.count('1')*count for vector,count in zip(vectors,counts)])
	cBase=0 #sum([vector.count('0')*count if '1' in vector else 0 for vector,count in zip(vectors,counts)])
	tBase=0 #sum([len(vector)*count for vector,count in zip(vectors,counts)])
	K4plus=0 #sum([count if len(vector)>=4 else 0 for vector,count in zip(vectors,counts)])
	nDR=0 #sum([count if len(set(record.vector))==2 and len(record.vector)>=4])
	nMR=0 #sum([count if '1' in vector and len(record.vector)>=4])
	for vector, count in zip(vectors, counts):
		nReads+=count
		mBase+=vector.count('1')*count
		if '1' in vector:
			cBase+=vector.count('0')*count
		tBase+=len(vector)*count
		if len(vector) >= 4:
			K4plus+=count
			if len(set(vector)) == 2:
				nDR+=count
			if '1' in vector:
				nMR+=count
	if K4plus > 0:
		PDR,CHALM=nDR / K4plus,nMR / K4plus
	else:
		PDR,CHALM=np.nan,np.nan
	return [nReads,mBase,cBase,tBase,
				K4plus,nDR,nMR,mBase / tBase,
				CHALM,PDR,
				mhl(vectors,counts),
				mbs(vectors, counts), #MBS
				cBase / tBase, #MCR
				entropy(vectors,counts),
				]

def open1(infile):
	if hasattr(infile, 'readline'):
		return infile
	if infile.endswith('.gz'):
		f = gzip.open(infile, 'rb')
	else:
		f = open(infile, 'r')
	return f

def stat_mhap(mhap_path=None,cpg_path=None,region=None,bed=None,
		 output=None):
	"""
	stat mhap file, for example:
		stat_mhap -m HBA_211015_H1930004_CB63_V1C_3C_1_P5-6-K9-G24.mC.mhap.gz -c ~/Ref/hg38/annotations/hg38_CpG.gz -r chr1:876345-877299
		# equal to:
		mhaptk stat --metrics MM MCR PDR CHALM MHL Entropy MBS --mhapPath HBA_211015_H1930004_CB63_V1C_3C_1_P5-6-K9-G24.mC.mhap.gz --cpgPath ~/Ref/hg38/annotations/hg38_CpG.gz --outputFile 2.txt --region chr1:10813-17492 && cat 2.txt
		# or: java version stat is not working.
		java -jar ~/Software/mHapSuite-2.0-alpha/target/mHapSuite-2.0-jar-with-dependencies.jar stat -metrics MM PDR CHALM MHL MCR MBS Entropy -mhapPath HBA_211015_H1930004_CB63_V1C_3C_1_P5-6-K9-G24.mC.mhap.gz -cpgPath ~/Ref/hg38/annotations/hg38_CpG.gz -outputFile 3.txt -region chr1:876345-877299
	Parameters
	----------
	mhap_path :
	cpg_path :
	region :
	bed :
	output :

	Returns
	-------

	"""
	assert not cpg_path is None
	header=['name','chrom','start','end','nReads','mBase','cBase','tBase',
				'K4plus','nDR','nMR','MM',
				'CHALM','PDR','MHL','MBS','MCR','entropy','nCPG']
	f_cpg = pysam.TabixFile(os.path.expanduser(cpg_path))
	f_mhap = pysam.TabixFile(os.path.expanduser(mhap_path))
	if not region is None:
		assert isinstance(region,str)
		chrom = region.split(':')[0]
		start, end = region.split(':')[1].split('-')
		start,end=int(start)-1,int(end)
		r=stat_single_region(f_mhap,chrom,start,end)
		cpgs = f_cpg.fetch(chrom, start=start, end=end)
		r.append(len(list(cpgs)))
		f_cpg.close()
		f_mhap.close()
		df=pd.DataFrame([[region,chrom,start,end]+r],columns=header)
		# print(df)
		return df
	if not bed is None:
		f_bed = open1(os.path.expanduser(bed))
		results=[]
		for line in f_bed.readlines():
			# print(line)
			line = line.decode('utf-8')
			values=line.split('\t')
			chrom=values[0]
			start=int(values[1])
			end=int(values[2])
			name=values[3]
			r = stat_single_region(f_mhap, chrom, start, end)
			cpgs = f_cpg.fetch(chrom, start=start, end=end)
			r.append(len(list(cpgs)))
			results.append([name,chrom,start,end]+r)
		f_cpg.close()
		f_bed.close()
		f_mhap.close()
		df=pd.DataFrame(results,columns=header)
		# df=df.loc[df.nReads > 0] #remove records with no coverage
		if output is None:
			# print(df)
			return df
		else:
			df.to_csv(os.path.expanduser(output),sep='\t',index=False)

if __name__ == "__main__":
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(serialize=lambda x:print(x) if not x is None else print(""))