### (1). Demultiplex
```shell
wget https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/gcp/smk/demultiplex.Snakefile
# Open an GCP VM machine and run the following code:
snakemake -s demultiplex.smk --use-conda \
                  --config gcp=True fq_dir="gs://mapping_example/fastq/test_fastq" outdir="test2" -j 8 \
                  --default-remote-prefix mapping_example \
                  --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote
# Or paste it into a yaml file and run using skypilot

# if Run GCP pipeline on local:
snakemake -s demultiplex.smk --use-conda \
                  --config fq_dir="/anvil/projects/x-mcb130189/Wubin/BICAN/test_pipeline/test_fastq" outdir="test2" -j 8
```
```text
name: demultiplex
workdir: .
num_nodes: 1
resources:
    cloud: gcp
    region: us-west1
    instance_type: n1-standard-8
    use_spot: True
    disk_size: 250
    disk_tier: 'medium'
    image_id: projects/ecker-wding/global/images/myimage

# file_mounts:
#   ~/Ref/hg38: ~/Ref/hg38

setup: |
  # pip install --upgrade pip
  # conda install -y -n base -c conda-forge -c bioconda mamba
  # mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/env.yaml
  mkdir -p ~/Ref && gsutil -m cp -r -n gs://wubin_ref/hg38 ~/Ref

run: |
  conda activate yap
  pip install git+https://github.com/DingWB/cemba_data
  snakemake -s demultiplex.Snakefile --use-conda \
                  --config gcp=True fq_dir="gs://mapping_example/fastq/test_fastq" outdir="test2" -j 8 \
                  --default-remote-prefix mapping_example \
                  --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote
```

```shell
sky spot launch -n demultiplex -y run_demultiplex.yaml
```

### (2). Merge lanes
```shell
wget https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/gcp/smk/merge_lanes.Snakefile

snakemake -s merge_lanes.smk \
                  --use-conda --config gcp=True outdir="test2" -j 8 \
                  --default-remote-prefix mapping_example \
                  --default-remote-provider GS --google-lifesciences-region us-west1 \
                  --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote
                  
# if Run GCP pipeline on local:
snakemake -s merge_lanes.smk --use-conda \
                  --config outdir="test2" -j 8
```


## 2 Run on GCP
```shell
yap-gcp yap_pipeline --fq_dir="gs://mapping_example/fastq/novaseq_fastq" \
--remote_prefix='mapping_example' --outdir='novaseq_mapping' --env_name='yap' \
--n_jobs1=16 --n_jobs2=60 \
--image="bican" --n_node 1 --disk_size1 300 --disk_size2 300 \
--demultiplex_template="~/Projects/BICAN/yaml/demultiplex.yaml" \
--mapping_template="~/Projects/BICAN/yaml/mapping.yaml" \
--genome="~/Ref/hg38/hg38_ucsc_with_chrL.fa" \
--hisat3n_dna_ref="~/Ref/hg38/hg38_ucsc_with_chrL" \
--mode='m3c' --bismark_ref='~/Ref/hg38/hg38_ucsc_with_chrL.bismark1' \
--chrom_size_path='~/Ref/hg38/hg38_ucsc.main.chrom.sizes' \
--aligner='hisat-3n' > run.sh
source run.sh

# demultiplex 
yap-gcp run_mapping --fastq_prefix="mapping/mCT" --gcp=False --config_path="mct_config.ini" --aligner='hisat-3n' --n_jobs=64 --print_only True | grep "^snakemake" > run_mct_mapping.sh
yap sbatch --project_name mapping --command_file_path run_mct_mapping.sh --queue shared --max_jobs 4 --dry_run --working_dir ./ --time_str 1
# or
split -n l/8 -d run_mct_mapping.sh run_mct_mapping.
for i in {0..7}; do 
  echo "run_mct_mapping.0${i}"
  Pbsgen -name mapping_${i} -c 64 -d 2 -m 90 -p shared -submit sh run_mct_mapping.0${i}
#  Pbsgen -name mapping_${i} -c 64 -d 2 -m 90 -p shared sh run_mct_mapping.0${i}
done;
```

# Testing pipeline
## 1.1. Make example fastq files
Randomly sampling 1000000 reads from 4 big fastq files

```shell
seqtk sample -s100 download/UWA7648_CX182024_Idg_1_P1-1-K15_22HC72LT3_S1_L001_R1_001.fastq.gz 1000000 | gzip > novaseq_fastq/UWA7648_CX182024_Idg_1_P1-1-K15_22HC72LT3_S1_L001_R1_001.fastq.gz
# | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" |
```


## 2.1 Run pipeline on GCP
```shell
yap-gcp get_demultiplex_skypilot_yaml > demultiplex.yaml # vim
# demultiplex: n1-highcpu-16
yap-gcp yap_pipeline --fq_dir="gs://mapping_example/fastq/novaseq_fastq" \
--remote_prefix='mapping_example' --outdir='novaseq_mapping' --env_name='yap' \
--n_jobs1=16 --n_jobs2=16 \
--image="bican" --n_node 1 --disk_size1 300 --disk_size2 300 \
--demultiplex_template="~/Projects/BICAN/yaml/demultiplex.yaml" \
--mapping_template="~/Projects/BICAN/yaml/mapping.yaml" \
--genome="~/Ref/hg38_Broad/hg38.fa" \
--hisat3n_dna_ref="~/Ref/hg38_Broad/hg38" \
--mode='m3c' --bismark_ref='~/Ref/hg38/hg38_ucsc_with_chrL.bismark1' \
--chrom_size_path='~/Ref/hg38_Broad/hg38.chrom.sizes' \
--aligner='hisat-3n' > run.sh
source run.sh
```


# Run Salk010 for test (comparing cost with Broad)
```shell
# salk10_test
## 1.1 Run demultiplex on GCP
yap-gcp get_demultiplex_skypilot_yaml > demultiplex.yaml # vim
# demultiplex: n1-highcpu-64
yap-gcp yap_pipeline --fq_dir="gs://mapping_example/fastq/salk10_test" \
--remote_prefix='bican' --outdir='salk010_test' --env_name='yap' \
--n_jobs1=16 --n_jobs2=64 \
--image="bican" --disk_size1 300 --disk_size2 500 \
--demultiplex_template="demultiplex.yaml" \
--mapping_template="mapping.yaml" \
--genome="~/Ref/hg38_Broad/hg38.fa" \
--hisat3n_dna_ref="~/Ref/hg38_Broad/hg38" \
--mode='m3c' --bismark_ref='~/Ref/hg38/hg38_ucsc_with_chrL.bismark1' \
--chrom_size_path='~/Ref/hg38_Broad/hg38.chrom.sizes' \
--aligner='hisat-3n' --n_node=2 > run.sh
	
source run.sh
  
# salk10
# if n2-highcpu-64, use 60 jobs
yap-gcp yap_pipeline --fq_dir="gs://nemo-tmp-4mxgixf-salk010/raw" \
--remote_prefix='bican' --outdir='salk010' --env_name='yap' \
--image="bican" --disk_size1 4096 --disk_size2 260 \
--n_jobs1 16 --n_jobs2 60 \
--demultiplex_template="~/Projects/BICAN/yaml/demultiplex.yaml" \
--mapping_template="~/Projects/BICAN/yaml/mapping.yaml" \
--genome="~/Ref/hg38_Broad/hg38.fa" \
--hisat3n_dna_ref="~/Ref/hg38_Broad/hg38" \
--mode='m3c' --bismark_ref='~/Ref/hg38/hg38_ucsc_with_chrL.bismark1' \
--chrom_size_path='~/Ref/hg38_Broad/hg38.chrom.sizes' \
--aligner='hisat-3n' --n_node 16 > run.sh
# source run.sh
```


## Add mhap
```shell
ava -jar ~/Software/mHapSuite-2.0-alpha/target/mHapSuite-2.0-jar-with-dependencies.jar convert -cpgPath ~/Ref/mm10/annotations/mm10_CpG.gz -inputFile 1.bam -outPutFile 1.mhap.gz
```

# run mapping from fastq
```shell
yap default-mapping-config --mode m3c --barcode_version V2 --genome "/anvil/projects/x-mcb130189/wenliang/hg38/hisat-3N/GRCh38.primary_assembly.genome_withchrL.fa" --chrom_size_path "/anvil/projects/x-mcb130189/wenliang/hg38/GRCh38.primary_withchrL_chrom_sizes" --hisat3n_dna_ref  "/anvil/projects/x-mcb130189/wenliang/hg38/hisat-3N/hg38" > m3c_config.ini

yap-gcp run_mapping --fastq_prefix="mapping" --gcp=False --config_path="m3c_config.ini" --aligner='hisat-3n' --n_jobs=4 --print_only=True

```

### mhap format
```text
chrom, start, end, vector, count, strand
chr1    17452   17562   11111   1       -
chr1    17478   17492   111     1       +
chr1    63627   63683   1111    1       -
chr1    91189   91267   110     1       +
```
```shell
java -jar ${mHapSuite} stat --metrics MM PDR CHALM MHL MCR MBS --cpgPath ${cpgPath} --mhapPath HBA_211129_H1930001_CX45_A46_3C_1_P1-1-K15-A2.mhap.gz -outputFile stat.out.tsv -region chr1:17452-17562

# chr1:17452-17562
chr     start   end     nReads  mBase   cBase   tBase   K4plus  nDR     nMR     nCPG    nPairs  MM      PDR     CHALM   MHL     MCR     MBS
chr1    17452   17562   2       8       0       8       1       0       1       5       0       NaN     NaN     NaN     NaN     NaN     NaN

# tabix HBA_211129_H1930001_CX45_A46_3C_1_P1-1-K15-A2.mhap.gz chr1:959919-959986
chr1    959919  959940  100100  1       +
chr1    959919  960019  100100000000000 1       +
chr1    959948  959986  00000   1       +
java -jar ${mHapSuite} stat --metrics MM PDR CHALM MHL MCR MBS Entropy R2 --cpgPath ${cpgPath} --mhapPath HBA_211129_H1930001_CX45_A46_3C_1_P1-1-K15-A2.mhap.gz -outputFile stat.out.tsv -region chr1:959919-959986
chr     start   end     nReads  mBase   cBase   tBase   K4plus  nDR     nMR     nCPG    nPairs  MM      PDR     CHALM   MHL     MCR     MBS
chr1    959919  959986  3       4       17      26      3       2       2       11      0       NaN     NaN     NaN     NaN     NaN     NaN

java -jar ${mHapSuite} stat --metrics MM PDR CHALM MHL MCR MBS Entropy R2 --cpgPath ${cpgPath} --mhapPath HBA_211129_H1930001_CX45_A46_3C_1_P1-1-K15-A2.mhap.gz -outputFile stat.out.tsv -region chr11:15964448-16741591
chr     start   end     nReads  mBase   cBase   tBase   K4plus  nDR     nMR     nCPG    nPairs  MM      PDR     CHALM   MHL     MCR     MBS     Entropy R2
chr11   15964448        16741591        204     262     23      341     12      1       6       4532    0       0.76832845      0.08333333      0.50000000      0.10091051      0.06744868      0.41692387      0.35904006      NaN
```

column	Description
- nReads	the number of mapped reads: sum of column `count`
- mBase	the methylated CpGs within mapped reads, sum of '1' for all records, including duplicatation.
- cBase	the number of unmethylated CpGs in reads with discordant methylation: No. of 0 in records with at least one 1
- tBase total number of CpGs within mapped reads, inlcuded duplicatation: sum([length(record.vector) for record in records]) 
- K4plus	the number of mapped reads with at least 4 CpGs
- nDR	the number of reads with discordant methylation status among all reads defined by K4plus: len([record.vector if len(set(record.vector))==2 and len(record.vector)>=4])
- nMR	the number of methylated reads among all reads defined by K4plus: len([record.vector if '1' in record.vector and len(record.vector)>=4])
- nCPG Total No. of unique CpG in hg38_CpG.gz
- MM: Mean Methylation: mBase / tBase
- PDR (Proportion of Discordant Reads): nDR / K4plus
- CHALM (Cellular Heterogeneity-Adjusted cLonal Methylation): ratio of methylated reads: nMR / K4plus
- MHL (methylated haplotype load): partitioned each segment into methylation haplotype blocks (MHBs). MHBs were defined as the genomic region in which the r2 value of two adjacent CpG sites is no less than 0.5. MHL=normalized fraction of methylated haplotypes at different lengths

```python
def mhl(vectors):
	mhl_dict={}
	for i in range(1,max([len(vector) for vector in vectors])+1):
		unreachable_n=i-1
		s=sum([len(vector)-unreachable_n for vector in vectors]) #sum
		count=sum([sum([1 if vector[j:j+i] == '1'*i else 0 for j in range(0,len(vector)-unreachable_n)]) for vector in vectors])
		mhl_dict[i]=count / s
	return sum([i*mhl_dict[i] for i in mhl_dict]) / sum([i for i in mhl_dict])
vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
vectors2=['1100']*8+['0011']*8
print(mhl(vectors1))
print(mhl(vectors2))
print(mhl(['0000']*16))
print(mhl(['1111']*16))

def methylation_entropy(vectors):
	df=pd.Series(vectors)
	return -1*df.value_counts(normalize=True).apply(lambda x:x*np.log2(x)).sum() / df.apply(len).max()
vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
vectors2=['1100']*8+['0011']*8
print(entropy(vectors1))
print(entropy(vectors2))
print(entropy(['0000']*16))
print(entropy(['1111']*16))
```

- MCR (methylation concurrence ratio): cBase / tBase; N=sum([vector.count('0') if '1' in vector else 0 for vector in vectors]); N is the sum of No of unmethylated CpGs in partially methylated reads
- MBS (methylation block score): N / nReads
```python
def mbs(vectors,counts):
	results=[]
	for vector,count in zip(vectors,counts): #for each read
	# find the longest successive methylated CpG
		for length in range(len(vector),-1,-1):
			if '1'*length in vector:
				break
		n=len(vector)
		results.append(length*length*count / (n*n))
		print(length,n,count)
	return sum(results) / len(results)
vectors1=['0000','0001','0010','0100','1000','0011','0101','1001','0110','1010','1100','0111','1011','1101','1110','1111']
vectors2=['1100']*8+['0011']*8
print(mbs(vectors1,[1]*len(vectors1)))
print(mbs(vectors2,[1]*len(vectors2)))
```
- R2: 
