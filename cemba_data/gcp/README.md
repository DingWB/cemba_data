# Run on GCP
```shell
yap-gcp get_demultiplex_skypilot_yaml > demultiplex.yaml # vim
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


# Run YAP pipeline on test datasets
### Download example fastq (cell level)
```shell
pip install pyfigshare
# setup the token: https://github.com/DingWB/pyfigshare?tab=readme-ov-file#1-setup-token

figshare download 26210798 -o yap_example -c 2 -f fastq
cd yap_example
mkdir -p bismark_mapping/bismark/fastq
mkdir -p hisat3n_mapping/hisat3n/fastq
cwd=$(pwd)
for fq in `ls ${cwd}/fastq/*.fq.gz | grep -v "trimmed"`; do
  file=$(basename ${fq})
  ln -s ${fq} ${cwd}/bismark_mapping/bismark/fastq/${file}
  ln -s ${fq} ${cwd}/hisat3n_mapping/hisat3n/fastq/${file}
done;
```

## Prepare mapping config files
```
yap default-mapping-config --mode m3c-mhap --barcode_version V2 --bismark_ref "~/Ref/hg38/hg38_ucsc_with_chrL.bismark1" --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --annotation_path "~/Ref/hg38/annotations/hg38_allc.gz" > m3c-mhap_config_bismark.ini

yap default-mapping-config --mode m3c --barcode_version V2 --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --hisat3n_dna_ref  "~/Ref/hg38/hg38_ucsc_with_chrL" > m3c_config_hisat3n.ini
# or 
yap default-mapping-config --mode m3c-mhap --barcode_version V2 --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --hisat3n_dna_ref  "~/Ref/hg38/hg38_ucsc_with_chrL" --annotation_path "~/Ref/hg38/annotations/hg38_allc.gz" > m3c-mhap_config_hisat3n.ini
```

## Run mapping
```shell
yap-gcp run_mapping --workd="bismark_mapping" --fastq_server="local" --gcp=False --config_path="m3c-mhap_config_bismark.ini" --aligner='bismark' --n_jobs=4 --print_only=True
cat bismark_mapping/snakemake/qsub/snakemake_cmd.txt #  add --notemp to keep all temporary files

yap-gcp run_mapping --workd="hisat3n_mapping" --fastq_server="local" --gcp=False --config_path="m3c-mhap_config_hisat3n.ini" --aligner='hisat3n' --n_jobs=4 --print_only=True
cat hisat3n_mapping/snakemake/qsub/snakemake_cmd.txt # sh to run
```


# Run yap-gcp on fastq stored on SRA / GEO or other ftp server
```shell
figshare download 26210798 -f HBA_snm3C_phenotype.tsv
```

download only two region (donor1 V1C and BNST)
```python
df=pd.read_csv("HBA_snm3C_phenotype.tsv",sep='\t')
df=df.loc[(df.source_name_ch1=='h1930001') & (df['brain region'].isin(['V1C','BNST']))]
# randomly  select 2 cells from each region:
df=df.groupby('brain region').sample(2)
for donor,df1 in df.groupby('source_name_ch1'):
    outdir=os.path.abspath(donor)
    for region,df2 in df1.groupby('brain region'):
        print(donor,region)
        outdir=os.path.join(donor,region)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        df2=df2.loc[:,['title','R1_ftp','R2_ftp']].set_index('title').stack().reset_index()
        df2.columns=['cell_id','read_type','fastq_path']
        df2.read_type=df2.read_type.apply(lambda x:x.split('_')[0])
        df2.to_csv(os.path.join(outdir,"CELL_IDS"),sep='\t',index=False)
```

## Run mapping (download cell fastq directly from ftp server and delete it after it is no longer needed)
```shell
yap-gcp run_mapping --workd="h1930001" --fastq_server='ftp' --gcp=False --config_path="m3c-mhap_config_hisat3n.ini" --aligner='hisat3n' --n_jobs=4 --total_memory_gb=20 --print_only=True
```


# Generate DAG graph
```shell
mamba install graphviz #required to run dot
snakemake --config fastq_server='ftp' --dag mhap/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.mhap.gz allc/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.allc.tsv.gz hic/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.hisat3n_dna.all_reads.3C.contact.tsv.gz > 1
dot -Tsvg 1 > snm3c_dag.svg
```