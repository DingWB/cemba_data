# Installation
## Create environment and install
```shell
conda install -y -n base -c conda-forge mamba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/env.yaml

# if failed, try:
# mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/env_greedy.yaml
conda activate yap

# conda env export > env_greedy.yaml
```
## To install this latest version:
```shell
pip install git+https://github.com/DingWB/cemba_data

# reinstall
pip uninstall -y cemba_data && pip install git+https://github.com/DingWB/cemba_data
```

# Documentation
## 1. Run on local (hisat-3n)
### (1). Make sure create the right environment

### (2). Generate config.ini
```shell
# m3c
yap default-mapping-config --mode m3c --barcode_version V2 --bismark_ref "~/Ref/mm10/mm10_ucsc_with_chrL.bismark1" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --hisat3n_dna_ref  "~/Ref/mm10/mm10_ucsc_with_chrL" > m3c_config.ini

#mC
yap default-mapping-config --mode mc --barcode_version V2 --bismark_ref "~/Ref/mm10/mm10_ucsc_with_chrL.bismark1" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --hisat3n_dna_ref  "~/Ref/mm10/mm10_ucsc_with_chrL" > mc_config.ini
# pay attention to the path of reference, should be the same as on the GCP if you are going to run the pipeline on GCP.    

# mct
# bismark & STAR for mct (bowtie2)
yap default-mapping-config --mode mct --barcode_version V2 --bismark_ref "~/Ref/mm10/mm10_ucsc_with_chrL.bismark2" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --gtf "~/Ref/mm10/annotations/gencode.vM23.annotation.gtf" --star_ref "~/Ref/mm10/star_ref" > mct_config.ini

# hisat-3n for mct    
yap default-mapping-config --mode mct --barcode_version V2 --hisat3n_dna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --hisat3n_rna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --gtf "~/Ref/mm10/annotations/gencode.vM23.annotation.gtf" > mct_config.ini
```

### (3). Demultiplex
```shell
# m3c
yap demultiplex --fastq_pattern "Pool_Remind1_m3c/*.fastq.gz" -o mapping/Pool_Remind1_m3c -j 16 --aligner hisat3n --config_path m3c_config.ini
# or
 yap-gcp run_demultiplex --fq_dir="Pool_Remind1_m3c" --outdir="mapping/Pool_Remind1_m3c" --gcp=False --n_jobs=16 --print_only=True 
 
# mc
 yap-gcp run_demultiplex --fq_dir="Pool_Remind1_mC" --outdir="mapping/Pool_Remind1_mC" --gcp=False --n_jobs=16 --print_only=True 
 
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

### (4). Run mapping
### Run on local computer or HPC
```shell
sh mapping/snakemake/qsub/snakemake_cmd.txt
# or
yap-gcp run_mapping --fastq_prefix="mapping/Pool_Remind1_m3c" --gcp=False --config_path="m3c_config.ini" --aligner='hisat-3n' --n_jobs=64 --print_only=True
# or bismark
yap-gcp run_mapping --fastq_prefix="mapping " --gcp=False --config_path="m3c_config.ini" --aligner='bismark' --n_jobs=64
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