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