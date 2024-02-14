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