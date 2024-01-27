[](http://www.network-science.de/ascii/)
<pre>
 **    **     **        *******
//**  **     ****      /**////**
 //****     **//**     /**   /**
  //**     **  //**    /*******
   /**    **********   /**////
   /**   /**//////**   /**
   /**   /**     /**   /**
   //    //      //    //
</pre>
[See Documentation](https://hq-1.gitbook.io/mc/)


# Install
To install this latest version:
```shell
pip install git+https://github.com/DingWB/cemba_data
```

# Documentation
1. Run on local
### (1). Make sure create the right environment
```shell
git clone https://github.com/DingWB/cemba_data.git
mamba env create -f cemba_data/env.yaml
pip install pysam==0.20.0
conda activate yap
```
Or directly read from http:
```shell
mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/env.yaml
pip install pysam==0.20.0
conda activate yap
```

### (2). Generate config.ini
```shell
yap default-mapping-config --mode m3c --barcode_version V2 --bismark_ref "~/Ref/hg38/hg38_ucsc_with_chrL.bismark1" \
      --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes"  \
      > config.ini
# pay attention to the path of reference, should be the same as on the GCP if you are going to run the pipeline on GCP.      
```
### (3). Demultiplex
```shell
yap demultiplex --fastq_pattern "test_fastq/*.gz" -o mapping -j 4 --aligner bismark --config_path config.ini

```
### (4). Run mapping
### Run on local computer or HPC
```shell
sh mapping/snakemake/qsub/snakemake_cmd.txt
```

## Run demultiplex on local and mapping on GCP
### (1). Demultiplex is the same
### (2). Mapping on GCP manually
```shell
scp mapping/AMB_220510_8wk_12D_13B_2_P3-5-A11/Snakefile highmem1:~/sky_workdir
scp -r mapping/AMB_220510_8wk_12D_13B_2_P3-5-A11/fastq highmem1:~/sky_workdir
# GCP
mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/env.yaml
mkdir -p ~/Ref && gsutil -m cp -r -n gs://wubin_ref/hg38 ~/Ref
prefix="mapping_example/mapping/test/AMB_220510_8wk_12D_13B_2_P3-6-A11"
snakemake --snakefile ~/sky_workdir/Snakefile -j 8 --default-resources mem_mb=100 --resources mem_mb=50000 --config gcp=True --default-remote-prefix ${prefix} --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote -np
```

### (3) Run on GCP automatically
```shell
wget https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/files/skypilot_template.yaml
# vim skypilot_template.yaml
# change image_id, --default-remote-prefix (such as mapping_example/mapping/test1/{outdir})
yap update-snakemake -o mapping -t skypilot_template.yaml
# spot
#sky spot launch -y mapping/snakemake/gcp/AMB_220510_8wk_12D_13B_2_P3-3-A11.yaml
cnda activate sky
sh mapping/snakemake/gcp/sky_spot.sh
```

## Run both demultiplex and mapping on GCP
### (1). Demultiplex
```shell
wget https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/files/gcp/demultiplex.Snakefile
# Open an GCP VM machine and run the following code:
snakemake --profile -s demultiplex.Snakefile --use-conda \
                  --config gcp=True fq_dir="gs://mapping_example/fastq/test_fastq" -j 8 \
                  --default-remote-prefix mapping_example \
                  --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote -np
# Or paste it into a yaml file and run using skypilot
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
wget https://raw.githubusercontent.com/DingWB/cemba_data/master/cemba_data/files/gcp/merge_lanes.Snakefile

snakemake -s merge_lanes.Snakefile \
                  --use-conda --config gcp=True outdir="test2" -j 8 \
                  --default-remote-prefix mapping_example \
                  --default-remote-provider GS --google-lifesciences-region us-west1 \
                  --default-remote-provider GS --google-lifesciences-region us-west1 --keep-remote
```