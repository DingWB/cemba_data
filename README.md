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
```shell
#yap-gcp get_demultiplex_skypilot_yaml > skypilot.yaml
yap-gcp prepare_demultiplex ----fq_dir gs://mapping_example/fastq/test_fastq \
              --remote_prefix mapping_example --outdir test_gcp \
              --env_name yap --n_jobs 96 --output run_demultiplex.yaml
# vim and change config in run_demultiplex.yaml
sky launch -y -n demultiplex run_demultiplex.yaml # Do Not use spot mode.
```