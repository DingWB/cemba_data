# Installation
## Create environment and install
```shell
conda install -y -n base -c conda-forge mamba
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

mamba env create -f https://raw.githubusercontent.com/DingWB/cemba_data/master/yap.yaml

conda activate yap
```
## To install this latest version:
```shell
pip install git+https://github.com/DingWB/cemba_data

# reinstall
pip uninstall -y cemba_data && pip install git+https://github.com/DingWB/cemba_data
```

# Documentation

The `yap-gcp` workflow has three steps: **(1) create a mapping config, (2) demultiplex, and (3) map**. The examples below use `m3c` and `hisat-3n`. Replace `m3c` with another supported mode, such as `mc` or `mct`, when needed. `hisat-3n` is the default aligner; to use Bismark, provide `--bismark_ref` when creating the config and set `--aligner=bismark` when running mapping.

## 1. Create a mapping config

Create the config once. Replace the example reference paths with your own:
```shell
yap default-mapping-config --mode m3c --barcode_version V2 \
  --hisat3n_dna_ref ~/Ref/hg38/hg38_ucsc_with_chrL \
  --genome_fasta ~/Ref/hg38/hg38_ucsc_with_chrL.fa \
  --chrom_size_path ~/Ref/hg38/hg38_ucsc.main.chrom.sizes \
  > m3c_config.ini
# for bismark, add: --bismark_ref "~/Ref/mm10/mm10_ucsc_with_chrL.bismark1"

# mct
# bismark & STAR for mct (bowtie2)
yap default-mapping-config --mode mct --barcode_version V2 --bismark_ref "~/Ref/mm10/mm10_ucsc_with_chrL.bismark2" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --gtf "~/Ref/mm10/annotations/gencode.vM23.annotation.gtf" --star_ref "~/Ref/mm10/star_ref" > mct_config.ini

# hisat-3n for mct    
yap default-mapping-config --mode mct --barcode_version V2 --hisat3n_dna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --hisat3n_rna_ref "~/Ref/mm10/mm10_ucsc_with_chrL" --genome "~/Ref/mm10/mm10_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/mm10/mm10_ucsc.nochrM.sizes" --gtf "~/Ref/mm10/annotations/gencode.vM23.annotation.gtf" > mct_config.ini

```

### Methylation output options

Set these options in `[output]` in the `.ini`, or pass them to `yap default-mapping-config`:

| Parameter | Options | Meaning |
| --- | --- | --- |
| `mc_format` | `allc`, `cz`, `both` | `allc` writes ALLC files; `cz` writes cytozip `.cz` files; `both` writes both. The config-generation command defaults to `allc`. |
| `reference_cz` | path to `.cz` | Required when `mc_format` is `cz` or `both`. Build it once per reference genome (command below), then set its path in the config. |
| `extract_mcg` | `True`, `False` | Generate the CGN-merged ALLC (`allc-CGN/`). Off by default; only applies when ALLC is generated. |

Build a reference `.cz` if using `mc_format=cz` or `both`:
```shell
czip build_ref -g ~/Ref/hg38/hg38_ucsc_with_chrL.fa \
  -O ~/Ref/hg38/hg38_ucsc_with_chrL.allc.cz \
  -s ~/Ref/hg38/hg38_ucsc.main.chrom.sizes -j 20
```

Then set options in `[output]` in the config:
```ini
[output]
mc_format = both # or cz if you don't want to generate allc format
reference_cz = ~/Ref/hg38/hg38_ucsc_with_chrL.allc.cz
extract_mcg = False
```
Alternatively, add `--mc_format both --reference_cz ~/Ref/hg38/hg38_ucsc_with_chrL.allc.cz` to the config-generation command in step 1. Use `--mc_format allc` for ALLC only. Add `--extract_mcg` to enable CGN extraction; it is not run by default.

## Local workflow
`yap-gcp` support both local and GCP

### 2. Demultiplex
```shell
yap-gcp run_demultiplex --fq_dir=/path/to/raw_fastq/Pool1 \
  --outdir=mapping/Pool1 --n_jobs=16 --print_only=True 
# --print_only only print the CMD, please run the printed CMD yourself
# suggestion: 16 CPUs with 120G memory
```

### 3. Mapping
Generate the local Snakemake commands, then run them:
```shell
yap-gcp run_mapping --workd mapping/Pool1 \
  --config_path=m3c_config.ini --aligner="hisat-3n" \
  --n_jobs=64 --print_only=True
sh mapping/Pool1/snakemake/qsub/snakemake_cmd.txt
```
For Bismark, use `--aligner=bismark`.

## GCP workflow

Upload raw FASTQ files to a GCS bucket. Configure SkyPilot and GCP credentials first; the reference paths in the config must be accessible to the GCP workers.

### 2. Demultiplex on GCP
Prepare and launch the demultiplex job. Its output will be `gs://my-bucket/m3c_run`:
```shell
yap-gcp prepare_demultiplex --fq_dir=gs://my-bucket/raw_fastq \
  --remote_prefix=my-bucket --outdir=m3c_run \
  --fastq_server=gcp --gcp=True --region=us-west1 \
  --env_name=yap --n_jobs=16
sky launch -y -n demultiplex run_demultiplex.yaml
```

### 3. Mapping on GCP
Prepare the mapping jobs from the demultiplexed GCS data, then launch them:
```shell
yap-gcp prepare_mapping --workd=gs://my-bucket/m3c_run \
  --config_path=m3c_config.ini --aligner=hisat-3n \
  --fastq_server=gcp --gcp=True --region=us-west1 \
  --env_name=yap --n_node=1 --n_jobs=64
sh mapping_gcp_tmp/run_mapping.sh
```
For Bismark, use `--aligner=bismark`.

## Start from BAM

Use `start_from_cell_bam` to skip demultiplexing and alignment and run post-BAM outputs on existing cell-level BAMs. The config references must match the genome used to align the BAMs.

Local example:
```shell
yap-gcp start_from_cell_bam --indir=/path/to/bams \
  --bam_pattern="*.hisat3n_dna.all_reads.name_sort.bam" \
  --output_dir mapping --config_path=m3c_config.ini \
  --n_jobs=16 --print_only=True
# --outputs will overwrite some config in .ini
sh mapping/snakemake/qsub/snakemake_cmd.txt
```

`--outputs` accepts `allc`, `cz`, `mcg`, and `hic`; for example, `allc,cz` produces both methylation formats. `mcg` enables CGN extraction and implies `allc`; `hic` requires the name-sorted BAM. Omit `--outputs` to follow the output settings in the config. For GCP, set `--gcp=True` and use an `--output_dir` GCS path that already contains `bam/<bam_pattern>`.
