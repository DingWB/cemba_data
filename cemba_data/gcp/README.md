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

figshare download 26210798 -o yap_example --cpu 2 --folder fastq
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
yap default-mapping-config --mode m3c --barcode_version V2 --bismark_ref "~/Ref/hg38/hg38_ucsc_with_chrL.bismark1" --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --annotation_path "~/Ref/hg38/annotations/hg38_allc.gz" > m3c_config_bismark.ini

yap default-mapping-config --mode m3c --barcode_version V2 --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --hisat3n_dna_ref  "~/Ref/hg38/hg38_ucsc_with_chrL" > m3c_config_hisat3n.ini
# or, to also generate mhap files (set generate_mhap True and provide annotation_path)
yap default-mapping-config --mode m3c --barcode_version V2 --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" --hisat3n_dna_ref  "~/Ref/hg38/hg38_ucsc_with_chrL" --annotation_path "~/Ref/hg38/annotations/hg38_allc.gz" --generate_mhap True > m3c_mhap_config_hisat3n.ini
```

## Run mapping
```shell
yap-gcp run_mapping --workd="bismark_mapping" --fastq_server="local" --gcp=False --config_path="m3c_config_bismark.ini" --aligner='bismark' --n_jobs=4 --print_only=True
cat bismark_mapping/snakemake/qsub/snakemake_cmd.txt #  add --notemp to keep all temporary files

yap-gcp run_mapping --workd="hisat3n_mapping" --fastq_server="local" --gcp=False --config_path="m3c_mhap_config_hisat3n.ini" --aligner='hisat3n' --n_jobs=4 --print_only=True
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
yap-gcp run_mapping --workd="h1930001" --fastq_server='ftp' --gcp=False --config_path="m3c_mhap_config_hisat3n.ini" --aligner='hisat3n' --n_jobs=4 --total_memory_gb=20 --print_only=True
```


# Generate DAG graph
```shell
mamba install graphviz #required to run dot
snakemake --config fastq_server='ftp' --dag mhap/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.mhap.gz allc/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.allc.tsv.gz hic/HBA_220218_H1930001_CX46_BNST_3C_1_P3-3-O5-K5.hisat3n_dna.all_reads.3C.contact.tsv.gz > 1
dot -Tsvg 1 > snm3c_dag.svg
```

# Predict Neuron Percentage
pip install scikit-learn==1.2.2
```python
def predict_neuron_percentage(
    model_path = '/gale/ddn/bican/snm3C/miseq/dev/rfc_model.hba-miseq.joblib',
    infile="results/MappingSummary-SALK090.csv.gz",
):
    import joblib
    import numpy as np
    import pandas as pd
    import sklearn # 1.2.2

    predictor = joblib.load(model_path)
    # debug
    df=pd.read_csv(infile,index_col=0)
    tmp = df[['mCGFrac','mCHFrac','mCCCFrac']].dropna()
    pred = pd.Series(predictor.predict(tmp), index = tmp.index)
    df['NeuN-alike'] = ~pred
    df.loc[df['FinalmCReads']<20, 'NeuN-alike'] = np.nan
    neun_ratios = df[['Plate','NeuN-alike']].pivot_table(index='Plate', columns='NeuN-alike', aggfunc=len, fill_value=0)
    tot_neun_ratio = neun_ratios[True].sum()/neun_ratios.sum().sum() #neun_ratios is a df,
    return tot_neun_ratio
```

# Changelog

## Fix: hisat-3n mapping-rate denominator (`cell_parser_hisat_summary`)

**File:** `cemba_data/hisat3n/stats_parser.py`

**What changed**

```python
# before
total_reads = report_dict['ReadPairsMappedInPE'] * 2 + report_dict['ReadsMappedInSE']
# after
total_reads = report_dict['ReadPairsMappedInPE'] * 2
```

`total_reads` is the denominator for `UniqueMappingRate`, `MultiMappingRate`,
and `OverallMappingRate`.

**Why**

The paired-end hisat-3n mapping runs in mixed mode. Mixed mode is the default
in hisat-3n/bowtie2 and is only disabled by passing `--no-mixed`; the pipeline
does not pass that flag, so mixed mode stays on. As a result the `--new-summary`
output has two blocks:

- `Total pairs: P`   -> `ReadPairsMappedInPE`
- `Total unpaired reads: U`  -> `ReadsMappedInSE`

In mixed mode, when a pair fails to align concordantly/discordantly, hisat-3n
splits it and re-tries each mate as an unpaired read. Therefore `U = 2 * Z`,
where `Z` is the number of pairs that failed to align as a pair. These unpaired
reads are **not new reads**; they are a subset of the original `P * 2` mates,
re-attempted individually.

The total number of distinct input mates satisfies the invariant:

```
P * 2 = (P - Z) * 2   [mates aligned as pairs]
      + U             [mates re-tried as unpaired], with U = 2 * Z
```

So the correct denominator is `P * 2`. The old formula `P * 2 + U` double-counts
the failed-pair mates (`2 * Z`), inflating the denominator and systematically
**deflating every mapping rate**.

**Impact by pipeline**

The fix lives in a single shared parser (`cell_parser_hisat_summary`), so it
applies everywhere that function is used:

- `mc`  -> DNA summary (`*.hisat3n_dna_summary.txt`)
- `mct` -> DNA summary + RNA summary (`*.hisat3n_rna_summary.txt`)
- `m3c` -> DNA summary (`*.hisat3n_dna_summary.txt`)

For `mc` / `mct`, most pairs map concordantly, so `Z` (and the inflation) is
small and the previous under-reporting was mild. For `m3c`, many read pairs span
a chromatin ligation junction and fail PE mapping, so `Z` is large and the
under-reporting was severe (e.g. `UniqueMappingRate` around 35%). Note that for
`m3c` the PE `UniqueMappingRate` is still expected to be low by design: the
reads that fail PE mapping are split at enzyme cut sites and recovered later by
the single-end split-read re-alignment (`*.hisat3n_dna_split_reads_summary.R1/R2.txt`),
which is independent of this denominator fix.

The single-end split-read parser (`cell_parser_hisat_se_summary`) uses
`Total reads` as its denominator and was never affected.
