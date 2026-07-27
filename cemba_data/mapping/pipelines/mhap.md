| `generate_mhap` | `False` (**default**) / `True` | also generate `mhap/<cell>.CG.mhap.gz` and `mhap/<cell>.CH.mhap.gz` |
| `annotation_path` | path to `*_allc.gz` | **required** when `generate_mhap = True` |
### 4. Also generate mhap files
Set `generate_mhap = True` and provide the `*_allc.gz` annotation:
```ini
[output]
generate_mhap = True
annotation_path = ~/Ref/hg38/annotations/hg38_allc.gz
```





```shell
# allc + allc-CGN + mhap, hisat-3n
yap default-mapping-config --mode m3c --barcode_version V2 \
  --hisat3n_dna_ref "~/Ref/hg38/hg38_ucsc_with_chrL" \
  --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" \
  --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" \
  --mc_format allc --extract_mcg True \
  --generate_mhap True --annotation_path "~/Ref/hg38/annotations/hg38_allc.gz" > m3c_config.ini
```