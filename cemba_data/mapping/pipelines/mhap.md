| `generate_mhap` | `False` (**default**) / `True` | also generate `mhap/<cell>.CG.mhap.gz` and `mhap/<cell>.CH.mhap.gz` |
| `mhap_cache_dir` | `None` (**default**) / a directory | cache the per-chromosome site flags of `reference_fasta` (1 byte / bp, ~3 GB for hg38). Without it every cell recomputes them (~20-30 s for hg38); the cache is shared by all cells of the same genome |
### 4. Also generate mhap files
Set `generate_mhap = True`; CpG / CH sites are taken from the genome fasta (`reference_fasta`)
and the chromosomes from `chrom_size_path`:
```ini
[output]
generate_mhap = True
; optional, recommended for many cells: cache the genome site flags once (~3 GB for hg38)
mhap_cache_dir = ~/Ref/hg38/mhap_cache
```

```shell
# allc + allc-CGN + mhap, hisat-3n
yap default-mapping-config --mode m3c --barcode_version V2 \
  --hisat3n_dna_ref "~/Ref/hg38/hg38_ucsc_with_chrL" \
  --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" \
  --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" \
  --mc_format allc --extract_mcg --generate_mhap --mhap_cache_dir ~/Ref/hg38/mhap_cache > m3c_config.ini
```

### mhap conventions (`cemba_data.mapping.pipelines.bam2mhap`)
Same format as mHapTools / mHapSuite: `chrom, start, end, haplotype, count, strand`, where
start / end are the 1-based positions of the first / last site and the haplotype has one
0/1 per consecutive site.

- CG (`*.CG.mhap.gz`): sites are the top-strand C of each CpG on both strands
  (bottom-strand reads are called on the G), so the files work with mHapTk / mHapSuite and a CpG list.
- CH (`*.CH.mhap.gz`): strand-specific sites; a bottom-strand site is at the bottom-strand C (top-strand G).
- Reads: unmapped, secondary, supplementary, QC-fail, duplicate and MAPQ < 10 reads are skipped;
  calls need base quality >= 20 and are placed with the CIGAR (soft clips, insertions, deletions).
- Calls are skipped when the read breaks the context (e.g. the G of a CpG read as A: SNP).
- Mates and m3C split parts of a fragment are merged, so overlapping bases are counted once;
  conflicting calls between them become no-calls.
- A record is cut where sites were not sequenced (between mates). A sequenced site without a call
  between two calls drops the fragment for CG (as mHapTools) and cuts it for CH.

Standalone use (one pass for CG and CH):
```python
from cemba_data.mapping.pipelines import bam2mhap
bam2mhap(bam_path='cell.deduped.bam', reference_fasta='hg38.fa', chroms='hg38.main.chrom.sizes',
         outputs={'CGN': 'cell.CG.mhap.gz', 'CHN': 'cell.CH.mhap.gz'},
         cache_dir='~/Ref/hg38/mhap_cache')  # optional
```
Or from BAMs only (no fastq / allc): `yap-gcp bam2mhap --bam_path ... --reference_fasta ...`, or for a
directory of cells `yap-gcp start_from_cell_bam --snakemake_template files/smk/hisat3n/bam2mhap.smk`
(reads `reference_fasta`, `chrom_size_path` and `mhap_cache_dir` from `mapping_config.ini`).
