| `generate_mhap` | `False` (**default**) / `True` | also generate `mhap/<cell>.CG.mhap.gz` and `mhap/<cell>.CH.mhap.gz` |
### 4. Also generate mhap files
Set `generate_mhap = True`; CpG / CH sites are taken from the genome fasta (`reference_fasta`)
and the chromosomes from `chrom_size_path`:
```ini
[output]
generate_mhap = True
```

```shell
# allc + allc-CGN + mhap, hisat-3n
yap default-mapping-config --mode m3c --barcode_version V2 \
  --hisat3n_dna_ref "~/Ref/hg38/hg38_ucsc_with_chrL" \
  --genome "~/Ref/hg38/hg38_ucsc_with_chrL.fa" \
  --chrom_size_path "~/Ref/hg38/hg38_ucsc.main.chrom.sizes" \
  --mc_format allc --extract_mcg --generate_mhap > m3c_config.ini
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
         outputs={'CGN': 'cell.CG.mhap.gz', 'CHN': 'cell.CH.mhap.gz'})
```
