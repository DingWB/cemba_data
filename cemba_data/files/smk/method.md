# Methods

## 1. Bisulfite non-conversion rate: lambda DNA spike-in and mCCC

### 1.1 Background

Bisulfite (or enzymatic) conversion deaminates unmethylated cytosine (C) to uracil (U), which is read as thymine (T) after PCR and sequencing, whereas 5-methylcytosine (5mC) is protected and still read as C. Incomplete conversion leaves unmethylated C as C, which looks like methylation and inflates the measured methylation level. This is especially important for non-CpG methylation (mCH), whose genomic level is low.

Two complementary estimates of the non-conversion rate are used:

1. **Lambda DNA spike-in (external control).** Unmethylated lambda phage DNA (Promega D1521) was treated in vitro with the CpG methyltransferase M.SssI (NEB), which methylates (ideally) every CpG. It was then spiked into each snm3C-seq reaction. Because the lambda sequence and its methylation state are known, any deviation from the expected state comes from the conversion chemistry:

   | Lambda site | Expected state | Meaning of the observed C fraction |
   |---|---|---|
   | CH (CA / CC / CT) | Unmethylated, should be 100% converted to T | **Non-conversion rate** |
   | CG | Methylated by M.SssI, should be 100% retained as C | Positive control; $1 - \text{mCG}$ = over-conversion (5mC wrongly deaminated) + incomplete M.SssI methylation |

2. **Genomic mCCC (internal proxy).** In mammalian genomes, cytosines in the CCC trinucleotide context are almost never methylated (DNMT3A mainly methylates CA). The methylation level at CCC therefore approximates the non-conversion rate. It is available in every cell at high coverage and needs no spike-in. It can, however, slightly overestimate the true non-conversion rate in cells with high mCH (e.g. neurons).

The lambda sequence must be included in the reference genome as a contig named **`chrL`** (e.g. `hg38_ucsc_with_chrL.fa`).

### 1.2 Lambda non-conversion rate

For each cell, ALLC records on `chrL` are retrieved (tabix), and each cytosine's context is reduced to a dinucleotide (CA / CC / CT / CG). Methylated base calls ($mC$) and total coverage ($Cov$) are summed per context. The non-conversion rate is

$$
r_{\lambda} = \frac{\sum_{c \in CH} mC_c}{\sum_{c \in CH} Cov_c}, \qquad CH \in \{CA, CC, CT\}
$$

The Bismark pipeline (`get_allc_lambda_frac` in `cemba_data/mapping/stats/utilities.py`) uses only the CY contexts (CC + CT):

$$
\text{LambdaCYFrac} = \frac{mC_{CC} + mC_{CT}}{Cov_{CC} + Cov_{CT}}, \qquad
\text{LambdaCYCov} = Cov_{CC} + Cov_{CT}
$$

The M.SssI-methylated CpGs on lambda give a conversion-specificity (positive) control:

$$
\text{mCG}_{\lambda} = \frac{mC_{CG}}{Cov_{CG}}, \qquad
\text{over-conversion} + \text{incomplete M.SssI methylation} \approx 1 - \text{mCG}_{\lambda}
$$

The hisat-3n pipeline (`cell_parser_allc_lambda` in `cemba_data/hisat3n/stats_parser.py`) reports `Lambda{X}mC`, `Lambda{X}Cov` and `Lambda{X}Frac` for $X \in \{CA, CC, CT, CH, CY, CG\}$, with $\text{Lambda}X\text{Frac} = mC_X / (Cov_X + 10^{-5})$. In NOMe mode, GpC sites (G immediately upstream) are excluded. The per-context CA / CC / CT fractions should agree; a higher CA fraction suggests that host-genome reads (carrying genuine mCA) are mis-mapped to `chrL`, and in that case LambdaCYFrac is the safer estimate. These columns need tabix-indexed ALLC and are not produced when `mc_format = cz`.

### 1.3 mCCC

Per-cell context counts (`allc/*.allc.tsv.gz.count.csv` or `cz/*.cz.count.csv`) are aggregated over all cytosines in the CCC context (HCCC in NOMe mode, to exclude GpC sites methylated by M.CviPI):

$$
\text{mCCCFrac} = \frac{\text{mCCCmC}}{\text{mCCCCov} + 10^{-5}}
$$

This is implemented in `cell_parser_allc_count` (`cemba_data/hisat3n/stats_parser.py`) for the hisat-3n pipeline, and in `generate_allc_stats` (`cemba_data/mapping/stats/utilities.py`, `mc_stat_feature = CHN CGN CCC`) for the Bismark/ALLCools pipeline.

### 1.4 Correcting methylation levels

Given a non-conversion rate $r$ (either $r_{\lambda}$ or mCCCFrac), the observed methylation fraction can be corrected as

$$
\text{mC}_{\text{corrected}} = \frac{\text{mC}_{\text{raw}} - r}{1 - r}
$$
