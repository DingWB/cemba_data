"""
BAM-only hisat-3n template: (re)generate allc / cz / allc-CGN / mhap / 3C contacts from
existing cell BAM files, without demultiplex, trimming or mapping. This is the default
template of `yap-gcp start_from_cell_bam`; it can also be run directly on a directory
that contains bam/{cell_id}...bam:

    snakemake -d <dir> -s from_bam.smk -j 16 \
        --config input_bam="bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam" mc_format=cz generate_mhap=True

What is generated is controlled by the mapping config (mapping_config.ini) or --config
(start_from_cell_bam writes its overrides at the top of the Snakefile):
  input_bam      BAM path pattern with {cell_id}, under bam/. A pattern ending with
                 ".deduped.bam" is used as is; any other BAM (e.g. the name-sorted
                 *.hisat3n_dna.all_reads.name_sort.bam kept by the m3c pipeline) is
                 position-sorted and deduplicated first.
  bam_deduped    override the ".deduped.bam" detection (True / False).
  mc_format      allc / cz / both / none ('none': no methylation output, e.g. mhap only).
  extract_mcg    True / False, allc-CGN merged ALLC (needs allc).
  generate_mhap  True / False.
  generate_hic   True / False, 3C contacts; needs the name-sorted (not deduplicated) BAM.
The post-BAM rules are reused from hisat3n.smk, so the outputs are identical to a full
mapping run; only MappingSummary.csv.gz is not produced (its mapping stats do not exist).
"""
import os, sys
import cemba_data
from snakemake.io import glob_wildcards
from cemba_data.hisat3n import read_mapping_config
PACKAGE_DIR = cemba_data.__path__[0]


def _coerce_bool(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ('true', '1', 'yes', 'y', 't', 'on')


# --config beats mapping_config.ini (same precedence as hisat3n_base.smk)
_ini_config = read_mapping_config()
def _cfg(key, default):
    return config.get(key, _ini_config.get(key, default))


INPUT_BAM = str(_cfg("input_bam", "bam/{cell_id}.hisat3n_dna.all_reads.name_sort.bam"))
if "{cell_id}" not in INPUT_BAM or not INPUT_BAM.startswith("bam/"):
    raise ValueError(f"input_bam must be a bam/ path pattern containing {{cell_id}}, got {INPUT_BAM!r}")
if "CELL_IDS" not in globals():  # written at the top of the Snakefile by start_from_cell_bam
    CELL_IDS = list(glob_wildcards(INPUT_BAM).cell_id)
    if len(CELL_IDS) == 0:
        raise ValueError(f"No BAM file matches {INPUT_BAM}")

BAM_DEDUPED = _coerce_bool(_cfg("bam_deduped", INPUT_BAM.endswith(".deduped.bam")))
NO_MC = str(_cfg("mc_format", "cz")).lower() in ("none", "no", "false", "")
if NO_MC:
    # placeholder accepted by hisat3n_base.smk; no methylation target is requested below
    config["mc_format"] = "allc"
GENERATE_HIC = _coerce_bool(_cfg("generate_hic", False))
if GENERATE_HIC and BAM_DEDUPED:
    raise ValueError("generate_hic needs the name-sorted BAM (*.all_reads.name_sort.bam), "
                     "contacts can not be called from a deduplicated, position-sorted BAM")
# required by hisat3n_base.smk, but no mapping is done here
config.setdefault("hisat3n_dna_reference", "not_used_from_bam")

include:
    os.path.join(PACKAGE_DIR, "files", "smk", "hisat3n_base.smk")


def get_hic_targets(cell_ids):
    if not GENERATE_HIC:
        return []
    return (expand("hic/{cell_id}.hisat3n_dna.all_reads.contact_stats.csv", cell_id=cell_ids)
            + expand("hic/{cell_id}.hisat3n_dna.all_reads.3C.contact.tsv.gz", cell_id=cell_ids)
            + expand("hic/{cell_id}.hisat3n_dna.all_reads.dedup_contacts.tsv.gz", cell_id=cell_ids))


TARGETS = ([] if NO_MC else get_methylation_targets(CELL_IDS)) \
          + get_mhap_targets(CELL_IDS) + get_hic_targets(CELL_IDS)
if len(TARGETS) == 0:
    raise ValueError("Nothing to generate: mc_format is 'none' and generate_mhap / generate_hic are False")
if not BAM_DEDUPED:
    TARGETS += expand("bam/{cell_id}.hisat3n_dna.all_reads.deduped.matrix.txt", cell_id=CELL_IDS)

# the first rule is the default target
rule summary:
    input:
        TARGETS
    run:
        shell(f"rm -rf {bam_dir}/temp")


module hisat3n:
    snakefile:
        os.path.join(PACKAGE_DIR, "files", "smk", "hisat3n.smk")
    config: config

if BAM_DEDUPED:
    DEDUP_BAM = INPUT_BAM
    # index next to the BAM (on GCP bam_dir is the local copy of the remote bam/ folder);
    # not temp, so an index that already existed next to the BAM is never removed
    DEDUP_BAI = bam_dir + "/" + os.path.basename(INPUT_BAM) + ".bai"

    rule index_input_bam:
        input:
            bam=DEDUP_BAM
        output:
            bai=local(DEDUP_BAI)
        shell:
            """
            samtools index {input.bam} {output.bai}
            """
else:
    DEDUP_BAM = local(bam_dir + "/{cell_id}.hisat3n_dna.all_reads.deduped.bam")
    DEDUP_BAI = bam_dir + "/{cell_id}.hisat3n_dna.all_reads.deduped.bam.bai"

    use rule sort_bam_by_pos from hisat3n with:
        input:
            bam=INPUT_BAM

    use rule dedup from hisat3n

    use rule index_bam from hisat3n

if GENERATE_HIC:
    use rule call_chromatin_contacts from hisat3n with:
        input:
            bam=INPUT_BAM

use rule unique_reads_allc from hisat3n with:
    input:
        bam=DEDUP_BAM,
        bai=local(DEDUP_BAI)

use rule unique_reads_cgn_extraction from hisat3n

use rule unique_reads_cz from hisat3n with:
    input:
        bam=DEDUP_BAM,
        bai=local(DEDUP_BAI)

use rule unique_reads_mhap from hisat3n with:
    input:
        bam=DEDUP_BAM,
        bai=local(DEDUP_BAI)
