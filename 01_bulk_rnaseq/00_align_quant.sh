#!/usr/bin/env bash
# =============================================================================
# 00_align_quant.sh -- bulk RNA-seq alignment and gene-level quantification.
#
# Steps: FastQC -> HISAT2 (splice-aware) -> sorted BAM -> featureCounts, producing
# the count matrix consumed by 01_deseq2_deg_enrich.R.
#
# Usage
#   ./00_align_quant.sh <fastq_dir> <reference_dir> <out_dir>
#
#   <fastq_dir>      contains <sample>_R1.fq.gz / <sample>_R2.fq.gz
#   <reference_dir>  contains hisat2 index prefix `genome.*` and `genes.gtf`
#   <out_dir>        receives BAMs, QC and all_samples.gene_counts.txt
#
# The sample list is read from the sample sheet (column 1), so the count matrix
# column order always matches the sheet used by 01_deseq2_deg_enrich.R:
#   SAMPLES=samples/BULK-BRAIN.tsv     ./00_align_quant.sh ...   (4 vs 4, 24 h)
#   SAMPLES=samples/BULK-NEU-basal.tsv ./00_align_quant.sh ...   (neutrophils, basal)
#   SAMPLES=samples/BULK-NEU-3h.tsv    ./00_align_quant.sh ...   (neutrophils, 3 h)
# The two neutrophil sheets are subsets of one counts matrix; run alignment once.
#
# NOTE ON STRANDEDNESS: `-s 2` assumes a reverse-stranded library (SMARTer Stranded
# Total RNA-Seq, used here). Using `-s 0` halves the usable reads and shifts every
# downstream number. Confirm with `infer_experiment.py` if you change kits.
# =============================================================================
set -euo pipefail

FASTQ_DIR=${1:?usage: 00_align_quant.sh <fastq_dir> <reference_dir> <out_dir>}
REF_DIR=${2:?usage: 00_align_quant.sh <fastq_dir> <reference_dir> <out_dir>}
OUT_DIR=${3:?usage: 00_align_quant.sh <fastq_dir> <reference_dir> <out_dir>}
SAMPLES=${SAMPLES:-samples/BULK-BRAIN.tsv}
THREADS=${THREADS:-16}

mkdir -p "${OUT_DIR}/bam" "${OUT_DIR}/qc"

# 1. Read QC -------------------------------------------------------------------
while IFS=$'\t' read -r sample _; do
    [ "${sample}" = "sample_id" ] && continue
    echo "[align] fastqc ${sample}"
    fastqc -t 4 -o "${OUT_DIR}/qc" \
        "${FASTQ_DIR}/${sample}_R1.fq.gz" "${FASTQ_DIR}/${sample}_R2.fq.gz"
done < "${SAMPLES}"

# 2. Alignment -----------------------------------------------------------------
while IFS=$'\t' read -r sample _; do
    [ "${sample}" = "sample_id" ] && continue
    echo "[align] hisat2 ${sample}"
    hisat2 -p "${THREADS}" -x "${REF_DIR}/genome" \
        -1 "${FASTQ_DIR}/${sample}_R1.fq.gz" \
        -2 "${FASTQ_DIR}/${sample}_R2.fq.gz" \
      | samtools sort -@ 8 -o "${OUT_DIR}/bam/${sample}.bam" -
    samtools index "${OUT_DIR}/bam/${sample}.bam"

    samtools flagstat "${OUT_DIR}/bam/${sample}.bam" > "${OUT_DIR}/qc/${sample}.flagstat"
done < "${SAMPLES}"

# 3. Gene-level counts ---------------------------------------------------------
# -p paired-end, --countReadPairs count fragments (not reads), -s 2 reverse-stranded
featureCounts -T "${THREADS}" -p --countReadPairs -s 2 \
    -a "${REF_DIR}/genes.gtf" \
    -o "${OUT_DIR}/all_samples.gene_counts.txt" \
    "${OUT_DIR}"/bam/*.bam

echo "[align] done -> ${OUT_DIR}/all_samples.gene_counts.txt"
echo "[align] verify strandedness: the assigned fraction should be > 0.6 in summary"
echo "        ${OUT_DIR}/all_samples.gene_counts.txt.summary"
