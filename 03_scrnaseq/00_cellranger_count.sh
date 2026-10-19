#!/usr/bin/env bash
# =============================================================================
# 00_cellranger_count.sh -- 10x demultiplexing and per-library count matrices.
#
# Produces one filtered feature-barcode matrix and one raw-matrix loom per
# library. 01_qc_integrate_annotate.py reads the loom; the annotated object it
# writes is what CellChat consumes in module 4.
#
# Usage
#   LIBS="NS1 YY1" ./00_cellranger_count.sh <transcriptome_ref> <out_dir>
#
#   <transcriptome_ref>  Cell Ranger reference, e.g. refdata-gex-mm10-2020-A
#   <out_dir>            receives cellranger/<LIB>/outs
#
# Expected input layout:
#   fastq/<LIB>/<LIB>_S*_L00*_R{1,2}_001.fastq.gz
#
# LIBRARY NAMES vs GENOTYPES. The two 10x libraries carry the legacy submission
# codes NS1 (WT) and YY1 (Gp6-KO). Cell Ranger ids stay as sequenced -- renaming a
# library after `cellranger count` invalidates every downstream path -- and the
# analysis labels are declared when the Python script is called:
#
#   --ctrl-loom cellranger/NS1/outs/*_raw.loom  --ctrl-label WT --ctrl-name WT
#   --case-loom cellranger/YY1/outs/*_raw.loom  --case-label KO --case-name Gp6-KO
#
# 01_qc_integrate_annotate.py builds cell barcodes as <barcode>_<label>; the label
# must match the one recorded with the run. Changing it later breaks the joins.
#
# LOOM OUTPUT. Module 3 reads `outs/raw_feature_bc_matrix.loom`. The Cell Ranger
# versions used for the paper (6.x-7.x) write it unless `--no-loom` is passed, so
# this script does not trade the matrix for speed. If your version no longer emits
# a loom, convert `outs/raw_feature_bc_matrix.h5` with `scanpy.read_10x_h5` before
# running 01 -- do not switch the downstream script to H5 silently.
#
# No BAM is written. The manuscript reports clustering, within-population DE,
# GO/KEGG and CellChat, none of which need one. Add `--create-bam=true` only if a
# BAM is required for another purpose.
# =============================================================================
set -euo pipefail

REF=${1:?usage: 00_cellranger_count.sh <transcriptome_ref> <out_dir>}
OUT=${2:?usage: 00_cellranger_count.sh <transcriptome_ref> <out_dir>}
LIBS=${LIBS:-"NS1 YY1"}
CORES=${CORES:-32}
MEM=${MEM:-128}

mkdir -p "${OUT}"

for lib in ${LIBS}; do
    echo "[10x] cellranger count ${lib}"
    cellranger count \
        --id="${lib}" \
        --transcriptome="${REF}" \
        --fastqs="${OUT}/../fastq/${lib}" \
        --sample="${lib}" \
        --localcores="${CORES}" \
        --localmem="${MEM}" \
        --output-dir="${OUT}/cellranger"
done

cat <<'EOF'
[10x] next: point 01_qc_integrate_annotate.py at the raw matrices
      cellranger/NS1/outs/*_raw.loom  ->  --ctrl-loom  (--ctrl-label WT)
      cellranger/YY1/outs/*_raw.loom  ->  --case-loom  (--case-label KO)
EOF
echo "[10x] done"
