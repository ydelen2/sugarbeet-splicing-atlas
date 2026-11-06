#!/bin/bash
#SBATCH --job-name=bv_salmon
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=24G
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
# Transcript quantification against the structurally filtered transcriptome (run after 04c_filter_index).
# Usage (from $NRDSTOR/bvas):
#   N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%20 workflow/06_salmon_quant.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-8}

i=${SLURM_ARRAY_TASK_ID}
skip_if_excluded "$i"
RUN=$(sample_field "$i" run)
IN="${FASTQ_DIR}/${RUN}"
OUT="${SALMON_DIR}/${RUN}"
IDX="${FILT_DIR}/salmon_index"
[ -d "$IDX" ] || { echo "salmon index not found: $IDX (run 04c_filter_index first)" >&2; exit 1; }
[ -f "$OUT/.salmon_ok" ] && { echo "already quantified $RUN"; exit 0; }
[ -f "$IN/.download_ok" ] && [ -s "$IN/${RUN}_1.fastq.gz" ] && [ -s "$IN/${RUN}_2.fastq.gz" ] \
  || { echo "FASTQ not available for $RUN (deleted after an earlier quantification? download it again)" >&2; exit 1; }

# Library type: detected automatically (-l A) unless metadata/curation/salmon_libtype.tsv sets it for
# this run; the detected type is written to lib_format_counts.json
LIBTYPE="A"
LT_FILE="${META_DIR}/curation/salmon_libtype.tsv"
if [ -f "$LT_FILE" ]; then
  OVR=$(awk -F'\t' -v r="$RUN" 'NR>1 && $1==r {print $2}' "$LT_FILE")
  [ -n "$OVR" ] && LIBTYPE="$OVR" && echo "library type set to $LIBTYPE by $LT_FILE"
fi

TMP=$(job_tmp 10 "salmon_${RUN}"); trap 'rm -rf "$TMP"' EXIT
salmon quant -i "$IDX" -l "$LIBTYPE" \
  -1 "$IN/${RUN}_1.fastq.gz" -2 "$IN/${RUN}_2.fastq.gz" \
  -p "$T" --gcBias --seqBias \
  -o "$TMP/out"
rm -rf "$OUT"; mkdir -p "$OUT"
cp -r "$TMP/out/." "$OUT/"
grep -E '"percent_mapped"|"num_processed"' "$OUT/aux_info/meta_info.json" || true
grep -E '"expected_format"|"compatible_fragment_ratio"' "$OUT/lib_format_counts.json" || true
touch "$OUT/.salmon_ok"

if [ "${DELETE_FASTQ_AFTER_QUANT}" = "1" ] && [ -f "${RUNRES_DIR}/${RUN}/.align_ok" ]; then
  rm -f "$IN"/*.fastq.gz
  echo "FASTQ deleted for $RUN"
fi
