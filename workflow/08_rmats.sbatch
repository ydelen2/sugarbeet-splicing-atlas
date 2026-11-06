#!/bin/bash
#SBATCH --job-name=bv_rmats
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48G
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
# Junction-count based differential splicing with rMATS, one array task per comparison.
# Uses the same filtered GTF as SUPPA2, so run it after 07_suppa has written $FINAL_GTF.
# Also needs the BAM files and the strandedness calls from 04a.
# Usage (from $NRDSTOR/bvas):
#   N=$(($(wc -l < metadata/comparisons_pilot.tsv)-1)); sbatch --array=1-${N} workflow/08_rmats.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-8}

CID=$(awk -F'\t' -v L="${SLURM_ARRAY_TASK_ID}" 'NR==L+1 {print $1}' "$COMPARISONS")
[ -n "$CID" ] || { echo "no comparison on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
OUT="${RMATS_DIR}/${CID}"
[ -f "$FINAL_GTF" ] || { echo "final GTF not found: $FINAL_GTF (07_suppa writes it)" >&2; exit 1; }
TMP=$(job_tmp 50 "rmats_${CID}"); trap 'rm -rf "$TMP"' EXIT

# Input lists are written to the job's own temp directory, so parallel array tasks never share files
IN="$TMP/inputs"
python "${SCRIPTS_DIR}/analysis/prepare_splicing_inputs.py" rmats --sheet "$SAMPLES" --comparisons "$COMPARISONS" \
  --bam-dir "$BAM_DIR" --runres-dir "$RUNRES_DIR" --out-dir "$IN" --only "$CID"
# shellcheck disable=SC1090
source "$IN/$CID/params.txt"

# STAR soft-clips read ends; rMATS skips clipped reads unless --allow-clipping is given (if supported).
# The help text is captured first: piping it into "grep -q" would fail under "set -o pipefail".
EXTRA=""
HELP=$(rmats.py --help 2>&1 || true)
case "$HELP" in *--allow-clipping*) EXTRA="--allow-clipping";; esac
NOVEL=""; [ "${RMATS_NOVELSS}" = "1" ] && NOVEL="--novelSS"
echo "rMATS options: readLength=$READ_LEN libType=$LIBTYPE extra='$EXTRA' novel='$NOVEL' out=$OUT"
rmats.py --b1 "$IN/$CID/b1.txt" --b2 "$IN/$CID/b2.txt" --gtf "$FINAL_GTF" \
  -t paired --readLength "$READ_LEN" --variable-read-length --libType "$LIBTYPE" \
  --nthread "$T" --od "$TMP/out" --tmp "$TMP/tmp" $EXTRA $NOVEL

mkdir -p "$OUT"
cp "$TMP"/out/*.txt "$OUT/"
cp "$IN/$CID/b1.txt" "$IN/$CID/b2.txt" "$IN/$CID/params.txt" "$OUT/"
echo "comparison $CID: sample_1 = $B1_LABEL, sample_2 = $B2_LABEL (IncLevelDifference = sample_1 minus sample_2)"
for f in "$OUT"/*.MATS.JCEC.txt; do
  awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; next}
              $c["FDR"]!="NA" && $c["FDR"]<0.05 && ($c["IncLevelDifference"]>=0.1 || $c["IncLevelDifference"]<=-0.1) {n++}
              END {printf "%s significant (FDR<0.05, |dPSI|>=0.1): %d\n", FILENAME, n}' "$f"
done
