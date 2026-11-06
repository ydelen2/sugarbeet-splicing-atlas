#!/bin/bash
#SBATCH --job-name=bv_suppa_diff
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# SUPPA2 differential splicing for comparisons added after 07_suppa has run: reuses the existing
# filtered GTF, event definitions and PSI matrices and only writes the group files and diffSplice
# results of the comparisons listed in $COMPARISONS (give a subset table through COMPARISONS_FILE).
# The expression filter of 07_suppa depends on the comparison groups; the job first reports how many
# transcripts the filter would add with the current full comparison table, so that the decision to
# keep the existing event universe is documented in the log.
# Usage (from $NRDSTOR/bvas):
#   COMPARISONS_FILE=metadata/comparisons_new.tsv ANALYSIS=full sbatch workflow/07b_suppa_diff.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
PREP="python ${SCRIPTS_DIR}/analysis/prepare_splicing_inputs.py"
EV="${SUPPA_DIR}/events"; PSI="${SUPPA_DIR}/psi"; GRP="${SUPPA_DIR}/groups"; DIFF="${SUPPA_DIR}/diff"
for f in "$FINAL_GTF" "$EV/all_events.ioe" "$EV/merged.ioi" "$PSI/all_events.psi" "$SUPPA_DIR/tpm_all.tsv"; do
  [ -f "$f" ] || { echo "missing $f (07_suppa writes it)" >&2; exit 1; }
done
ISO_PSI=$(ls "$PSI"/transcripts*isoform.psi | awk 'NR==1')
FULL_COMPARISONS="${META_DIR}/comparisons_all.tsv"
echo "comparisons in this run: $(awk 'END{print NR-1}' "$COMPARISONS") from $COMPARISONS"

echo "== 0. effect of the added comparisons on the expression filter (report only) =="
TMP=$(job_tmp 5 "suppa_diff"); trap 'rm -rf "$TMP"' EXIT
python "${SCRIPTS_DIR}/analysis/filter_transcripts.py" expression --gtf "${FILT_DIR}/structural.gtf" \
  --map "${FILT_DIR}/transcript_map.tsv" --tpm "$SUPPA_DIR/tpm_all.tsv" \
  --min-tpm "$MIN_TPM" --min-runs "$MIN_RUNS" --criterion "$EXPR_FILTER" \
  --sheet "$SAMPLES" --comparisons "$FULL_COMPARISONS" --out "$TMP/final_full.gtf" > "$TMP/filter.log" 2>&1
tx() { awk -F'\t' '$3=="transcript"' "$1" | wc -l; }
echo "transcripts in the GTF used for the events: $(tx "$FINAL_GTF"); with the current full comparison table: $(tx "$TMP/final_full.gtf")"
awk -F'\t' '$3=="transcript" {match($9, /transcript_id "[^"]+"/); print substr($9, RSTART+15, RLENGTH-16)}' "$FINAL_GTF" | sort > "$TMP/a.txt"
awk -F'\t' '$3=="transcript" {match($9, /transcript_id "[^"]+"/); print substr($9, RSTART+15, RLENGTH-16)}' "$TMP/final_full.gtf" | sort > "$TMP/b.txt"
echo "transcripts only in the current-table filter: $(comm -13 "$TMP/a.txt" "$TMP/b.txt" | wc -l); only in the event GTF: $(comm -23 "$TMP/a.txt" "$TMP/b.txt" | wc -l)"

echo "== 1. group files =="
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$SUPPA_DIR/tpm_all.tsv" --suffix tpm --out-dir "$GRP/events"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$PSI/all_events.psi" --suffix psi --out-dir "$GRP/events"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$SUPPA_DIR/tpm_all.tsv" --suffix tpm --out-dir "$GRP/isoforms"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$ISO_PSI" --suffix psi --out-dir "$GRP/isoforms"

echo "== 2. differential splicing (empirical method, gene-level multiple testing correction) =="
while IFS=$'\t' read -r -u 3 CID REF TEST; do
  mkdir -p "$DIFF/$CID"
  suppa.py diffSplice -m empirical -gc -i "$EV/all_events.ioe" \
    -p "$GRP/events/$CID/$REF.psi" "$GRP/events/$CID/$TEST.psi" \
    -e "$GRP/events/$CID/$REF.tpm" "$GRP/events/$CID/$TEST.tpm" \
    -o "$DIFF/$CID/events"
  suppa.py diffSplice -m empirical -gc -i "$EV/merged.ioi" \
    -p "$GRP/isoforms/$CID/$REF.psi" "$GRP/isoforms/$CID/$TEST.psi" \
    -e "$GRP/isoforms/$CID/$REF.tpm" "$GRP/isoforms/$CID/$TEST.tpm" \
    -o "$DIFF/$CID/isoforms"
  awk -F'\t' -v c="$CID" 'NR==1 {print c ": dpsi header: " $0; next}
       {n++} $2!="nan" && $3!="nan" {t++} $2!="nan" && $3!="nan" && ($2>=0.1 || $2<=-0.1) && $3<0.05 {s++}
       END {printf "%s: events %d, testable %d, significant (|dPSI|>=0.1, p<0.05) %d\n", c, n, t, s}' "$DIFF/$CID/events.dpsi"
done 3< <(awk -F'\t' 'NR>1 {print $1 "\t" $5 "\t" $6}' "$COMPARISONS")
echo "SUPPA2 diffSplice done for the listed comparisons: $DIFF"
