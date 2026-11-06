#!/bin/bash
#SBATCH --job-name=bv_suppa
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Splicing events, PSI values and differential splicing with SUPPA2 (run after 06_salmon_quant).
# Usage (from $NRDSTOR/bvas):  sbatch workflow/07_suppa.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
PREP="python ${SCRIPTS_DIR}/analysis/prepare_splicing_inputs.py"
GTF="$FINAL_GTF"
EV="${SUPPA_DIR}/events"; PSI="${SUPPA_DIR}/psi"; GRP="${SUPPA_DIR}/groups"; DIFF="${SUPPA_DIR}/diff"
mkdir -p "$EV" "$PSI" "$GRP" "$DIFF"
[ -f "${FILT_DIR}/structural.gtf" ] || { echo "run 04c_filter_index first" >&2; exit 1; }

echo "== 0. transcript TPM matrix and expression filter =="
$PREP tpm --sheet "$SAMPLES" --salmon-dir "$SALMON_DIR" --out "$SUPPA_DIR/tpm_all.tsv"
python "${SCRIPTS_DIR}/analysis/filter_transcripts.py" expression --gtf "${FILT_DIR}/structural.gtf" \
  --map "${FILT_DIR}/transcript_map.tsv" --tpm "$SUPPA_DIR/tpm_all.tsv" \
  --min-tpm "$MIN_TPM" --min-runs "$MIN_RUNS" --criterion "$EXPR_FILTER" \
  --sheet "$SAMPLES" --comparisons "$COMPARISONS" --out "$GTF"

echo "== 1. events from the filtered transcriptome =="
# Local events: SE, alternative 5' and 3' sites (SS), MX, RI, alternative first and last exons (FL)
suppa.py generateEvents -i "$GTF" -o "$EV/merged" -f ioe -e SE SS MX RI FL
suppa.py generateEvents -i "$GTF" -o "$EV/merged" -f ioi
awk 'FNR==1 && NR!=1 {next} {print}' "$EV"/merged_*_strict.ioe > "$EV/all_events.ioe"
for f in "$EV"/merged_*_strict.ioe; do echo "$(basename "$f") $(awk 'END{print NR-1}' "$f")"; done
echo "all events: $(awk 'END{print NR-1}' "$EV/all_events.ioe")"

echo "== 2. PSI per event and per isoform =="
# The TPM matrix still holds the transcripts removed by the expression filter; SUPPA2 only uses the
# transcripts listed in the event files, which come from the filtered GTF.
suppa.py psiPerEvent -i "$EV/all_events.ioe" -e "$SUPPA_DIR/tpm_all.tsv" -o "$PSI/all_events"
suppa.py psiPerIsoform -g "$GTF" -e "$SUPPA_DIR/tpm_all.tsv" -o "$PSI/transcripts"
ISO_PSI=$(ls "$PSI"/transcripts*isoform.psi | awk 'NR==1')

echo "== 3. group files per comparison =="
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$SUPPA_DIR/tpm_all.tsv" --suffix tpm --out-dir "$GRP/events"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$PSI/all_events.psi" --suffix psi --out-dir "$GRP/events"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$SUPPA_DIR/tpm_all.tsv" --suffix tpm --out-dir "$GRP/isoforms"
$PREP split --sheet "$SAMPLES" --comparisons "$COMPARISONS" --matrix "$ISO_PSI" --suffix psi --out-dir "$GRP/isoforms"

echo "== 4. differential splicing (empirical method, gene-level multiple testing correction) =="
# Conditions are passed in the order reference, test; the header of each .dpsi file names the
# two conditions and shows the direction of dPSI.
# The comparison list is read through file descriptor 3, so commands inside the loop cannot consume it
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
  # events with |dPSI| >= 0.1 and p < 0.05
  awk -F'\t' -v c="$CID" 'NR==1 {print c ": dpsi header: " $0; next}
       {n++} $2!="nan" && $3!="nan" {t++} $2!="nan" && $3!="nan" && ($2>=0.1 || $2<=-0.1) && $3<0.05 {s++}
       END {printf "%s: events %d, testable %d, significant (|dPSI|>=0.1, p<0.05) %d\n", c, n, t, s}' "$DIFF/$CID/events.dpsi"
done 3< <(awk -F'\t' 'NR>1 {print $1 "\t" $5 "\t" $6}' "$COMPARISONS")
echo "SUPPA2 done: $SUPPA_DIR"
