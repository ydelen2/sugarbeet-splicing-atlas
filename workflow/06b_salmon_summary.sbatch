#!/bin/bash
#SBATCH --job-name=bv_salmon_sum
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=48G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# After 06_salmon_quant: transcript TPM matrix, Salmon QC table, and a scan of expression-filter
# thresholds used to choose MIN_TPM and MIN_RUNS before 07_suppa.
# Usage (from $NRDSTOR/bvas):  ANALYSIS=full sbatch workflow/06b_salmon_summary.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
SUM="${RES_DIR}/summary"; mkdir -p "$SUM" "$SUPPA_DIR"
A="${SCRIPTS_DIR}/analysis"

echo "== Salmon QC =="
python "$A/salmon_qc.py" --sheet "$SAMPLES" --salmon-dir "$SALMON_DIR" --runres-dir "$RUNRES_DIR" \
  --out "$SUM/salmon_qc.tsv"

echo; echo "== transcript TPM matrix =="
# the sheet given to the matrix builder excludes use_das=excluded runs, which have no Salmon output
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="use_das") u=i; print; next} $u!="excluded"' "$SAMPLES" \
  > "$SUM/samples_quantified.tsv"
python "$A/prepare_splicing_inputs.py" tpm --sheet "$SUM/samples_quantified.tsv" --salmon-dir "$SALMON_DIR" \
  --out "$SUPPA_DIR/tpm_all.tsv"

echo; echo "== expression filter scan =="
python "$A/expression_filter_scan.py" --tpm "$SUPPA_DIR/tpm_all.tsv" --map "${FILT_DIR}/transcript_map.tsv" \
  --sheet "$SAMPLES" --comparisons "$COMPARISONS" --out "$SUM/expression_filter_scan.tsv"
echo "tables: $SUM/salmon_qc.tsv  $SUM/expression_filter_scan.tsv"
