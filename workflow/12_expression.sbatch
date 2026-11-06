#!/bin/bash
#SBATCH --job-name=bv_expression
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Differential gene expression per comparison from the STAR gene counts (PyDESeq2) and its overlap with the
# high-confidence differential splicing events. Needs the outputs of 03_align (counts, strandedness calls of
# 04a) and 09_summarize (high-confidence events, splicing-regulator families).
# PyDESeq2 is installed once into the project environment with
#   conda install -p "$ENV_DIR" -c conda-forge -c bioconda pydeseq2
# Usage (from $NRDSTOR/bvas):
#   ANALYSIS=full sbatch workflow/12_expression.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-8}
SUMMARY_DIR="${SUMMARY_DIR:-${RES_DIR}/summary}"
OUT="${RES_DIR}/expression"; mkdir -p "$OUT"
python -c "import pydeseq2; print('pydeseq2', pydeseq2.__version__)"
FAM="${RES_DIR}/domains/splicing_regulators_families.tsv"
FAM_OPT=""; [ -s "$FAM" ] && FAM_OPT="--families $FAM"
python "${SCRIPTS_DIR}/analysis/differential_expression.py" --sheet "$SAMPLES" --comparisons "$COMPARISONS" \
  --runres-dir "$RUNRES_DIR" --events "${SUMMARY_DIR}/event_concordance_high_confidence.tsv" \
  --suppa-dir "$SUPPA_DIR" --map "${FILT_DIR}/transcript_map.tsv" $FAM_OPT \
  --out-prefix "${OUT}/de" --threads "$T"
echo "done: ${OUT}/de_summary.tsv"
