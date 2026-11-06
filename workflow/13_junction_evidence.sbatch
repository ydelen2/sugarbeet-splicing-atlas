#!/bin/bash
#SBATCH --job-name=bv_junction_evidence
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Independent evidence for the junctions and splicing events:
#  1. the two public PacBio full-length cDNA runs of the reference genotype KWS2320 (metadata/longread_runs.tsv)
#     are downloaded, checked against their ENA MD5 sums, aligned with minimap2 in spliced mode and compared
#     with the annotated, pooled and expressed-set junctions and with the high-confidence events
#     (scripts/analysis/longread_junctions.py);
#  2. the junctions of the expressed transcript set that are neither annotated nor pooled are looked up in the
#     STAR SJ.out.tab of every aligned run (scripts/analysis/junction_support.py).
# Needs minimap2 in the project environment (install once):
#   conda install -p "$ENV_DIR" -c conda-forge -c bioconda minimap2
# Usage (from $NRDSTOR/bvas):
#   ANALYSIS=full sbatch workflow/13_junction_evidence.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-16}
command -v minimap2 >/dev/null || { echo "minimap2 not found in $ENV_DIR (see the header of this script)" >&2; exit 1; }
OUT="${RES_DIR}/junction_evidence"; LR="${BULK}/longread"; mkdir -p "$OUT" "$LR"
echo "minimap2 $(minimap2 --version); $(samtools --version | head -n 1)"

echo "== 1. long-read download =="
tail -n +2 "${META_DIR}/longread_runs.tsv" | while IFS=$'\t' read -r run study inst strat lib nreads nbases url md5 bytes; do
  dest="$LR/${run}.fastq.gz"
  if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then echo "ok (cached) $run"; continue; fi
  ok=0
  for attempt in 1 2 3 4; do
    rm -f "$dest" "$dest.aria2"
    aria2c -q -x 2 -s 2 --file-allocation=none --max-tries=3 --retry-wait=60 --connect-timeout=60 -d "$LR" -o "${run}.fastq.gz" "$url" || true
    if [ -s "$dest" ] && [ ! -e "$dest.aria2" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then ok=1; break; fi
    sleep $((attempt * 120))
  done
  [ "$ok" = 1 ] || { echo "download failed: $run" >&2; exit 1; }
  echo "downloaded $run ($(stat -c %s "$dest") bytes, MD5 ok)"
done

echo "== 2. spliced alignment and comparison =="
EVENTS="--events main=${RES_DIR}/summary/events_high_confidence_annotated.tsv"
[ -s "${RES_DIR}/summary_tissue/events_high_confidence_annotated.tsv" ] && \
  EVENTS="$EVENTS --events tissue=${RES_DIR}/summary_tissue/events_high_confidence_annotated.tsv"
# -ub: the cDNA reads are not stranded; -G 100k: the maximum intron length used for STAR
zcat "$LR"/*.fastq.gz \
  | minimap2 -ax splice -ub -G 100k --secondary=no -t "$T" "$GENOME_FA" - 2> "$OUT/minimap2.log" \
  | samtools view -F 0x904 -q 1 - \
  | python "${SCRIPTS_DIR}/analysis/longread_junctions.py" \
      --ref-junctions "${MERGE_DIR}/ref_junctions.tsv" --pooled-junctions "${MERGE_DIR}/pooled_junctions.tsv" \
      --gtf "$FINAL_GTF" --ioe "${SUPPA_DIR}/events/all_events.ioe" $EVENTS --out-prefix "$OUT/longread"

echo "== 3. Illumina support of the expressed-set junctions outside the pooled table =="
python "${SCRIPTS_DIR}/analysis/junction_support.py" --gtf "$FINAL_GTF" \
  --ref-junctions "${MERGE_DIR}/ref_junctions.tsv" --pooled-junctions "${MERGE_DIR}/pooled_junctions.tsv" \
  --runres-dir "$RUNRES_DIR" --sheet "$SAMPLES" --out-prefix "$OUT/support"
echo "done: $OUT"
