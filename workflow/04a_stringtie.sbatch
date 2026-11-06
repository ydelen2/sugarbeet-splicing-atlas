#!/bin/bash
#SBATCH --job-name=bv_stringtie
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
# Usage (from $NRDSTOR/bvas):
#   N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%30 workflow/04a_stringtie.sbatch
# Reference-guided assembly per run. Strandedness inferred from STAR ReadsPerGene (col3 vs col4).
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
i=${SLURM_ARRAY_TASK_ID}
skip_if_excluded "$i"
RUN=$(sample_field "$i" run)
BAM="${BAM_DIR}/${RUN}/${RUN}.Aligned.sortedByCoord.out.bam"
KEEP="${RUNRES_DIR}/${RUN}"
OUT="${ST_DIR}/${RUN}"; mkdir -p "$OUT"
[ -f "$KEEP/.align_ok" ] && [ -f "$BAM" ] || { echo "no alignment for $RUN" >&2; exit 1; }

RPG="$KEEP/${RUN}.ReadsPerGene.out.tab"
read fwd rev < <(tail -n +5 "$RPG" | awk '{f+=$3; r+=$4} END{print f+0, r+0}')
if   [ "$fwd" -gt $((rev*4)) ]; then LIB="--fr"; STR="forward"
elif [ "$rev" -gt $((fwd*4)) ]; then LIB="--rf"; STR="reverse"
else LIB=""; STR="unstranded"; fi
echo -e "${RUN}\t${STR}\t${fwd}\t${rev}" > "$KEEP/strandedness.tsv"

TMP=$(job_tmp 5 "st_${RUN}"); trap 'rm -rf "$TMP"' EXIT
stringtie "$BAM" -p 4 $LIB -G "$GENES_GTF" -o "$TMP/${RUN}.gtf" -A "$TMP/${RUN}.gene_abund.tab" \
  -f 0.05 -m 200 -c 2 -j 2
cp "$TMP/${RUN}.gtf" "$TMP/${RUN}.gene_abund.tab" "$OUT/"
touch "$OUT/.stringtie_ok"
