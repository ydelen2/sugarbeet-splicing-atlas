#!/bin/bash
#SBATCH --job-name=bv_align
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
# Usage (from $NRDSTOR/bvas):
#   N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%20 workflow/03_align.sbatch
# Per run: fastp -> STAR 2-pass -> sorted BAM (bulk storage) + small outputs copied to NRDStor.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-8}

i=${SLURM_ARRAY_TASK_ID}
skip_if_excluded "$i"
RUN=$(sample_field "$i" run)
IN="${FASTQ_DIR}/${RUN}"
OUT="${BAM_DIR}/${RUN}"                 # bulk: BAM
KEEP="${RUNRES_DIR}/${RUN}"             # persistent: junctions, counts, logs
mkdir -p "$OUT" "$KEEP" "$QC_DIR/fastp"
[ -f "$IN/.download_ok" ] || { echo "FASTQ not ready for $RUN" >&2; exit 1; }
[ -f "$KEEP/.align_ok" ] && [ -f "$OUT/${RUN}.Aligned.sortedByCoord.out.bam" ] && { echo "already aligned $RUN"; exit 0; }

# Temp space scales with the run: trimmed FASTQ, STAR sorting bins and the BAM are all written there.
# A 34 Gb run produces about 25 GiB of trimmed FASTQ and a 20 GiB BAM. Four times the sequenced Gb plus
# 20 GB is a generous bound; if the node-local disk is smaller, job_tmp falls back to network scratch.
NEED_GB=$(awk -v g="$(sample_field "$i" gb)" 'BEGIN {printf "%d", g*4 + 20}')
TMP=$(job_tmp "$NEED_GB" "$RUN"); trap 'rm -rf "$TMP"' EXIT
echo "temp space requested: ${NEED_GB} GB"
echo "tmp: $TMP"

# 1) fastp
fastp -w "$T" \
  -i "$IN/${RUN}_1.fastq.gz" -I "$IN/${RUN}_2.fastq.gz" \
  -o "$TMP/${RUN}_1.fq.gz" -O "$TMP/${RUN}_2.fq.gz" \
  --detect_adapter_for_pe --trim_poly_g --cut_right --cut_mean_quality 20 --length_required 50 \
  -j "$QC_DIR/fastp/${RUN}.fastp.json" -h "$QC_DIR/fastp/${RUN}.fastp.html" 2> "$KEEP/fastp.log"

# 2) STAR 2-pass, written to local tmp first, then moved (avoids heavy random I/O on network storage)
STAR --runThreadN "$T" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$TMP/${RUN}_1.fq.gz" "$TMP/${RUN}_2.fq.gz" --readFilesCommand zcat \
  --outFileNamePrefix "$TMP/${RUN}." \
  --outTmpDir "$TMP/star_tmp" \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM NM MD jM jI XS \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 \
  --alignIntronMin 20 --alignIntronMax 100000 --alignMatesGapMax 100000 \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSJtype Standard \
  --quantMode GeneCounts \
  --outBAMsortingThreadN 4 --limitBAMsortRAM 40000000000

samtools index -@ "$T" "$TMP/${RUN}.Aligned.sortedByCoord.out.bam"
samtools flagstat -@ "$T" "$TMP/${RUN}.Aligned.sortedByCoord.out.bam" > "$TMP/${RUN}.flagstat"

# 3) persist small outputs on NRDStor, BAM to bulk storage
cp "$TMP/${RUN}.SJ.out.tab" "$TMP/${RUN}.ReadsPerGene.out.tab" "$TMP/${RUN}.Log.final.out" \
   "$TMP/${RUN}.Log.out" "$TMP/${RUN}.flagstat" "$KEEP/"
cp "$TMP/${RUN}._STARpass1/SJ.out.tab" "$KEEP/${RUN}.SJ.pass1.out.tab" 2>/dev/null || true
mv "$TMP/${RUN}.Aligned.sortedByCoord.out.bam" "$TMP/${RUN}.Aligned.sortedByCoord.out.bam.bai" "$OUT/"

grep -E "Uniquely mapped reads %|Number of splices: Total|Mismatch rate" "$KEEP/${RUN}.Log.final.out"
touch "$KEEP/.align_ok"
# FASTQ files are kept here; 06_salmon_quant reads them again and deletes them if DELETE_FASTQ_AFTER_QUANT=1
