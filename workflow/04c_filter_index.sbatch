#!/bin/bash
#SBATCH --job-name=bv_filter
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Structural filter of the merged transcriptome (class codes, exon number) and Salmon index on the result.
# Run after 04b_merge. Usage (from $NRDSTOR/bvas):  sbatch workflow/04c_filter_index.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
mkdir -p "$FILT_DIR"
TMP=$(job_tmp 20 filter); trap 'rm -rf "$TMP"' EXIT

python "${SCRIPTS_DIR}/analysis/filter_transcripts.py" structural \
  --merged "${MERGE_DIR}/merged.gtf" --annotated "${MERGE_DIR}/gffcmp.annotated.gtf" --out-dir "$TMP/f" \
  --drop-seq-prefix "$ORGANELLE_PREFIX"

cd "$TMP/f"
gffread -w structural.transcripts.fa -g "$GENOME_FA" structural.gtf
echo "transcript sequences: $(grep -c '^>' structural.transcripts.fa)"
awk '/^>/ {print substr($1, 2)}' "$GENOME_FA" > decoys.txt
cat structural.transcripts.fa "$GENOME_FA" > gentrome.fa
# --keepDuplicates: by default Salmon drops transcripts whose sequence is identical to another one
# (on chromosomes and unplaced scaffolds). They would then have no expression value
# and SUPPA2 would compute PSI for their genes without them.
salmon index -t gentrome.fa -d decoys.txt -i salmon_index -k 31 -p 8 --keepDuplicates > salmon_index.log 2>&1
rm -f gentrome.fa
grep -i "duplicate" salmon_index.log || echo "no duplicate warning in the index log"
awk '/^>/ {n++} END {print "index built; decoy sequences: " n}' "$GENOME_FA"

rm -rf "${FILT_DIR}/salmon_index"
cp -r structural.gtf structural_report.tsv transcript_map.tsv structural.transcripts.fa decoys.txt \
      salmon_index salmon_index.log "$FILT_DIR/"
echo "filter done -> $FILT_DIR"
