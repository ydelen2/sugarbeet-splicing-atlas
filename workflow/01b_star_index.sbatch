#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Usage (from $NRDSTOR/bvas):  sbatch workflow/01b_star_index.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env

# Build in temp, then copy to NRDStor (index is read-only afterwards)
TMP=$(job_tmp 40 star_index); trap 'rm -rf "$TMP"' EXIT
# genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1); 6.5e8 bp -> 13
STAR --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir "$TMP/idx" \
  --outTmpDir "$TMP/star_tmp" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GENES_GTF" \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 13
rm -rf "$STAR_INDEX"; mkdir -p "$STAR_INDEX"
cp -r "$TMP/idx/." "$STAR_INDEX/"
du -sh "$STAR_INDEX"
echo "STAR index done: $STAR_INDEX"
