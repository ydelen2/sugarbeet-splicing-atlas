#!/bin/bash
# Reference for the full analysis: RefBeet-3.0 plus the sugar beet chloroplast and mitochondrial genomes.
# RefBeet-3.0 has no organelle sequences, so organelle reads could otherwise be placed on nuclear
# organelle-derived insertions. The organelle sequences carry no gene models; transcripts assembled on
# them are removed by the structural filter (sequence names start with $ORGANELLE_PREFIX).
#   chloroplast   GenBank KR230391.1  (Beta vulgaris subsp. vulgaris cultivar KWS2320, 149,722 bp)
#   mitochondrion RefSeq  NC_002511.2 (Beta vulgaris subsp. vulgaris, 368,801 bp)
# Run on a login node from $NRDSTOR/bvas:  bash workflow/01c_add_organelles.sh
# Then:  ANALYSIS=full sbatch workflow/01b_star_index.sbatch
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export ANALYSIS=full
source "$HERE/config.sh"
activate_env
cd "$REF_DIR"
[ -s "$NUCLEAR_FA" ] || { echo "run 01a_get_reference.sh first" >&2; exit 1; }

EF="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
[ -s organelle_cp.raw.fa ] || curl -sf "${EF}KR230391.1"  -o organelle_cp.raw.fa
sleep 1
[ -s organelle_mt.raw.fa ] || curl -sf "${EF}NC_002511.2" -o organelle_mt.raw.fa

check() { # file expected_length label
  local n len
  n=$(grep -c '^>' "$1"); len=$(awk '!/^>/ {s+=length($0)} END {print s}' "$1")
  echo "$3: $n sequence(s), $len bp; header: $(awk 'NR==1' "$1")"
  [ "$n" = "1" ] && [ "$len" = "$2" ] || { echo "unexpected content in $1" >&2; exit 1; }
}
check organelle_cp.raw.fa 149722 chloroplast
check organelle_mt.raw.fa 368801 mitochondrion

{ cat "$NUCLEAR_FA"
  awk -v p="${ORGANELLE_PREFIX}cp_" '/^>/ {print ">" p substr($1, 2); next} {print}' organelle_cp.raw.fa
  awk -v p="${ORGANELLE_PREFIX}mt_" '/^>/ {print ">" p substr($1, 2); next} {print}' organelle_mt.raw.fa
} > "$GENOME_FA"
samtools faidx "$GENOME_FA"
awk -F'\t' -v p="$ORGANELLE_PREFIX" '{n++; s+=$2} index($1, p)==1 {print "added: " $1 "\t" $2}
     END {print "sequences: " n "  total: " s " bp"}' "${GENOME_FA}.fai"
md5sum organelle_cp.raw.fa organelle_mt.raw.fa "$(basename "$GENOME_FA")" > organelles.md5
echo "Next: ANALYSIS=full sbatch workflow/01b_star_index.sbatch"
