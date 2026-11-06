#!/bin/bash
#SBATCH --job-name=bv_origin
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Where do reads that STAR could not place on RefBeet-3.0 come from?
# A subsample of read pairs from selected runs is aligned with HISAT2 to a small index made of
#   cp_  sugar beet chloroplast, KWS2320 (GenBank KR230391.1)
#   mt_  sugar beet mitochondrion (RefSeq NC_002511.2)
#   ssc_ Sclerotinia sclerotiorum strain 1980 (RefSeq assembly GCF_000146945.2)
# and the share of reads assigned to each source is reported. The script also lists RefBeet-3.0
# sequences whose size or name suggests an organelle genome.
# Usage (from $NRDSTOR/bvas):  sbatch workflow/diagnostics/read_origin.sbatch RUN1 RUN2 ...
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
T=${SLURM_NTASKS_PER_NODE:-8}
NPAIRS=200000
[ "$#" -gt 0 ] || { echo "give run accessions as arguments" >&2; exit 1; }

D="${RES_DIR}/diagnostics/read_origin"; mkdir -p "$D/idx" "$D/sam"; cd "$D"

echo "== RefBeet-3.0 sequences that could be organelle genomes =="
awk -F'\t' 'tolower($1) ~ /chl|plast|mito|^cp|^mt|chrc|chrm/ || ($2 > 100000 && $2 < 500000) {print $1 "\t" $2}' \
  "${GENOME_FA}.fai" > refbeet_candidate_organelle.tsv
echo "candidates (name hint or 100-500 kb): $(awk 'END{print NR}' refbeet_candidate_organelle.tsv)"
awk -F'\t' 'tolower($1) ~ /chl|plast|mito|^cp|^mt|chrc|chrm/' refbeet_candidate_organelle.tsv

EF="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
SSC="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/945/GCF_000146945.2_ASM14694v2/GCF_000146945.2_ASM14694v2_genomic.fna.gz"
[ -s cp.fa ]  || curl -sf "${EF}KR230391.1"  -o cp.fa
[ -s mt.fa ]  || curl -sf "${EF}NC_002511.2" -o mt.fa
[ -s ssc.fa ] || curl -sf "$SSC" | gunzip -c > ssc.fa
for f in cp mt ssc; do
  printf "%s: %s sequences, %s bp, first header: %s\n" "$f" "$(grep -c '^>' $f.fa)" \
    "$(awk '!/^>/{n+=length($0)} END{print n}' $f.fa)" "$(awk '/^>/{print; exit}' $f.fa)"
done
if [ ! -s idx/ref.1.ht2 ]; then
  { awk '/^>/{print ">cp_" substr($1,2); next} {print}' cp.fa
    awk '/^>/{print ">mt_" substr($1,2); next} {print}' mt.fa
    awk '/^>/{print ">ssc_" substr($1,2); next} {print}' ssc.fa; } > ref.fa
  hisat2-build -p "$T" ref.fa idx/ref > idx/build.log 2>&1
fi

printf "run\tstudy\tcondition\ttissue\treads\tchloroplast_pct\tmitochondrion_pct\tsclerotinia_pct\tnot_in_these_pct\tstar_uniq_pct\n" > read_origin.tsv
for RUN in "$@"; do
  line=$(awk -F'\t' -v r="$RUN" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["run"]==r {print $c["study"] "\t" $c["condition"] "\t" $c["tissue"]}' "$SAMPLES")
  [ -n "$line" ] || { echo "$RUN not in $SAMPLES, skipped"; continue; }
  IN="${FASTQ_DIR}/${RUN}"
  [ -s "$IN/${RUN}_1.fastq.gz" ] || { echo "FASTQ missing for $RUN, skipped"; continue; }
  seqkit head -n "$NPAIRS" "$IN/${RUN}_1.fastq.gz" > sub_1.fq
  seqkit head -n "$NPAIRS" "$IN/${RUN}_2.fastq.gz" > sub_2.fq
  hisat2 -p "$T" -x idx/ref -1 sub_1.fq -2 sub_2.fq --no-unal -k 1 2> "sam/${RUN}.hisat2.log" \
    | samtools view -F 0x904 - | cut -f3 > "sam/${RUN}.targets"
  uniq=$(awk -F'|' '/Uniquely mapped reads %/{gsub(/[ \t%]/,"",$2); print $2}' "${RUNRES_DIR}/${RUN}/${RUN}.Log.final.out")
  awk -v run="$RUN" -v meta="$line" -v n=$((2*NPAIRS)) -v u="$uniq" '
      /^cp_/ {cp++} /^mt_/ {mt++} /^ssc_/ {ss++}
      END {printf "%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n", run, meta, n,
           100*cp/n, 100*mt/n, 100*ss/n, 100*(n-cp-mt-ss)/n, u}' "sam/${RUN}.targets" >> read_origin.tsv
done
rm -f sub_1.fq sub_2.fq
column -t -s $'\t' read_origin.tsv
echo "table: $D/read_origin.tsv"
