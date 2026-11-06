#!/bin/bash
#SBATCH --job-name=bv_domains
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Protein domain scan of the BeetSet-3 proteome (Pfam-A with the gathering thresholds, HMMER) and
# DIAMOND searches against the Arabidopsis thaliana UniProt reference proteome in both directions,
# the inputs of the splicing-regulator family analysis. Reference files are downloaded once into
# ref/pfam and ref/proteomes; the Pfam release and the download dates are written to the log.
# Usage (from $NRDSTOR/bvas):  ANALYSIS=full sbatch workflow/10_domains.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
module load compiler/gcc/11 openmpi/4.1          # toolchain the tool modules are built against (module spider hmmer/3.4)
module load hmmer/3.4 diamond/2.1
T="${SLURM_CPUS_ON_NODE:-16}"

PFAM_DIR="${REF_DIR}/pfam"; PROT_DIR="${REF_DIR}/proteomes"; OUT="${RES_DIR}/domains"
mkdir -p "$PFAM_DIR" "$PROT_DIR" "$OUT"
BEET_AA="${REF_DIR}/BeetSet-3.genes.2604.aa"
ATH_ID="UP000006548"                                   # Arabidopsis thaliana reference proteome
ATH_AA="${PROT_DIR}/athaliana_${ATH_ID}.fa"

echo "== versions =="; hmmsearch -h | grep -m1 'HMMER'; diamond version

echo; echo "== Pfam-A =="
if [ ! -s "$PFAM_DIR/Pfam-A.hmm" ]; then
  cd "$PFAM_DIR"
  curl -sSL -o Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
  curl -sSL -o Pfam.version.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
  gzip -t Pfam-A.hmm.gz && gunzip -f Pfam-A.hmm.gz && gunzip -f Pfam.version.gz
  date -u +"downloaded %Y-%m-%d" > download_date.txt
  cd "$PROJ"
fi
cat "$PFAM_DIR/Pfam.version" "$PFAM_DIR/download_date.txt"
[ -s "$PFAM_DIR/Pfam-A.hmm.h3m" ] || hmmpress "$PFAM_DIR/Pfam-A.hmm"
echo "profiles: $(grep -c '^NAME' "$PFAM_DIR/Pfam-A.hmm")"

echo; echo "== Arabidopsis proteome =="
if [ ! -s "$ATH_AA" ]; then
  curl -sSL -o "$ATH_AA" "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3A${ATH_ID}%29"
  date -u +"downloaded %Y-%m-%d from rest.uniprot.org, query proteome:${ATH_ID}" > "${ATH_AA}.source.txt"
fi
cat "${ATH_AA}.source.txt"; echo "sequences: $(grep -c '>' "$ATH_AA")"

echo; echo "== beet proteome =="
# stop codon characters at the end of the sequences are removed for the search tools
sed 's/\*$//' "$BEET_AA" > "$OUT/beet_proteins.fa"
echo "sequences: $(grep -c '>' "$OUT/beet_proteins.fa"); first header: $(head -n 1 "$OUT/beet_proteins.fa")"

echo; echo "== hmmsearch Pfam-A, gathering thresholds =="
hmmsearch --cpu "$T" --cut_ga --noali --domtblout "$OUT/pfam.domtblout" "$PFAM_DIR/Pfam-A.hmm" \
  "$OUT/beet_proteins.fa" > /dev/null
echo "domain hits: $(grep -vc '^#' "$OUT/pfam.domtblout"); proteins with a hit: $(grep -v '^#' "$OUT/pfam.domtblout" | awk '{print $1}' | sort -u | wc -l)"

echo; echo "== DIAMOND, both directions =="
diamond makedb --in "$ATH_AA" -d "$OUT/ath" --quiet
diamond makedb --in "$OUT/beet_proteins.fa" -d "$OUT/beet" --quiet
FMT="qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore"
diamond blastp -q "$OUT/beet_proteins.fa" -d "$OUT/ath" -o "$OUT/beet_vs_ath.tsv" --very-sensitive \
  --evalue 1e-5 --max-target-seqs 10 --threads "$T" --outfmt 6 $FMT --quiet
diamond blastp -q "$ATH_AA" -d "$OUT/beet" -o "$OUT/ath_vs_beet.tsv" --very-sensitive \
  --evalue 1e-5 --max-target-seqs 10 --threads "$T" --outfmt 6 $FMT --quiet
echo "$FMT" > "$OUT/diamond_columns.txt"
echo "beet queries with a hit: $(cut -f1 "$OUT/beet_vs_ath.tsv" | sort -u | wc -l); ath queries with a hit: $(cut -f1 "$OUT/ath_vs_beet.tsv" | sort -u | wc -l)"
echo "outputs: $OUT/pfam.domtblout  $OUT/beet_vs_ath.tsv  $OUT/ath_vs_beet.tsv"
