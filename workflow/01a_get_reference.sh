#!/bin/bash
# Download RefBeet-3.0 and BeetSet-3, index the genome and convert the gene models to GTF.
# Run on a login node from $NRDSTOR/bvas:  bash workflow/01a_get_reference.sh
# Safe to rerun: files already downloaded are skipped.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"

LOGFILE="$LOG_DIR/get_reference.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
echo "log: $LOGFILE"
trap 'echo "FAILED at line $LINENO (exit $?)"' ERR

activate_env
cd "$REF_DIR"

BASE="https://bvseq.boku.ac.at/Genome/Download/RefBeet-3.0"
for f in README.txt RefBeet-3.0.fa.gz BeetSet-3.genes.2604.gff BeetSet-3.genes.2604.cds \
         BeetSet-3.genes.2604.mrna BeetSet-3.genes.2604.aa BeetSet-3.genes.2604.functional \
         RefBeet-3.0.ncRNA.annotation; do
  [ -s "$f" ] || wget -c -q --show-progress "${BASE}/${f}"
done
md5sum README.txt RefBeet-3.0.fa.gz BeetSet-3.genes.2604.* RefBeet-3.0.ncRNA.annotation > downloaded.md5

[ -s RefBeet-3.0.fa ] || pigz -dk RefBeet-3.0.fa.gz
samtools faidx RefBeet-3.0.fa

# Note: awk reads its whole input, so these summaries never trigger SIGPIPE under "set -o pipefail"
echo "== chromosomes =="
awk -F'\t' '$1 ~ /^0_Bvchr/ {printf "%s\t%d\n", $1, $2}' RefBeet-3.0.fa.fai
awk -F'\t' '{n++; s+=$2} $1 ~ /^0_Bvchr/ {c++; cs+=$2}
  END{printf "sequences: %d  total: %d bp  chromosome-level: %d (%d bp, %.1f%%)  other: %d (%d bp)\n",
      n, s, c, cs, 100*cs/s, n-c, s-cs}' RefBeet-3.0.fa.fai

echo "== GFF feature types =="
awk -F'\t' '!/^#/ && NF>=9 {n[$3]++} END{for(k in n) printf "%9d %s\n", n[k], k}' "$GENES_GFF" | sort -k2
echo "== GFF attribute examples (first gene and mRNA) =="
awk -F'\t' '!/^#/ && $3=="gene" && !g {print $9; g=1} !/^#/ && $3=="mRNA" && !m {print $9; m=1} g && m {exit}' "$GENES_GFF" || true

# GFF3 -> GTF (STAR, StringTie, rMATS and SUPPA2 all read GTF)
gffread "$GENES_GFF" -T -o "$GENES_GTF"

echo "== GTF summary =="
awk -F'\t' '$3=="exon" {
    t=""; g="";
    if (match($9, /transcript_id "[^"]+"/)) t=substr($9, RSTART+15, RLENGTH-16);
    if (match($9, /gene_id "[^"]+"/))       g=substr($9, RSTART+9,  RLENGTH-10);
    ex++; tx[t]=1; if (g=="") nog[t]=1; else tg[t]=g
  }
  END{
    nt=0; for (t in tx) nt++;
    nm=0; for (t in nog) nm++;
    for (t in tg) gt[tg[t]]++;
    ng=0; for (g in gt) { ng++; dist[gt[g]]++ }
    printf "exons: %d\ntranscripts: %d\ngenes (gene_id): %d\ntranscripts without gene_id: %d\n", ex, nt, ng, nm;
    printf "transcripts per gene: %.3f\n", (ng>0 ? (nt-nm)/ng : 0);
    print "distribution (transcripts per gene: number of genes)";
    for (k in dist) printf "  %d: %d\n", k, dist[k]
  }' "$GENES_GTF"
echo "published BeetSet-3 values (Dohm et al. 2025): 28271 genes, 32609 transcripts"

GLEN=$(awk '{s+=$2} END{print s}' RefBeet-3.0.fa.fai)
echo "genome length: $GLEN"
echo "Next: sbatch workflow/01b_star_index.sbatch"
