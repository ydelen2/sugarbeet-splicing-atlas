#!/bin/bash
#SBATCH --job-name=bv_merge
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# Usage (from $NRDSTOR/bvas):  sbatch workflow/04b_merge.sbatch
# Merge per-run assemblies, compare to BeetSet-3 and pool junctions. All outputs on NRDStor.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env
MERGE="${RES_DIR}/merged"; mkdir -p "$MERGE"
TMP=$(job_tmp 20 merge); trap 'rm -rf "$TMP"' EXIT
cd "$TMP"

# Runs used for transcript merging and junction pooling: those with a STAR unique mapping rate of at
# least MERGE_MIN_UNIQ percent (config.sh). Low-rate runs stay aligned and keep their role in the
# comparisons; they only do not contribute to the reference transcriptome. The selection is written
# to merge_runs.tsv (run, unique_pct, used).
[ -s "${QC_DIR}/star_summary.tsv" ] || { echo "run 05_qc_summary.sh first: ${QC_DIR}/star_summary.tsv missing" >&2; exit 1; }
awk -F'\t' -v min="$MERGE_MIN_UNIQ" 'NR==1 {print "run\tunique_pct\tused"; next}
     {print $1 "\t" $3 "\t" (($3+0 >= min) ? "yes" : "no")}' "${QC_DIR}/star_summary.tsv" > merge_runs.tsv
awk -F'\t' 'NR>1 && $3=="yes" {print $1}' merge_runs.tsv | sort > used_runs.txt
echo "runs with unique mapping >= ${MERGE_MIN_UNIQ}%: $(awk 'END{print NR}' used_runs.txt) of $(awk 'END{print NR-1}' merge_runs.tsv)"
echo "excluded from merging:"; awk -F'\t' 'NR>1 && $3=="no" {printf "  %s\t%s%%\n", $1, $2}' merge_runs.tsv

: > gtf_list.txt
while read -r r; do
  f="${ST_DIR}/${r}/${r}.gtf"
  [ -s "$f" ] && echo "$f" >> gtf_list.txt || echo "no assembly for $r (skipped)"
done < used_runs.txt
echo "assemblies merged: $(wc -l < gtf_list.txt)"
stringtie --merge -p 8 -G "$GENES_GTF" -o merged.gtf -m 200 -c 1 -F 0.5 -T 1 -f 0.05 gtf_list.txt
gffcompare -r "$GENES_GTF" -R -o gffcmp merged.gtf
echo "== gffcompare class codes (= known, j novel isoform of known gene, u intergenic, i intronic, x antisense) =="
{ printf "class\tmulti_exon\tsingle_exon\n"
  awk -F'\t' '
    $3=="transcript" { if (match($9, /class_code "[^"]+"/)) c=substr($9, RSTART+12, RLENGTH-13); else c="none";
                       match($9, /transcript_id "[^"]+"/); t=substr($9, RSTART+15, RLENGTH-16); code[t]=c }
    $3=="exon"       { match($9, /transcript_id "[^"]+"/); t=substr($9, RSTART+15, RLENGTH-16); ex[t]++ }
    END { for (t in code) { k=code[t]; if (ex[t]>1) m[k]++; else s[k]++; all[k]=1 }
          for (k in all) printf "%s\t%d\t%d\n", k, m[k], s[k] }' gffcmp.annotated.gtf | sort -t$'\t' -k2,2nr
} | tee class_codes.txt

# The Salmon index is built in 04c_filter_index on the structurally filtered transcriptome.

# Reference junctions from the BeetSet-3 GTF: the intron between consecutive exons of each transcript
awk -F'\t' '$3=="exon" {match($9, /transcript_id "[^"]+"/); t=substr($9, RSTART+15, RLENGTH-16); print t "\t" $1 "\t" $4 "\t" $5}' "$GENES_GTF" \
  | sort -k1,1 -k3,3n \
  | awk -F'\t' '$1==pt && $2==pc {print $2 "\t" pe+1 "\t" $3-1} {pt=$1; pc=$2; pe=$4}' | sort -u > ref_junctions.tsv
echo "reference junctions (BeetSet-3): $(awk 'END{print NR}' ref_junctions.tsv)"

# Pooled junctions from the final (second-pass) SJ.out.tab files of the runs used for merging: at least
# JUNCTION_MIN_READS uniquely mapped reads in a run, and seen in at least JUNCTION_MIN_RUN_FRACTION of
# those runs (never fewer than 2). Column 6 of SJ.out.tab is not used, because in two-pass mode it also
# marks first-pass junctions as annotated; annotation status is taken from the GTF instead.
N_USED=$(awk 'END{print NR}' used_runs.txt)
MIN_RUNS_J=$(awk -v n="$N_USED" -v f="$JUNCTION_MIN_RUN_FRACTION" 'BEGIN {m=int(n*f + 0.999999); print (m<2 ? 2 : m)}')
echo "junction pooling: >= ${JUNCTION_MIN_READS} reads in a run, present in >= ${MIN_RUNS_J} of ${N_USED} runs"
while read -r r; do
  f="${RUNRES_DIR}/${r}/${r}.SJ.out.tab"
  [ -s "$f" ] && awk -v m="$JUNCTION_MIN_READS" '$7>=m {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' "$f"
done < used_runs.txt | sort | uniq -c | awk -v m="$MIN_RUNS_J" '$1>=m {print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1}' > pooled_tmp.tsv
{ printf "chrom\tintron_start\tintron_end\tstrand_code\tmotif_code\tn_runs\tstatus\n"
  awk -F'\t' 'NR==FNR {ref[$1 "\t" $2 "\t" $3]=1; next}
              {print $0 "\t" ((($1 "\t" $2 "\t" $3) in ref) ? "annotated" : "novel")}' ref_junctions.tsv pooled_tmp.tsv
} > pooled_junctions.tsv
awk -F'\t' 'NR>1 {n[$7]++} END{for (k in n) print "pooled junctions, " k ": " n[k]}' pooled_junctions.tsv

awk -F'\t' 'NR>1 && $7=="novel" {n[($6>=50)?"50+":($6>=20)?"20-49":($6>=10)?"10-19":($6>=5)?"5-9":"<5"]++}
     END {for (k in n) print "novel junctions supported by " k " runs: " n[k]}' pooled_junctions.tsv | sort
cp merged.gtf gffcmp.* class_codes.txt ref_junctions.tsv pooled_junctions.tsv merge_runs.tsv gtf_list.txt "$MERGE/"
echo "merge done -> $MERGE"
