#!/bin/bash
# Summarize fastp and STAR results into one table and a MultiQC report. Run on a login node after 03_align.
#
# Junction annotation is taken from the FIRST-pass SJ.out.tab (saved as <run>.SJ.pass1.out.tab). In STAR
# two-pass mode the second pass treats junctions found in the first pass as annotated, so the
# "Annotated (sjdb)" figure in Log.final.out and column 6 of the final SJ.out.tab are always close to 100%
# (STAR issue #863). In the first pass, column 6 = 1 means the junction is in the GTF.
set -euo pipefail
source "$(dirname "$0")/config.sh"
activate_env
cd "$QC_DIR"

printf "run\tinput_reads\tuniq_pct\tmulti_pct\tmismatch_rate\tjunction_reads_uniq\tannotated_junction_read_pct\tnovel_junctions_3reads\n" > star_summary.tsv
for d in "${RUNRES_DIR}"/*/; do
  r=$(basename "$d"); L="$d/${r}.Log.final.out"; P1="$d/${r}.SJ.pass1.out.tab"
  [ -f "$L" ] || continue
  star=$(awk -F'|' '
    /Number of input reads/{ir=$2} /Uniquely mapped reads %/{u=$2} /% of reads mapped to multiple loci/{m=$2}
    /Mismatch rate per base/{mm=$2}
    END{gsub(/[ \t%]/,"",u); gsub(/[ \t%]/,"",m); gsub(/[ \t]/,"",ir); gsub(/[ \t%]/,"",mm);
        printf "%s\t%s\t%s\t%s", ir, u, m, mm}' "$L")
  if [ -f "$P1" ]; then
    sj=$(awk '{t+=$7; if ($6==1) a+=$7; else if ($7>=3) n++}
              END{printf "%d\t%.1f\t%d", t, (t>0 ? 100*a/t : 0), n}' "$P1")
  else
    sj="NA\tNA\tNA"
  fi
  printf "%s\t%s\t%b\n" "$r" "$star" "$sj" >> star_summary.tsv
done
awk 'NR<=50' star_summary.tsv | column -t
echo "runs summarized: $(awk 'END{print NR-1}' star_summary.tsv)"
multiqc -f -o "$QC_DIR/multiqc" "$QC_DIR/fastp" "$RUNRES_DIR" >/dev/null 2>&1 && echo "MultiQC: $QC_DIR/multiqc/multiqc_report.html"
