#!/bin/bash
# Report, and optionally resubmit, the runs of the sample sheet that have not finished a step.
#   ANALYSIS=full bash workflow/retry_missing.sh STEP [max_parallel] [--submit]
# STEP is one of: download, align, stringtie, salmon. Without --submit the script only reports.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"
STEP="${1:?give a step: download, align, stringtie or salmon}"
PAR="${2:-10}"
SUBMIT="${3:-}"

case "$STEP" in
  download)  DIR="$FASTQ_DIR";  MARK=".download_ok";  SCRIPT="02_download.sbatch";     JOBNAME="ena_dl" ;;
  align)     DIR="$RUNRES_DIR"; MARK=".align_ok";     SCRIPT="03_align.sbatch";        JOBNAME="bv_align" ;;
  stringtie) DIR="$ST_DIR";     MARK=".stringtie_ok"; SCRIPT="04a_stringtie.sbatch";   JOBNAME="bv_stringtie" ;;
  salmon)    DIR="$SALMON_DIR"; MARK=".salmon_ok";    SCRIPT="06_salmon_quant.sbatch"; JOBNAME="bv_salmon" ;;
  *) echo "unknown step '$STEP'" >&2; exit 1 ;;
esac

TOTAL=$(awk 'END {print NR-1}' "$SAMPLES")
# runs marked use_das=excluded are never counted as missing
IDX=$(awk -F'\t' -v d="$DIR" -v m="$MARK" \
      'NR==1 {for (i=1; i<=NF; i++) {if ($i=="run") c=i; if ($i=="use_das") u=i}; next}
       u && $u=="excluded" {next}
       { f = d "/" $c "/" m; if ((getline line < f) < 0) printf "%s%d", (n++ ? "," : ""), NR-1; close(f) }' "$SAMPLES")
if [ -z "$IDX" ]; then echo "[$ANALYSIS] $STEP: all $TOTAL runs done"; exit 0; fi
N=$(echo "$IDX" | tr ',' '\n' | awk 'END {print NR}')
echo "[$ANALYSIS] $STEP: $N of $TOTAL runs not done"
echo "$IDX" | tr ',' '\n' | awk -v s="$SAMPLES" 'BEGIN {while ((getline l < s) > 0) {k++; split(l, a, "\t"); run[k-1]=a[1] "\t" a[2]}}
     NR<=15 {print "  task " $1 "\t" run[$1]} END {if (NR>15) print "  ... and " NR-15 " more"}'

[ "$SUBMIT" = "--submit" ] || { echo "add --submit to resubmit these tasks"; exit 0; }
if squeue -u "$USER" -h -n "$JOBNAME" | grep -q .; then
  echo "$JOBNAME jobs are still running or pending; wait for them to finish first" >&2
  exit 1
fi
cd "$PROJ"
JOB=$(ANALYSIS="$ANALYSIS" sbatch --parsable --array="${IDX}%${PAR}" "workflow/$SCRIPT")
echo "resubmitted as job $JOB with at most $PAR tasks in parallel"
