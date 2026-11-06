#!/bin/bash
# Resubmit the download for every run of the sample sheet that has no .download_ok marker yet.
# Run on a login node after the main download array has finished:
#   ANALYSIS=full bash workflow/02b_retry_downloads.sh [max_parallel_tasks]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"
PAR="${1:-4}"

if squeue -u "$USER" -h -n ena_dl | grep -q .; then
  echo "download jobs are still running or pending; wait for them to finish first" >&2
  exit 1
fi

# runs marked use_das=excluded are never counted as missing
IDX=$(awk -F'\t' -v d="$FASTQ_DIR" \
      'NR==1 {for (i=1; i<=NF; i++) {if ($i=="run") c=i; if ($i=="use_das") u=i}; next}
       u && $u=="excluded" {next}
       { f = d "/" $c "/.download_ok"; if ((getline line < f) < 0) printf "%s%d", (n++ ? "," : ""), NR-1; close(f) }' "$SAMPLES")
TOTAL=$(awk 'END {print NR-1}' "$SAMPLES")
if [ -z "$IDX" ]; then echo "all $TOTAL runs are downloaded"; exit 0; fi
N=$(echo "$IDX" | tr ',' '\n' | awk 'END {print NR}')
echo "runs still missing: $N of $TOTAL"
cd "$PROJ"
JOB=$(ANALYSIS="$ANALYSIS" sbatch --parsable --array="${IDX}%${PAR}" workflow/02_download.sbatch)
echo "resubmitted as job $JOB with at most $PAR tasks in parallel"
