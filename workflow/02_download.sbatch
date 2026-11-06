#!/bin/bash
#SBATCH --job-name=ena_dl
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=4G
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
# Usage (from $NRDSTOR/bvas):
#   N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%8 workflow/02_download.sbatch
# Downloads one run (2 FASTQ) per task from ENA over HTTPS into bulk storage and verifies MD5.
# If compute nodes have no outbound internet, use 02_download_xfer.sh on swan-xfer instead.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
activate_env

i=${SLURM_ARRAY_TASK_ID}
skip_if_excluded "$i"
RUN=$(sample_field "$i" run)
FQ1=$(sample_field "$i" fq1); FQ2=$(sample_field "$i" fq2)
M1=$(sample_field "$i" md5_1); M2=$(sample_field "$i" md5_2)
OUT="${FASTQ_DIR}/${RUN}"; mkdir -p "$OUT"; cd "$OUT"

# ENA refuses connections when too many are open at once. Use 2 connections per file, wait between
# attempts with increasing delays, and keep partial files so an interrupted download resumes.
fetch() { # url md5 dest
  local url="$1" md5="$2" dest="$3" attempt
  if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then echo "ok (cached) $dest"; return 0; fi
  # A control file without its data file is left over from an earlier failed job: start clean.
  [ -e "$dest" ] || rm -f "$dest.aria2"
  for attempt in 1 2 3 4 5 6; do
    # aria2c only, no wget fallback: aria2c can leave a full-size file with missing pieces after a
    # broken connection, and "wget -c" would then report it as fully retrieved.
    aria2c -x 2 -s 2 --continue=true --file-allocation=none --max-tries=3 --retry-wait=60 \
           --connect-timeout=60 --console-log-level=warn --summary-interval=0 -o "$dest" "$url" || true
    if [ -s "$dest" ] && [ ! -e "$dest.aria2" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then
      echo "ok $dest (attempt $attempt)"; return 0
    fi
    if [ -s "$dest" ] && [ ! -e "$dest.aria2" ]; then
      echo "attempt $attempt: complete file with wrong md5, deleting $dest"; rm -f "$dest"
    else
      echo "attempt $attempt: download incomplete for $dest"
    fi
    sleep $((attempt * 120))
  done
  echo "FAILED $url" >&2; return 1
}
# spread the start of the tasks so they do not all connect in the same second
sleep $(( (SLURM_ARRAY_TASK_ID * 7) % 60 ))
fetch "$FQ1" "$M1" "${RUN}_1.fastq.gz"
fetch "$FQ2" "$M2" "${RUN}_2.fastq.gz"
touch "${OUT}/.download_ok"
