#!/bin/bash
# Alternative downloader for the transfer node (swan-xfer.unl.edu) or login node, 2 parallel runs.
# Usage:  nohup bash 02_download_xfer.sh metadata/samples_pilot.tsv > logs/download_xfer.log 2>&1 &
set -euo pipefail
source "$(dirname "$0")/config.sh"
SHEET="${1:-$SAMPLES}"
command -v aria2c >/dev/null 2>&1 || { module load anaconda; conda activate "$ENV_DIR"; }

dl_one() {
  local run="$1" fq1="$2" fq2="$3" m1="$4" m2="$5"
  local out="${FASTQ_DIR}/${run}"; mkdir -p "$out"
  for pair in "$fq1|$m1|${run}_1.fastq.gz" "$fq2|$m2|${run}_2.fastq.gz"; do
    IFS='|' read -r url md5 dest <<< "$pair"
    if [ -s "$out/$dest" ] && [ "$(md5sum "$out/$dest" | cut -d' ' -f1)" = "$md5" ]; then continue; fi
    # Same logic as 02_download.sbatch: aria2c only, 2 connections, waits between attempts, and a file
    # with a leftover .aria2 control file is never accepted as complete.
    [ -e "$out/$dest" ] || rm -f "$out/$dest.aria2"
    local ok=0 a
    for a in 1 2 3 4 5 6; do
      aria2c -q -x 2 -s 2 --continue=true --file-allocation=none --max-tries=3 --retry-wait=60 \
             --connect-timeout=60 -d "$out" -o "$dest" "$url" || true
      if [ -s "$out/$dest" ] && [ ! -e "$out/$dest.aria2" ]; then
        if [ "$(md5sum "$out/$dest" | cut -d' ' -f1)" = "$md5" ]; then ok=1; break; fi
        rm -f "$out/$dest"
      fi
      sleep $((a * 120))
    done
    [ "$ok" = 1 ] || { echo "FAILED $run $dest"; return 1; }
  done
  touch "$out/.download_ok"; echo "done $run"
}
export -f dl_one; export FASTQ_DIR

# Columns are looked up by name, so adding columns to the sample sheet does not break this script
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
            {print $c["run"] "\t" $c["fq1"] "\t" $c["fq2"] "\t" $c["md5_1"] "\t" $c["md5_2"]}' "$SHEET" | \
  xargs -P 2 -L 1 bash -c 'dl_one "$0" "$1" "$2" "$3" "$4"'
