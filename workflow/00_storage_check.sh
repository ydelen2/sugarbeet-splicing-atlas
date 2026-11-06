#!/bin/bash
# Run on a Swan login node BEFORE anything else:  bash workflow/00_storage_check.sh
# Shows quotas/usage for $WORK and $NRDSTOR and a quick sequential write/read speed test on both.
set -uo pipefail

echo "== variables =="
echo "WORK    = ${WORK:-<unset>}"
echo "NRDSTOR = ${NRDSTOR:-<unset>}"
[ -n "${NRDSTOR:-}" ] && [ -d "$NRDSTOR" ] || echo "!! NRDStor not mounted here. Request access: https://hcc.unl.edu/docs/handling_data/data_storage/nrdstor/"

echo; echo "== group usage / quota (hcc-du) =="
command -v hcc-du >/dev/null && hcc-du || echo "hcc-du not found"

echo; echo "== df =="
df -h "${WORK:-/work}" "${NRDSTOR:-/mnt/nrdstor}" 2>/dev/null

speed() { # dir label
  local d="$1/.bvas_speedtest.$$"
  mkdir -p "$d" || return
  echo; echo "== $2 : write 4 GB =="
  sync; dd if=/dev/zero of="$d/test.bin" bs=1M count=4096 conv=fdatasync 2>&1 | tail -1
  echo "== $2 : read 4 GB =="
  dd if="$d/test.bin" of=/dev/null bs=1M 2>&1 | tail -1
  echo "== $2 : create 2000 small files =="
  /usr/bin/time -f "%e s" bash -c "for i in \$(seq 2000); do echo x > $d/f\$i; done" 2>&1 | tail -1
  rm -rf "$d"
}
[ -n "${WORK:-}" ] && speed "$WORK" "WORK"
[ -n "${NRDSTOR:-}" ] && [ -d "$NRDSTOR" ] && speed "$NRDSTOR" "NRDSTOR"

echo; echo "BULK_ON_NRDSTOR=1 in workflow/config.sh keeps FASTQ and BAM files on NRDSTOR. It needs about 6 TB of free quota"
echo "and a write speed comparable to WORK; otherwise set BULK_ON_NRDSTOR=0."
