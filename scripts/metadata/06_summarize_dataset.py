"""Summarize the sample sheets: runs, BioProjects, compressed FASTQ size and sequenced bases.

Usage:
  python scripts/metadata/06_summarize_dataset.py
"""
import collections
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
META = os.path.join(REPO, "metadata")


def load(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    raw = {r["run_accession"]: r for r in load(os.path.join(META, "raw", "ena_runs_raw.tsv"))}
    total_bytes = sum(sum(float(x) for x in r["fastq_bytes"].split(";") if x) for r in raw.values())
    print(f"all transcriptomic runs: {len(raw)}  compressed FASTQ: {total_bytes / 1e9:.1f} GB")

    for sheet in ("samples_pilot.tsv", "samples_all.tsv"):
        rows = load(os.path.join(META, sheet))
        fq = sum(sum(float(x) for x in raw[r["run"]]["fastq_bytes"].split(";") if x) for r in rows)
        bases = sum(float(raw[r["run"]]["base_count"]) for r in rows)
        print(f"\n{sheet}: runs={len(rows)}  bioprojects={len(set(r['study'] for r in rows))}  "
              f"compressed FASTQ={fq / 1e9:.1f} GB  bases={bases / 1e9:.1f} Gb")
        by = collections.Counter((r["study"], r["condition"]) for r in rows)
        for (study, cond), n in sorted(by.items()):
            print(f"  {study}\t{cond}\t{n}")


if __name__ == "__main__":
    main()
