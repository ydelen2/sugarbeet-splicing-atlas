"""Select ENA runs suitable for splice junction analysis.

Criteria (all must hold):
  library_strategy is RNA-Seq or ssRNA-seq
  library_layout is PAIRED
  mean read length per mate >= 90 bp (base_count / read_count / 2)
  sequenced bases per run >= 2 Gb

Usage:
  python scripts/metadata/02_filter_runs.py [--raw metadata/raw]
Output:
  <raw>/samples_candidate.tsv
"""
import argparse
import collections
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MIN_READ_LEN = 90
MIN_GB = 2.0
OUT_COLS = ["run_accession", "study_accession", "sample_accession", "library_layout", "library_strategy",
            "instrument_model", "read_len", "gb", "first_public", "study_title", "sample_title", "sample_alias",
            "experiment_title", "sample_description", "tissue_type", "dev_stage", "cultivar", "fastq_ftp",
            "fastq_md5"]


def read_len(r):
    rc = float(r["read_count"] or 0)
    bc = float(r["base_count"] or 0)
    if rc == 0:
        return 0
    mean = bc / rc
    return round(mean / 2) if r["library_layout"] == "PAIRED" else round(mean)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(REPO, "metadata", "raw"))
    args = ap.parse_args()

    with open(os.path.join(args.raw, "ena_runs_raw.tsv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    keep, dropped = [], collections.Counter()
    for r in rows:
        r["read_len"] = read_len(r)
        r["gb"] = round(float(r["base_count"] or 0) / 1e9, 2)
        if r["library_strategy"] not in ("RNA-Seq", "ssRNA-seq"):
            dropped["library_strategy"] += 1
        elif r["library_layout"] != "PAIRED":
            dropped["single_end"] += 1
        elif r["read_len"] < MIN_READ_LEN:
            dropped[f"read_len<{MIN_READ_LEN}"] += 1
        elif r["gb"] < MIN_GB:
            dropped[f"bases<{MIN_GB}Gb"] += 1
        else:
            keep.append(r)

    keep.sort(key=lambda r: (r["study_accession"], r["run_accession"]))
    out = os.path.join(args.raw, "samples_candidate.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(keep)

    print(f"input runs: {len(rows)}  kept: {len(keep)}  dropped: {dict(dropped)}")
    by = collections.Counter(r["study_accession"] for r in keep)
    for study, n in by.most_common():
        print(f"  {study}\t{n}")


if __name__ == "__main__":
    main()
