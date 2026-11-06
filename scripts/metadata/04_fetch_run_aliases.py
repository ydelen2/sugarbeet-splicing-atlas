"""Fetch submitter aliases (run, experiment, library, sample) for candidate runs from ENA.

Several BioProjects encode the experimental condition only in these aliases
(for example "4_C_10w_Tap_GT1_Rep1" or "G1_30hpi_Hokkaido_C").

Usage:
  python scripts/metadata/04_fetch_run_aliases.py [--raw metadata/raw]
Output:
  <raw>/run_aliases.tsv
"""
import argparse
import csv
import os
import urllib.parse
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIELDS = ["run_accession", "study_accession", "run_alias", "experiment_alias", "experiment_title",
          "library_name", "sample_alias", "sample_title", "submitted_ftp"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(REPO, "metadata", "raw"))
    args = ap.parse_args()

    with open(os.path.join(args.raw, "samples_candidate.tsv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    studies = sorted(set(r["study_accession"] for r in rows))

    alias = {}
    for st in studies:
        params = urllib.parse.urlencode({
            "result": "read_run", "query": f'study_accession="{st}"', "fields": ",".join(FIELDS),
            "format": "tsv", "limit": 0,
        })
        txt = urllib.request.urlopen("https://www.ebi.ac.uk/ena/portal/api/search?" + params, timeout=120)
        for r in csv.DictReader(txt.read().decode("utf-8").splitlines(), delimiter="\t"):
            alias[r["run_accession"]] = r

    out = os.path.join(args.raw, "run_aliases.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            if r["run_accession"] in alias:
                w.writerow(alias[r["run_accession"]])
    print(f"{len(studies)} studies queried; aliases written to {out}")


if __name__ == "__main__":
    main()
