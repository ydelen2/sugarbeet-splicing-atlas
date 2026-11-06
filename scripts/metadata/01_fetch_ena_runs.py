"""Download the run-level metadata of all public Beta vulgaris transcriptomic runs from ENA.

Query: tax_tree(3555) AND library_source="TRANSCRIPTOMIC" (ENA Portal API, result=read_run).
ENA is a live archive, so the retrieval date is written next to the table.

Usage:
  python scripts/metadata/01_fetch_ena_runs.py [--outdir metadata/raw]
Output:
  <outdir>/ena_runs_raw.tsv
  <outdir>/ena_runs_raw.retrieved.txt
"""
import argparse
import datetime
import os
import urllib.parse
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FIELDS = [
    "run_accession", "study_accession", "secondary_study_accession", "sample_accession",
    "experiment_accession", "tax_id", "scientific_name", "instrument_platform", "instrument_model",
    "library_layout", "library_strategy", "library_source", "library_selection", "read_count",
    "base_count", "first_public", "study_title", "sample_title", "sample_alias", "experiment_title",
    "sample_description", "tissue_type", "dev_stage", "cultivar", "isolate", "fastq_ftp",
    "fastq_bytes", "fastq_md5",
]
QUERY = 'tax_tree(3555) AND library_source="TRANSCRIPTOMIC"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(REPO, "metadata", "raw"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    params = urllib.parse.urlencode({
        "result": "read_run", "query": QUERY, "fields": ",".join(FIELDS),
        "format": "tsv", "limit": 0,
    })
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + params
    data = urllib.request.urlopen(url, timeout=300).read().decode("utf-8")
    data = data.replace("\r\n", "\n")

    out = os.path.join(args.outdir, "ena_runs_raw.tsv")
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(data)
    n = data.count("\n") - 1
    with open(os.path.join(args.outdir, "ena_runs_raw.retrieved.txt"), "w", encoding="utf-8", newline="\n") as fh:
        now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
        fh.write(f"retrieved_utc\t{now}\n")
        fh.write(f"query\t{QUERY}\nruns\t{n}\nurl\t{url}\n")
    print(f"{n} runs written to {out}")


if __name__ == "__main__":
    main()
