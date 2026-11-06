"""Fetch BioSample attributes (tissue, treatment, cultivar, replicate, ...) for candidate runs from ENA.

Usage:
  python scripts/metadata/03_fetch_biosample_attributes.py [--raw metadata/raw]
Output:
  <raw>/samples_candidate_attrs.tsv  (one row per run; attributes joined as "key=value | key=value")
"""
import argparse
import collections
import csv
import os
import time
import urllib.request
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BATCH = 20
RETRIES = 6


def fetch_samples(accessions):
    """Return the list of SAMPLE elements for the given accessions."""
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(accessions)
    for attempt in range(1, RETRIES + 1):
        try:
            data = urllib.request.urlopen(url, timeout=300).read()
            return ET.fromstring(data).findall("SAMPLE")
        except Exception as exc:  # the ENA browser API times out occasionally
            wait = 15 * attempt
            print(f"  attempt {attempt}/{RETRIES} failed for {accessions[0]} ({exc}); waiting {wait} s")
            time.sleep(wait)
    if len(accessions) > 1:
        print(f"  falling back to single requests for {len(accessions)} samples")
        out = []
        for acc in accessions:
            out.extend(fetch_samples([acc]))
        return out
    raise RuntimeError(f"could not fetch {accessions[0]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(REPO, "metadata", "raw"))
    args = ap.parse_args()

    with open(os.path.join(args.raw, "samples_candidate.tsv"), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = sorted(set(r["sample_accession"] for r in rows))
    print(f"unique samples: {len(samples)}")

    attrs = {}
    for i in range(0, len(samples), BATCH):
        for s in fetch_samples(samples[i:i + BATCH]):
            d = collections.OrderedDict()
            t = s.find("TITLE")
            d["TITLE"] = t.text if t is not None and t.text else ""
            for a in s.findall(".//SAMPLE_ATTRIBUTE"):
                tag, val = a.find("TAG"), a.find("VALUE")
                if tag is not None and val is not None and tag.text and not tag.text.startswith("ENA-"):
                    d[tag.text] = val.text or ""
            ids = [s.get("accession")] + [x.text for x in s.findall(".//EXTERNAL_ID")] + \
                  [x.text for x in s.findall(".//PRIMARY_ID")]
            for k in ids:
                if k:
                    attrs[k] = d
        print(f"  {min(i + BATCH, len(samples))}/{len(samples)}")

    out = os.path.join(args.raw, "samples_candidate_attrs.tsv")
    missing = 0
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "study_accession", "sample_accession", "read_len", "gb", "attributes"])
        for r in rows:
            d = attrs.get(r["sample_accession"])
            if d is None:
                missing += 1
                d = {}
            w.writerow([r["run_accession"], r["study_accession"], r["sample_accession"], r["read_len"], r["gb"],
                        " | ".join(f"{k}={v}" for k, v in d.items())])
    print(f"written {out}; runs without BioSample record: {missing}")


if __name__ == "__main__":
    main()
