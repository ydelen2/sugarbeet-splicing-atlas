"""Salmon quality summary for every run of a sample sheet.

For each run: Salmon mapping rate, the library type Salmon detected, and the strandedness inferred
earlier from STAR gene counts (04a_stringtie writes it to per_run/<run>/strandedness.tsv). The two
library calls are compared (STAR unstranded = IU, reverse = ISR, forward = ISF for paired-end data).

Usage
  python salmon_qc.py --sheet samples_all.tsv --salmon-dir results_full/salmon \
         --runres-dir results_full/per_run --out results_full/summary/salmon_qc.tsv
"""
import argparse
import collections
import csv
import json
import os

STAR_TO_SALMON = {"unstranded": "IU", "reverse": "ISR", "forward": "ISF"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--salmon-dir", required=True)
    ap.add_argument("--runres-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.sheet, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("use_das") != "excluded"]

    out = []
    for r in rows:
        run = r["run"]
        d = os.path.join(a.salmon_dir, run)
        rec = dict(run=run, study=r["study"], condition=r["condition"], use_das=r["use_das"],
                   percent_mapped="", salmon_libtype="", libtype_setting="", star_strand="", agree="")
        try:
            with open(os.path.join(d, "aux_info", "meta_info.json"), encoding="utf-8") as fh:
                rec["percent_mapped"] = round(float(json.load(fh)["percent_mapped"]), 2)
            with open(os.path.join(d, "lib_format_counts.json"), encoding="utf-8") as fh:
                rec["salmon_libtype"] = json.load(fh)["expected_format"]
        except (OSError, KeyError, ValueError):
            rec["salmon_libtype"] = "missing"
        sf = os.path.join(a.runres_dir, run, "strandedness.tsv")
        if os.path.isfile(sf):
            with open(sf, encoding="utf-8") as fh:
                rec["star_strand"] = fh.read().split("\t")[1]
        exp = STAR_TO_SALMON.get(rec["star_strand"], "")
        rec["agree"] = "yes" if exp and exp == rec["salmon_libtype"] else "no"
        # with a fixed library type (-l ISR) Salmon reports that type as "expected_format" as well
        try:
            with open(os.path.join(d, "cmd_info.json"), encoding="utf-8") as fh:
                lt = json.load(fh).get("libType", "")
            rec["libtype_setting"] = "unknown" if not lt else ("automatic" if lt == "A" else f"fixed {lt}")
        except (OSError, ValueError):
            rec["libtype_setting"] = ""
        out.append(rec)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    cols = ["run", "study", "condition", "use_das", "percent_mapped", "salmon_libtype", "libtype_setting",
            "star_strand", "agree"]
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    mapped = [x["percent_mapped"] for x in out if x["percent_mapped"] != ""]
    print(f"runs: {len(out)}  with Salmon output: {len(mapped)}")
    if mapped:
        s = sorted(mapped)
        print(f"percent mapped: median {s[len(s) // 2]:.1f}, min {s[0]:.1f}, max {s[-1]:.1f}, "
              f"below 50: {sum(v < 50 for v in s)}")
    print(f"library call agrees with STAR: {sum(x['agree'] == 'yes' for x in out)} of {len(out)}")
    for x in out:
        if x["agree"] != "yes":
            print(f"  disagree: {x['run']}\t{x['study']}\tSalmon {x['salmon_libtype']}\tSTAR {x['star_strand']}")
    print("\nlibrary type per study (Salmon call: runs):")
    by = collections.defaultdict(collections.Counter)
    for x in out:
        by[x["study"]][x["salmon_libtype"]] += 1
    for st in sorted(by):
        print(f"  {st}\t" + "  ".join(f"{k}: {v}" for k, v in sorted(by[st].items())))


if __name__ == "__main__":
    main()
