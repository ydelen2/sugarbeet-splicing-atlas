"""Long-read observation of both forms of splicing events, compared at equal long-read depth.

The chance of seeing both forms of an event in long reads grows with the expression of the gene, and the
high-confidence events sit in better covered genes because rMATS requires 20 reads in every replicate. The
comparison is therefore stratified by long-read depth: for every SE, A3, A5 and MX event of the SUPPA2 universe,
depth is the number of long reads on the junctions of either form (from longread_junctions.tsv). Within
depth strata the share of events with every junction of both forms observed is compared between the
high-confidence events and the other events of the universe, and the strata are combined with the
Mantel-Haenszel common odds ratio and test. Events with a depth of one read cannot show both forms and form
a stratum without information, which is left out of the test.

The high-confidence events also have, by definition, an intermediate inclusion level in at least one comparison,
and both forms of an event near PSI 0 or 1 are rarely seen whatever its coverage. With --psi (the SUPPA2
psiPerEvent table of the universe) the comparison is restricted to events whose mean PSI over the runs with a
defined value lies within --psi-min and --psi-max, so that the other events are matched on inclusion level as
well as on depth.

Usage
  python longread_event_depth.py --junctions junction_evidence/longread_junctions.tsv --ioe suppa/events/all_events.ioe \
         --events summary/events_high_confidence_annotated.tsv --out junction_evidence/longread_event_depth.tsv
  python longread_event_depth.py ... --psi suppa/psi/all_events.psi --psi-min 0.1 --psi-max 0.9 \
         --out junction_evidence/longread_event_depth_psi.tsv
"""
import argparse
import csv
import gzip
import os
import sys

import numpy as np
from scipy.stats import chi2
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from event_forms import parse_event  # noqa: E402

BINS = [(1, 1), (2, 4), (5, 9), (10, 19), (20, 49), (50, 99), (100, None)]


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--ioe", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--psi", help="SUPPA2 psiPerEvent table; restricts the events to a mean PSI range")
    ap.add_argument("--psi-min", type=float, default=0.1)
    ap.add_argument("--psi-max", type=float, default=0.9)
    a = ap.parse_args()

    keep = None
    if a.psi:
        keep = set()
        with opener(a.psi) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                vals = [float(x) for x in f[1:] if x not in ("nan", "NaN", "", "NA")]
                vals = [v for v in vals if v == v]
                if vals and a.psi_min <= sum(vals) / len(vals) <= a.psi_max:
                    keep.add(f[0])
        print(f"events with mean PSI in [{a.psi_min}, {a.psi_max}]: {len(keep)}")

    with opener(a.junctions) as fh:
        sup = {(r["chrom"], int(r["intron_start"]), int(r["intron_end"])): int(r["reads"])
               for r in csv.DictReader(fh, delimiter="\t")}
    with opener(a.events) as fh:
        hc = {r["suppa_event"] for r in csv.DictReader(fh, delimiter="\t")}
    rows = []
    with opener(a.ioe) as fh:
        fh.readline()
        for line in fh:
            ev = line.split("\t")[2]
            if ev.split(";", 1)[1].split(":", 1)[0] not in ("SE", "A3", "A5", "MX"):
                continue
            if keep is not None and ev not in keep:
                continue
            chrom, strand, etype, f1, f2, regions = parse_event(ev)
            k1 = [(chrom, x + 1, y - 1) for x, y in f1]
            k2 = [(chrom, x + 1, y - 1) for x, y in f2]
            depth = sum(sup.get(k, 0) for k in set(k1 + k2))
            if depth:
                rows.append((ev in hc, depth, all(k in sup for k in k1) and all(k in sup for k in k2)))

    out, tables = [], []
    for lo, hi in BINS:
        sel = [r for r in rows if r[1] >= lo and (hi is None or r[1] <= hi)]
        t = [sum(1 for r in sel if r[0] and r[2]), sum(1 for r in sel if r[0] and not r[2]),
             sum(1 for r in sel if not r[0] and r[2]), sum(1 for r in sel if not r[0] and not r[2])]
        informative = t[0] + t[2] > 0 and t[1] + t[3] > 0 and t[0] + t[1] > 0 and t[2] + t[3] > 0
        if informative:
            tables.append(np.array([[t[0], t[1]], [t[2], t[3]]]))
        out.append(dict(depth=f"{lo}-{hi if hi else 'max'}", hc_both=t[0], hc_events=t[0] + t[1],
                        hc_share=round(t[0] / (t[0] + t[1]), 4) if t[0] + t[1] else "",
                        other_both=t[2], other_events=t[2] + t[3],
                        other_share=round(t[2] / (t[2] + t[3]), 4) if t[2] + t[3] else "",
                        in_test="yes" if informative else "no (no event with both forms)"))
    st = StratifiedTable(tables)
    test = st.test_null_odds(correction=True)
    lo_ci, hi_ci = st.oddsratio_pooled_confint()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        psi_note = f"mean PSI {a.psi_min}-{a.psi_max}; " if a.psi else ""
        fh.write(f"# {psi_note}detected SE, A3, A5 and MX events {len(rows)}, high-confidence {sum(1 for r in rows if r[0])}; "
                 f"Mantel-Haenszel common odds ratio {st.oddsratio_pooled:.2f} (95% CI {lo_ci:.2f}-{hi_ci:.2f}), "
                 f"chi-square {test.statistic:.1f}, p {chi2.sf(test.statistic, 1):.2g}; "
                 f"homogeneity of odds ratios p {st.test_equal_odds().pvalue:.2g}\n")
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(open(a.out, encoding="utf-8").read())


if __name__ == "__main__":
    main()
