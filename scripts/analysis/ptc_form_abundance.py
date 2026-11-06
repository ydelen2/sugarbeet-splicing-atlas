"""Abundance of the premature-stop form of the high-confidence events across all quantified runs.

If transcripts with a premature termination codon are degraded by nonsense-mediated decay, the form of an event
that carries the premature stop should be the less abundant form more often than a form without one. The
annotated form of an event is, by construction, the one that ends at the annotated stop codon (nmd_features.py), so
the comparison is made among the non-annotated forms only: for every coding event with both forms reconstructed,
the share of the non-annotated form in the event PSI (SUPPA2 psiPerEvent, PSI = share of form 1) is summarized
over all runs in which the PSI is defined, and events whose non-annotated form carries a premature stop are
compared with events whose non-annotated form does not. The test is a two-sided Mann-Whitney U test on the
per-event median shares, by event type, with retained introns kept apart because intron-retaining transcripts
largely escape the decay pathway in Arabidopsis.

Outputs (<out-prefix>_...)
  events.tsv    per event: non-annotated form, premature stop in it, runs with PSI, median and mean share
  summary.tsv   per event type and group: n, median of the per-event median shares, share of events in which
                the non-annotated form is the minor form (median share below 0.5), Mann-Whitney p

Usage
  python ptc_form_abundance.py --nmd results_full/summary/events_high_confidence_nmd.tsv \
         --consequences results_full/summary/events_high_confidence_consequences.tsv \
         --psi results_full/suppa/psi/all_events.psi --out-prefix results_full/summary/ptc_form_abundance
"""
import argparse
import csv
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intron_features import mann_whitney  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmd", required=True, help="events_high_confidence_nmd.tsv")
    ap.add_argument("--consequences", required=True, help="events_high_confidence_consequences.tsv (cds_context)")
    ap.add_argument("--psi", required=True, help="SUPPA2 psiPerEvent table of the event universe")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-runs", type=int, default=20, help="runs with a defined PSI needed to keep an event")
    a = ap.parse_args()

    with open(a.consequences, encoding="utf-8") as fh:
        context = {r["suppa_event"]: r["cds_context"] for r in csv.DictReader(fh, delimiter="\t")}
    wanted = {}                                        # event -> (type, non-annotated form, ptc in it)
    with open(a.nmd, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["annotated_form"] not in ("1", "2") or r["ptc_form"] not in ("1", "2", "none"):
                continue
            if context.get(r["suppa_event"]) != "CDS":
                continue
            alt = "2" if r["annotated_form"] == "1" else "1"
            wanted[r["suppa_event"]] = (r["event_type"], alt, "PTC" if r["ptc_form"] == alt else "no PTC")
    print(f"coding events with both forms reconstructed: {len(wanted)}")

    rows = []
    with open(a.psi, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        n_runs = len(header) - 1 if header[0] == "" else len(header)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ev = f[0]
            if ev not in wanted:
                continue
            vals = []
            for x in f[1:]:
                try:
                    v = float(x)
                except ValueError:
                    continue
                if not math.isnan(v):
                    vals.append(v)
            if len(vals) < a.min_runs:
                continue
            etype, alt, ptc = wanted[ev]
            share = [v if alt == "1" else 1.0 - v for v in vals]
            rows.append(dict(suppa_event=ev, event_type=etype, non_annotated_form=alt, ptc_in_non_annotated=ptc,
                             runs=len(vals), median_share=round(statistics.median(share), 4),
                             mean_share=round(statistics.fmean(share), 4)))
    print(f"PSI table: {n_runs} runs; events kept: {len(rows)}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    with open(a.out_prefix + "_events.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    out = []
    groups = [("SE", ["SE"]), ("A3", ["A3"]), ("A5", ["A5"]), ("MX", ["MX"]), ("SE+A3+A5+MX", ["SE", "A3", "A5", "MX"]),
              ("RI", ["RI"])]
    for label, types in groups:
        sub = [r for r in rows if r["event_type"] in types]
        ptc = [r["median_share"] for r in sub if r["ptc_in_non_annotated"] == "PTC"]
        non = [r["median_share"] for r in sub if r["ptc_in_non_annotated"] == "no PTC"]
        p = mann_whitney(ptc, non)
        for name, v in (("PTC", ptc), ("no PTC", non)):
            out.append(dict(event_type=label, non_annotated_form=name, n=len(v),
                            median_of_median_share=(round(statistics.median(v), 4) if v else ""),
                            minor_form_fraction=(round(sum(1 for x in v if x < 0.5) / len(v), 4) if v else ""),
                            mann_whitney_p=(f"{p:.3g}" if p == p else "")))
    with open(a.out_prefix + "_summary.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# share of the non-annotated form in the event PSI, median over runs with a defined PSI; "
                 "Mann-Whitney p compares the PTC and no-PTC events of the same type\n")
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for r in out:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
