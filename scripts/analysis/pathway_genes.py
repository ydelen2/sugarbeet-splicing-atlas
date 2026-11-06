"""Genes of one pathway among the genes with high-confidence splicing events, by annotation terms.

The pathway is defined by a term table (metadata/curation/aba_pathway_terms.tsv for abscisic acid) with
the columns field (description or preferred_name), pattern (case-insensitive regular expression applied to
the whole eggNOG field), component and note. A gene belongs to the pathway when any pattern matches its
eggNOG description or preferred name in the BeetSet-3 functional annotation (best scoring transcript per
gene, as in annotate_events.py).

Reported: every pathway gene with a high-confidence event (events, types, comparisons, groups, direction),
the number of pathway genes among all annotated genes and among the genes with a testable SUPPA2 event
(background), and a hypergeometric test of the pathway share among the genes with events against that
background, overall and per condition group.

Usage
  python pathway_genes.py --functional ref/BeetSet-3.genes.2604.functional --terms metadata/curation/aba_pathway_terms.tsv \
         --events results_full/summary/events_high_confidence_annotated.tsv --suppa-dir results_full/suppa \
         --map results_full/filtered/transcript_map.tsv --out-prefix results_full/summary/aba_pathway
"""
import argparse
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_events import read_functional  # noqa: E402
from splicing_regulators import testable_genes  # noqa: E402
from stats_util import hypergeom_upper  # noqa: E402


def read_terms(path):
    terms = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            terms.append((r["field"], re.compile(r["pattern"], re.I), r["component"]))
    return terms


def match(terms, description, preferred):
    hits = []
    for field, pat, comp in terms:
        text = description if field == "description" else preferred
        if pat.search(text or ""):
            hits.append((comp, pat.pattern))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functional", required=True)
    ap.add_argument("--terms", required=True)
    ap.add_argument("--events", required=True, help="events_high_confidence_annotated.tsv")
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--comparisons", help="comparison table that defines the testable-gene background")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    terms = read_terms(a.terms)
    func = read_functional(a.functional)
    pathway = {}
    for gene, (_, desc, pref, _) in func.items():
        hits = match(terms, desc, pref)
        if hits:
            pathway[gene] = (desc, pref, hits)
    bg = testable_genes(a.suppa_dir, a.map, a.comparisons)
    print(f"annotated genes {len(func)}, pathway genes {len(pathway)}; testable background {len(bg)}, "
          f"pathway genes in it {len(bg & set(pathway))}")

    events = collections.defaultdict(list)
    fg_all, fg_group = set(), collections.defaultdict(set)
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for g in r["ref_gene"].split(","):
                if not g:
                    continue
                fg_all.add(g)
                for grp in r["groups"].split(","):
                    fg_group[grp].add(g)
                if g in pathway:
                    events[g].append(r)

    rows = []
    for g in sorted(events):
        desc, pref, hits = pathway[g]
        for r in events[g]:
            rows.append(dict(
                gene=g, preferred_name=pref, description=desc,
                component=";".join(sorted({c for c, _ in hits})), matched_pattern=";".join(p for _, p in hits),
                suppa_event=r["suppa_event"], event_type=r["event_type"], n_comparisons=r["n_comparisons"],
                comparisons=r["comparisons"], groups=r["groups"], direction=r["direction"],
                suppa_dpsi_mean=r["suppa_dpsi_mean"], rmats_dpsi_mean=r["rmats_dpsi_mean"]))
    with open(a.out_prefix + "_events.tsv", "w", encoding="utf-8", newline="\n") as fh:
        cols = list(rows[0].keys()) if rows else ["gene"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    tests = []
    bg_path = bg & set(pathway)
    for name, fg in [("all", fg_all)] + sorted(fg_group.items()):
        f = fg & bg
        k = len(f & bg_path)
        expected = len(f) * len(bg_path) / len(bg) if bg else 0.0
        tests.append(dict(set=name, genes_with_events=len(f), pathway_genes_with_events=k,
                          background_genes=len(bg), pathway_genes_in_background=len(bg_path),
                          expected=round(expected, 2), enrichment=("" if expected == 0 else round(k / expected, 2)),
                          p_hypergeometric=f"{hypergeom_upper(k, len(f), len(bg_path), len(bg)):.3g}"))
    with open(a.out_prefix + "_enrichment.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tests[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tests)

    print(f"pathway genes with high-confidence events: {len(events)} genes, {len(rows)} events")
    for t in tests:
        print(f"  {t['set']}: {t['pathway_genes_with_events']} of {t['genes_with_events']} genes with events, "
              f"expected {t['expected']}, p {t['p_hypergeometric']}")
    print(f"tables: {a.out_prefix}_events.tsv  {a.out_prefix}_enrichment.tsv")


if __name__ == "__main__":
    main()
