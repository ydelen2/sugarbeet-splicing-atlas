"""Enrichment of high-confidence events in gene sets, controlled for the number of testable events per gene.

The hypergeometric tests of splicing_regulators.py and go_enrichment.py treat every testable gene as having the
same chance of carrying a high-confidence event. A gene with many testable events has more chances to be called,
and the splicing-regulator genes are expressed in most tissues and carry several events each. This script repeats
the test as a permutation: for every gene in the background (reference genes with at least one testable SUPPA2
event in the comparisons of the comparison table) the number of distinct testable events over all comparisons is
counted; genes are binned by that number (1, 2, 3-4, 5-8, 9-16, 17 or more); random gene sets are drawn with the
same number of genes per bin as the set under test; and the share of draws in which at least as many genes carry
a high-confidence event as in the real set is the permutation p-value.

Gene sets tested: the SR genes and all family genes from splicing_regulators_families.tsv (column das_events
non-empty marks the genes with events), the genes annotated with a GO term (--functional with --go-term, for
example GO:0008380 RNA splicing, read from the eggNOG-mapper table as in go_enrichment.py), and any sets given
with --gene-set name=file (one gene per line). The GO sets are restricted to the background genes that carry at
least one GO term, as in the enrichment test.

Outputs: <out> with one line per set: genes, genes with events, median testable events per gene in the set and in
the background, expected genes with events under the matched null, permutation p.

Usage
  python matched_permutation.py --suppa-dir results_full/suppa --comparisons metadata/comparisons_all.tsv \
         --map results_full/filtered/transcript_map.tsv --events results_full/summary/events_high_confidence_annotated.tsv \
         --families results_full/domains/splicing_regulators_families.tsv \
         --functional ref/BeetSet-3.genes.2604.functional --go-term GO:0008380 --go-term GO:0000398 \
         --out results_full/summary/matched_permutation.tsv
"""
import argparse
import collections
import csv
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from go_enrichment import read_go  # noqa: E402
from splicing_regulators import comparison_ids  # noqa: E402

BINS = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, None)]


def bin_of(n):
    for i, (lo, hi) in enumerate(BINS):
        if n >= lo and (hi is None or n <= hi):
            return i
    return len(BINS) - 1


def testable_events_per_gene(suppa_dir, map_path, comparisons):
    ref_of = collections.defaultdict(set)
    with open(map_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ref_gene"]:
                ref_of[r["gene"]].add(r["ref_gene"])
    events_of = collections.defaultdict(set)
    for c in comparison_ids(comparisons):
        path = os.path.join(suppa_dir, "diff", c, "events.dpsi")
        with open(path, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                ev, d, p = line.rstrip("\n").split("\t")[:3]
                if d != "nan" and p != "nan":
                    for g in ev.split(";", 1)[0].split("_and_"):
                        for rg in ref_of.get(g, ()):
                            events_of[rg].add(ev)
    return {g: len(v) for g, v in events_of.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--events", required=True, help="events_high_confidence_annotated.tsv")
    ap.add_argument("--families", required=True, help="splicing_regulators_families.tsv")
    ap.add_argument("--gene-set", action="append", default=[], help="name=file with one gene per line")
    ap.add_argument("--functional", help="eggNOG-mapper table of BeetSet-3 (GO terms per transcript)")
    ap.add_argument("--go-term", action="append", default=[], help="GO id whose annotated genes form a set")
    ap.add_argument("--draws", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    n_test = testable_events_per_gene(a.suppa_dir, a.map, a.comparisons)
    with_event = set()
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for g in r["ref_gene"].split(","):
                if g:
                    with_event.add(g)
    background = set(n_test)
    print(f"background {len(background)} genes, {len(background & with_event)} with a high-confidence event")

    sets = {}
    with open(a.families, encoding="utf-8") as fh:
        fam = list(csv.DictReader(fh, delimiter="\t"))
    sets["SR"] = {r["gene"] for r in fam if r["family"] == "SR"}
    sets["all families"] = {r["gene"] for r in fam}
    for item in a.gene_set:
        name, path = item.split("=", 1)
        with open(path, encoding="utf-8") as fh:
            sets[name] = {l.strip() for l in fh if l.strip()}

    go_bg = None
    if a.go_term:
        gos, _, _ = read_go(a.functional)
        go_bg = {g for g in background if gos.get(g)}
        for term in a.go_term:
            sets[term] = {g for g in go_bg if term in gos[g]}
        print(f"GO background: {len(go_bg)} genes with a GO term")
    by_bin = collections.defaultdict(list)
    for g in background:
        by_bin[bin_of(n_test[g])].append(g)
    by_bin_go = collections.defaultdict(list)
    for g in (go_bg or ()):
        by_bin_go[bin_of(n_test[g])].append(g)
    rng = random.Random(a.seed)
    rows = []
    bg_median = statistics.median(n_test[g] for g in background)
    for name, genes in sets.items():
        pool = by_bin_go if name in a.go_term else by_bin
        genes = genes & (go_bg if name in a.go_term else background)
        k = len(genes & with_event)
        need = collections.Counter(bin_of(n_test[g]) for g in genes)
        hits, total = 0, 0
        for _ in range(a.draws):
            draw = 0
            for b, m in need.items():
                draw += sum(1 for g in rng.sample(pool[b], m) if g in with_event)
            total += draw
            hits += draw >= k
        rows.append(dict(set=name, genes=len(genes), with_event=k,
                         median_testable_events=statistics.median(n_test[g] for g in genes) if genes else "",
                         background_median_testable_events=bg_median,
                         expected_matched=round(total / a.draws, 2), permutation_p=f"{(hits + 1) / (a.draws + 1):.2g}",
                         bins=";".join(f"{BINS[b][0]}-{BINS[b][1] or 'max'}:{m}" for b, m in sorted(need.items()))))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# {a.draws} random gene sets matched on the number of testable events per gene (bins 1, 2, 3-4, 5-8, "
                 f"9-16, 17+); permutation p = (draws with at least as many genes with events + 1) / (draws + 1)\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(open(a.out, encoding="utf-8").read())


if __name__ == "__main__":
    main()
