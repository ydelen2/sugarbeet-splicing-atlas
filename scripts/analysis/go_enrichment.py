"""GO term enrichment of the genes with high-confidence differential splicing events.

Gene annotation: the eggNOG-mapper table of BeetSet-3 (one row per transcript, GO terms in the GOs
column); a gene carries the union of the terms of its transcripts. The eggNOG lists already contain
the ancestor terms of every annotation (the root terms are present in annotated rows), so no
propagation is done here; the .obo file is only used for the term names and namespaces.

Foreground sets: the genes with a high-confidence event in each condition group (the groups column of
events_high_confidence_annotated.tsv, split on commas) and in all groups together. Background: the
reference genes with at least one testable SUPPA2 event in any comparison. Background and foregrounds
are restricted to genes with at least one GO term, and the counts are reported.

Test: hypergeometric upper tail for every term with at least --min-bg background genes and --min-fg
foreground genes, computed in log space; Benjamini-Hochberg over the terms tested in the same
foreground set.

Usage
  python go_enrichment.py --functional ref/BeetSet-3.genes.2604.functional --obo ref/go/go-basic.obo \
      --events results_full/summary/events_high_confidence_annotated.tsv --suppa-dir results_full/suppa \
      --map results_full/filtered/transcript_map.tsv --out results_full/summary/go_enrichment.tsv
"""
import argparse
import collections
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splicing_regulators import testable_genes  # noqa: E402
from event_concordance import bh  # noqa: E402

ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_upper_log(k, n, K, N):
    """P(X >= k) for n draws with K successes among N, summed in log space."""
    denom = log_comb(N, n)
    total = 0.0
    for x in range(k, min(n, K) + 1):
        total += math.exp(log_comb(K, x) + log_comb(N - K, n - x) - denom)
    return min(1.0, total)


def read_go(functional):
    """{gene: set of GO ids}; also returns how many annotated rows carry a root term."""
    gos = collections.defaultdict(set)
    rows, with_root = 0, 0
    with open(functional, encoding="utf-8") as fh:
        header = None
        for line in fh:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r.get("GOs", "-") in ("-", ""):
                continue
            terms = set(r["GOs"].split(","))
            rows += 1
            with_root += bool(terms & ROOTS)
            gos[r["query"].rsplit(".", 1)[0]] |= terms
    return gos, rows, with_root


def read_obo(path):
    names = {}
    if not path or not os.path.isfile(path):
        return names
    cur = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                if cur.get("id"):
                    names[cur["id"]] = (cur.get("name", ""), cur.get("namespace", ""))
                cur = {}
            elif ":" in line and not line.startswith("["):
                key, val = line.split(":", 1)
                if key in ("id", "name", "namespace") and key not in cur:
                    cur[key] = val.strip()
    if cur.get("id"):
        names[cur["id"]] = (cur.get("name", ""), cur.get("namespace", ""))
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functional", required=True)
    ap.add_argument("--obo")
    ap.add_argument("--events", required=True)
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--comparisons", help="comparison table that defines the testable-gene background")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-bg", type=int, default=5)
    ap.add_argument("--min-fg", type=int, default=2)
    ap.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args()

    gos, rows, with_root = read_go(a.functional)
    names = read_obo(a.obo)
    print(f"annotated transcripts: {rows}, with a root term: {with_root}; genes with GO terms: {len(gos)}; term names: {len(names)}")

    bg_all = testable_genes(a.suppa_dir, a.map, a.comparisons)
    bg = {g for g in bg_all if g in gos}
    fg = collections.defaultdict(set)
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = [g for g in r["ref_gene"].split(",") if g]
            for grp in r["groups"].split(",") + ["all"]:
                fg[grp].update(genes)
    print(f"background: {len(bg_all)} testable genes, {len(bg)} with GO terms")

    term_bg = collections.defaultdict(set)
    for g in bg:
        for t in gos[g]:
            term_bg[t].add(g)

    out = []
    for grp in sorted(fg):
        f = {g for g in fg[grp] if g in bg}
        term_fg = collections.defaultdict(set)
        for g in f:
            for t in gos[g]:
                term_fg[t].add(g)
        tested = [t for t in term_bg if len(term_bg[t]) >= a.min_bg and len(term_fg.get(t, ())) >= a.min_fg]
        pvals = [hypergeom_upper_log(len(term_fg[t]), len(f), len(term_bg[t]), len(bg)) for t in tested]
        adj = bh(pvals)
        sig = 0
        for t, p, q in zip(tested, pvals, adj):
            if q < a.fdr:
                sig += 1
            out.append(dict(set=grp, go=t, name=names.get(t, ("", ""))[0], namespace=names.get(t, ("", ""))[1],
                            foreground_genes=len(f), background_genes=len(bg), term_foreground=len(term_fg[t]),
                            term_background=len(term_bg[t]), expected=round(len(f) * len(term_bg[t]) / len(bg), 2),
                            fold=round(len(term_fg[t]) * len(bg) / (len(f) * len(term_bg[t])), 2), p=p, fdr=q))
        print(f"{grp}: {len(fg[grp])} genes with events, {len(f)} in the GO background; terms tested {len(tested)}, "
              f"FDR < {a.fdr}: {sig}")

    out.sort(key=lambda d: (d["set"], d["fdr"], d["p"]))
    cols = list(out[0].keys()) if out else []
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
