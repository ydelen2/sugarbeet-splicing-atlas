"""Annotate the high-confidence events with reference gene models and eggNOG descriptions.

Input events: event_concordance_high_confidence.tsv (one row per event and comparison, any event type).
Output: one row per event with
  event_type      SE, A3, A5, RI or MX
  ref_gene        the BeetSet-3 gene whose span contains all coordinates of the event on the same strand
                  (several candidates are kept comma separated when more than one gene qualifies)
  description, preferred_name, pfams   from the eggNOG-mapper table, taken from the transcript of
                  that gene with the highest score; "-" values in the table are written as empty
  n_comparisons, comparisons, groups   where the event was found (group = comparison id up to the first "_",
                  with the temperature suffix of the cold comparisons removed, so cold4C and cold0C are one group)
  direction       "higher" or "lower" inclusion in the test group when every comparison agrees, else "mixed"
  suppa_dpsi_mean, rmats_dpsi_mean

Usage
  python annotate_events.py --events results_full/summary/event_concordance_high_confidence.tsv \
         --gff ref/BeetSet-3.genes.2604.gff --functional ref/BeetSet-3.genes.2604.functional \
         --out results_full/summary/events_high_confidence_annotated.tsv
"""
import argparse
import collections
import csv
import re


def read_genes(gff):
    """{gene_id: (chrom, strand, start, end)} from gene lines."""
    genes = {}
    with open(gff, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
            if m:
                genes[m.group(1)] = (f[0], f[6], int(f[3]), int(f[4]))
    return genes


def read_functional(path):
    """{gene_id: (score, description, preferred_name, pfams)} keeping the best scoring transcript."""
    best = {}
    with open(path, encoding="utf-8") as fh:
        header = None
        for line in fh:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            gene = r["query"].rsplit(".", 1)[0]
            try:
                score = float(r["score"])
            except ValueError:
                continue
            if gene not in best or score > best[gene][0]:
                best[gene] = (score, *(("" if r[c] == "-" else r[c]) for c in ("Description", "Preferred_name", "PFAMs")))
    return best


def parse_event(ev):
    """Sequence, strand and the span of all coordinates of a SUPPA2 event id (any event type)."""
    gene, body = ev.split(";", 1)
    parts = body.split(":")
    chrom, strand = parts[1], parts[-1]
    coords = [int(x) for f in parts[2:-1] for x in f.split("-")]
    return chrom, strand, min(coords), max(coords)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--functional", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes = read_genes(a.gff)
    func = read_functional(a.functional)

    rows = collections.defaultdict(list)
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["suppa_event"]].append(r)

    out = []
    for ev, rs in rows.items():
        chrom, strand, lo, hi = parse_event(ev)
        cands = [g for g in rs[0]["ref_genes"].split(",") if g]
        hits = [g for g in cands if g in genes and genes[g][0] == chrom and genes[g][1] == strand
                and genes[g][2] <= lo and hi <= genes[g][3]]
        if not hits:
            hits = cands
        ann = [func.get(g, (0, "", "", "")) for g in hits]
        signs = {("higher" if float(r["suppa_dpsi"]) > 0 else "lower") for r in rs}
        comps = sorted(r["comparison"] for r in rs)
        out.append(dict(
            suppa_event=ev, event_type=rs[0].get("event_type", ev.split(";", 1)[1].split(":")[0]),
            merged_gene=rs[0]["gene"], ref_gene=",".join(hits),
            ref_gene_candidates=",".join(cands),
            description=" | ".join(x[1] for x in ann), preferred_name=" | ".join(x[2] for x in ann),
            pfams=" | ".join(x[3] for x in ann),
            n_comparisons=len(rs), comparisons=",".join(comps),
            groups=",".join(sorted({re.sub(r"\d+C$", "", c.split("_")[0]) for c in comps})),
            direction=signs.pop() if len(signs) == 1 else "mixed",
            suppa_dpsi_mean=round(sum(float(r["suppa_dpsi"]) for r in rs) / len(rs), 4),
            rmats_dpsi_mean=round(sum(float(r["rmats_dpsi"]) for r in rs) / len(rs), 4)))

    out.sort(key=lambda d: (-d["n_comparisons"], d["suppa_event"]))
    cols = list(out[0].keys()) if out else []
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    multi = sum(1 for d in out if "," in d["ref_gene"])
    none = sum(1 for d in out if not d["ref_gene"])
    print(f"events: {len(out)}; single reference gene: {len(out) - multi - none}; several: {multi}; none: {none}")
    print(f"events without eggNOG description: {sum(1 for d in out if not d['description'])}")
    by_t = collections.Counter(d["event_type"] for d in out)
    print("events by type: " + ", ".join(f"{k}: {by_t[k]}" for k in sorted(by_t)))
    by_n = collections.Counter(d["n_comparisons"] for d in out)
    print("events by number of comparisons: " + ", ".join(f"{k}: {by_n[k]}" for k in sorted(by_n)))
    print(f"mixed direction across comparisons: {sum(1 for d in out if d['direction'] == 'mixed')}")
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
