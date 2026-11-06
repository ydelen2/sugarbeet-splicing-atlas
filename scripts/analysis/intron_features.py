"""Length and GC content of the alternative regions of the high-confidence events against the annotation.

Three sets are compared for retained introns and for skipped exons:
  annotated   all introns (GFF intron features) or all internal exons (exons of an mRNA except its first and
              last) of BeetSet-3
  tested      the alternative region of every event of that type in the SUPPA2 event file (the event universe)
  consensus   the alternative region of the high-confidence events of that type
For every region the length and the GC fraction of the genomic sequence are computed, and the consensus set
is compared with the other two sets by a two-sided Mann-Whitney U test (normal approximation with tie
correction).

Usage
  python intron_features.py --gff ref/BeetSet-3.genes.2604.gff --genome ref/RefBeet-3.0.fa \
      --ioe results_full/suppa/events/all_events.ioe --events results_full/summary/events_high_confidence_annotated.tsv \
      --out-prefix results_full/summary/region_features
Outputs <out-prefix>_regions.tsv (one row per region) and <out-prefix>_summary.tsv.
"""
import argparse
import collections
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from event_consequences import Fasta, regions  # noqa: E402


def gff_regions(path):
    """Annotated introns and internal exons: [(chrom, start, end, strand)]."""
    introns, exons_by_mrna = [], collections.defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "intron":
                introns.append((f[0], int(f[3]), int(f[4]), f[6]))
            elif f[2] == "exon":
                parent = next((kv.split("=", 1)[1] for kv in f[8].split(";") if kv.startswith("Parent=")), "")
                exons_by_mrna[parent].append((f[0], int(f[3]), int(f[4]), f[6]))
    internal = set()
    for ex in exons_by_mrna.values():
        ex.sort(key=lambda e: e[1])
        for e in ex[1:-1]:
            internal.add(e)
    return sorted(set(introns)), sorted(internal)


def event_regions(ids):
    out = []
    for ev in ids:
        chrom, strand, etype, regs = regions(ev)
        if etype in ("RI", "SE") and regs:
            out.append((etype, (chrom, regs[0][0], regs[0][1], strand)))
    return out


def mann_whitney(x, y):
    """Two-sided p-value of the Mann-Whitney U test, normal approximation with tie and continuity corrections."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return float("nan")
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks, i, tie_term = [0.0] * len(allv), 0, 0.0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j + 2) / 2
        for k in range(i, j + 1):
            ranks[k] = r
        t = j - i + 1
        tie_term += t ** 3 - t
        i = j + 1
    r1 = sum(r for r, (v, g) in zip(ranks, allv) if g == 0)
    u1 = r1 - n1 * (n1 + 1) / 2
    n = n1 + n2
    mean = n1 * n2 / 2
    var = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1)))
    if var == 0:
        return 1.0
    z = (abs(u1 - mean) - 0.5) / math.sqrt(var)          # continuity correction, as in scipy's asymptotic method
    return math.erfc(max(z, 0.0) / math.sqrt(2))


def median(v):
    s = sorted(v)
    n = len(s)
    return (s[n // 2 - 1] + s[n // 2]) / 2 if n % 2 == 0 else s[n // 2]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--ioe", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    fasta = Fasta(a.genome)
    introns, exons = gff_regions(a.gff)
    with open(a.ioe, encoding="utf-8") as fh:
        fh.readline()
        tested = event_regions(line.split("\t")[2] for line in fh)
    with open(a.events, encoding="utf-8") as fh:
        consensus = event_regions(r["suppa_event"] for r in csv.DictReader(fh, delimiter="\t"))

    sets = {("RI", "annotated"): introns, ("SE", "annotated"): exons,
            ("RI", "tested"): [r for t, r in tested if t == "RI"], ("SE", "tested"): [r for t, r in tested if t == "SE"],
            ("RI", "consensus"): [r for t, r in consensus if t == "RI"], ("SE", "consensus"): [r for t, r in consensus if t == "SE"]}
    feats = {}
    with open(a.out_prefix + "_regions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("type\tset\tchrom\tstart\tend\tstrand\tlength\tgc\n")
        for (etype, name), regs in sets.items():
            vals = []
            for chrom, s, e, strand in regs:
                if chrom not in fasta.idx:
                    continue
                seq = fasta.fetch(chrom, s, e)
                if not seq:
                    continue
                gc = sum(1 for c in seq if c in "GC") / len(seq)
                vals.append((e - s + 1, gc))
                fh.write(f"{etype}\t{name}\t{chrom}\t{s}\t{e}\t{strand}\t{e - s + 1}\t{gc:.4f}\n")
            feats[(etype, name)] = vals
    rows = []
    for etype in ("RI", "SE"):
        cons = feats[(etype, "consensus")]
        for name in ("annotated", "tested", "consensus"):
            v = feats[(etype, name)]
            row = dict(type=etype, set=name, n=len(v), median_length=median([x[0] for x in v]) if v else "",
                       median_gc=round(median([x[1] for x in v]), 4) if v else "")
            if name != "consensus" and v and cons:
                row["p_length_vs_consensus"] = f"{mann_whitney([x[0] for x in cons], [x[0] for x in v]):.3g}"
                row["p_gc_vs_consensus"] = f"{mann_whitney([x[1] for x in cons], [x[1] for x in v]):.3g}"
            else:
                row["p_length_vs_consensus"] = row["p_gc_vs_consensus"] = ""
            rows.append(row)
    with open(a.out_prefix + "_summary.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(r[k]) for k in rows[0]))
    print(f"tables: {a.out_prefix}_regions.tsv  {a.out_prefix}_summary.tsv")


if __name__ == "__main__":
    main()
