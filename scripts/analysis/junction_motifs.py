"""Splice-site motifs and intron lengths of the novel junctions, and the junction support of the events.

STAR removes non-canonical unannotated junctions during alignment, so every novel junction is canonical by
construction; the informative comparison is the share of the three canonical classes (GT-AG, GC-AG, AT-AC)
and the intron length, between annotated junctions, all pooled novel junctions (by the number of runs that
support them) and the novel junctions used by the expressed transcript set. For the high-confidence events
the script reports whether every junction of an event is annotated, annotated or pooled, or whether at least
one junction lies below the pooling threshold.

Inputs may be plain or gzip-compressed. STAR motif codes: 0 non-canonical, 1 and 2 GT-AG, 3 and 4 GC-AG,
5 and 6 AT-AC.

Usage
  python junction_motifs.py --pooled-junctions merged/pooled_junctions.tsv --ref-junctions merged/ref_junctions.tsv \
         --gtf filtered/final.gtf --events main=summary/events_high_confidence_annotated.tsv \
         --out junction_evidence/junction_motifs.txt
"""
import argparse
import collections
import csv
import gzip
import os
import re
import statistics

MOTIF = {"0": "non-canonical", "1": "GT-AG", "2": "GT-AG", "3": "GC-AG", "4": "GC-AG", "5": "AT-AC", "6": "AT-AC"}
CLASSES = ("GT-AG", "GC-AG", "AT-AC", "non-canonical")


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def describe(rows, label):
    n = len(rows)
    c = collections.Counter(MOTIF[r["motif_code"]] for r in rows)
    lengths = [int(r["intron_end"]) - int(r["intron_start"]) + 1 for r in rows]
    return (f"{label}: {n} junctions; " + ", ".join(f"{k} {c[k]} ({100 * c[k] / n:.2f}%)" for k in CLASSES)
            + f"; intron length median {statistics.median(lengths):.0f} nt, at most 100 nt "
              f"{100 * sum(1 for x in lengths if x <= 100) / n:.1f}%, longer than 10 kb {100 * sum(1 for x in lengths if x > 10000) / n:.2f}%")


def gtf_junctions(path):
    exons, chrom = collections.defaultdict(list), {}
    with opener(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "exon":
                t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                exons[t].append((int(f[3]), int(f[4])))
                chrom[t] = f[0]
    js = set()
    for t, ex in exons.items():
        ex.sort()
        for i in range(len(ex) - 1):
            js.add((chrom[t], ex[i][1] + 1, ex[i + 1][0] - 1))
    return js


def event_junctions(ev):
    """Junctions (chrom, first intron base, last intron base) of a SUPPA2 event id."""
    parts = ev.split(";", 1)[1].split(":")
    chrom = parts[1]
    return [(chrom, int(x) + 1, int(y) - 1) for f in parts[2:-1] if "-" in f for x, y in [f.split("-")]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pooled-junctions", required=True)
    ap.add_argument("--ref-junctions", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--events", action="append", default=[], help="label=path of an events_high_confidence_annotated.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with opener(a.pooled_junctions) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    pooled = {(r["chrom"], int(r["intron_start"]), int(r["intron_end"])): r for r in rows}
    with opener(a.ref_junctions) as fh:
        ref = {(f[0], int(f[1]), int(f[2])) for f in (l.rstrip("\n").split("\t") for l in fh)}
    lines = [describe([r for r in rows if r["status"] == "annotated"], "annotated junctions (pooled table)"),
             describe([r for r in rows if r["status"] == "novel"], "novel junctions (pooled table)")]
    novel = [r for r in rows if r["status"] == "novel"]
    for lo, hi in ((9, 19), (20, 49), (50, 99), (100, 10 ** 9)):
        sel = [r for r in novel if lo <= int(r["n_runs"]) <= hi]
        if sel:
            lines.append("  " + describe(sel, f"novel, seen in {lo}-{hi if hi < 10 ** 9 else 'max'} runs"))
    js = gtf_junctions(a.gtf)
    used = [pooled[j] for j in js if j not in ref and j in pooled]
    lines.append(f"expressed-set junctions: {len(js)}; annotated {sum(1 for j in js if j in ref)}, pooled novel {len(used)}, "
                 f"neither annotated nor pooled {sum(1 for j in js if j not in ref and j not in pooled)}")
    lines.append("  " + describe(used, "pooled novel junctions used by the expressed set"))
    lines.append(f"  median number of supporting runs of these junctions: {statistics.median(int(r['n_runs']) for r in used):.0f}")
    for item in a.events:
        label, path = item.split("=", 1)
        with open(path, encoding="utf-8") as fh:
            evs = list(csv.DictReader(fh, delimiter="\t"))
        c = collections.Counter()
        for r in evs:
            ej = event_junctions(r["suppa_event"])
            if all(j in ref for j in ej):
                c["all junctions annotated"] += 1
            elif all(j in ref or j in pooled for j in ej):
                c["all junctions annotated or pooled"] += 1
            else:
                c["at least one junction below the pooling threshold"] += 1
        lines.append(f"{label} high-confidence events: {len(evs)}; " + ", ".join(f"{k} {v}" for k, v in c.items()))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
