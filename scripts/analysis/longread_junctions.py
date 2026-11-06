"""Independent long-read evidence for splice junctions and splicing events.

Reads a SAM stream of spliced long-read alignments (minimap2 -ax splice, primary alignments only) from
standard input and collects the introns of every read from the N operations of its CIGAR string. An intron
counts only when the aligned blocks on both sides, up to the neighboring intron or the read end, cover at
least --min-anchor reference bases, which removes most spurious introns of error-prone reads. Junction
coordinates are 1-based first and last intron bases, as in STAR SJ.out.tab and ref_junctions.tsv.

Reported (the read count of a junction counts reads, and subreads of one molecule may be counted more than
once, so presence of a junction in at least one read is the main criterion):
  junctions    the share of annotated junctions, pooled novel junctions (by the number of Illumina runs that
               support them) and junctions of the expressed transcript set that occur in the long reads
  events       for SE, A3, A5 and MX events, the events of which at least one junction occurs in the long
               reads (detected) and, among those, the events of which every junction of both forms occurs
               (both forms observed); for RI only the spliced junction can be observed. The same statistics
               are computed for every set given with --events and for the whole SUPPA2 event universe (--ioe),
               which serves as the reference for the coverage of the long reads.

Usage
  zcat reads.fastq.gz | minimap2 -ax splice -ub --secondary=no genome.fa - | samtools view -F 0x904 -q 1 - | \
    python longread_junctions.py --ref-junctions merged/ref_junctions.tsv --pooled-junctions merged/pooled_junctions.tsv \
      --gtf filtered/final.gtf --ioe suppa/events/all_events.ioe \
      --events main=summary/events_high_confidence_annotated.tsv --out-prefix junction_evidence/longread
"""
import argparse
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from event_forms import parse_event  # noqa: E402

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
TYPES = ("SE", "A3", "A5", "MX", "RI")


def read_introns(cigar, pos, min_anchor):
    """Introns (start, end) of one alignment whose flanking aligned blocks cover >= min_anchor bases."""
    ref = pos
    blocks, introns = [0], []           # aligned reference bases between introns
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op in "M=XD":
            blocks[-1] += n
            ref += n
        elif op == "N":
            introns.append((ref, ref + n - 1))
            blocks.append(0)
            ref += n
    return [iv for i, iv in enumerate(introns) if blocks[i] >= min_anchor and blocks[i + 1] >= min_anchor]


def read_junction_table(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        first = fh.readline().rstrip("\n").split("\t")
        header = first if first[0] == "chrom" else None
        if header is None:
            out[(first[0], int(first[1]), int(first[2]))] = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), int(f[2]))
            out[key] = dict(zip(header, f)) if header else None
    return out


def gtf_junctions(path):
    exons = collections.defaultdict(list)
    chrom = {}
    with open(path, encoding="utf-8") as fh:
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


def event_status(ev, observed):
    """(type, detected, both forms observed) of a SUPPA2 event given the observed junctions."""
    chrom, strand, etype, f1, f2, regions = parse_event(ev)
    to_key = lambda j: (chrom, j[0] + 1, j[1] - 1)       # SUPPA2 pairs are exon end - exon start
    k1, k2 = [to_key(j) for j in f1], [to_key(j) for j in f2]
    seen = [k in observed for k in k1 + k2]
    detected = any(seen)
    if etype == "RI":
        return etype, detected, None
    both = all(k in observed for k in k1) and all(k in observed for k in k2)
    return etype, detected, both


def share(a, b):
    return f"{a} of {b} ({100 * a / b:.1f}%)" if b else f"{a} of 0"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-junctions", required=True)
    ap.add_argument("--pooled-junctions", required=True)
    ap.add_argument("--gtf", required=True, help="expressed transcript set (final.gtf)")
    ap.add_argument("--ioe", required=True, help="SUPPA2 event file of the event universe")
    ap.add_argument("--events", action="append", default=[], help="label=path of an events_high_confidence_annotated.tsv")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-anchor", type=int, default=10)
    a = ap.parse_args()

    reads = spliced = 0
    support = collections.Counter()
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.split("\t", 10)
        reads += 1
        introns = read_introns(f[5], int(f[3]), a.min_anchor)
        if introns:
            spliced += 1
        for s, e in set(introns):
            support[(f[2], s, e)] += 1
    observed = set(support)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    with open(a.out_prefix + "_junctions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("chrom\tintron_start\tintron_end\treads\n")
        for (c, s, e), n in sorted(support.items()):
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")

    ref = read_junction_table(a.ref_junctions)
    pooled = read_junction_table(a.pooled_junctions)
    novel = {k: int(v["n_runs"]) for k, v in pooled.items() if v and v.get("status") == "novel"}
    gtf = gtf_junctions(a.gtf)
    lines = [f"alignments read {reads}, with at least one intron {spliced}; distinct junctions {len(observed)}, "
             f"in at least two reads {sum(1 for n in support.values() if n >= 2)}",
             f"annotated junctions observed: {share(len(observed & set(ref)), len(ref))}",
             f"pooled novel junctions observed: {share(sum(1 for k in novel if k in observed), len(novel))}"]
    for lo, hi in ((9, 19), (20, 49), (50, 99), (100, 10 ** 9)):
        sel = [k for k, n in novel.items() if lo <= n <= hi]
        lines.append(f"  novel, seen in {lo}-{hi if hi < 10 ** 9 else 'max'} runs: {share(sum(1 for k in sel if k in observed), len(sel))}")
    cls = {"annotated": [j for j in gtf if j in ref], "pooled novel": [j for j in gtf if j not in ref and j in novel],
           "neither annotated nor pooled": [j for j in gtf if j not in ref and j not in novel]}
    for k, v in cls.items():
        lines.append(f"expressed-set junctions, {k}: {share(sum(1 for j in v if j in observed), len(v))}")
    lr_only = [j for j in observed if j not in ref and j not in novel and j not in gtf]
    lines.append(f"long-read junctions in none of the Illumina sets: {len(lr_only)}, of them in at least two reads "
                 f"{sum(1 for j in lr_only if support[j] >= 2)}")

    sets = [("event universe", None)] + [tuple(x.split("=", 1)) for x in a.events]
    rows = []
    for label, path in sets:
        if path is None:
            evs = []
            with open(a.ioe, encoding="utf-8") as fh:
                fh.readline()
                for line in fh:
                    ev = line.split("\t")[2]
                    if ev.split(";", 1)[1].split(":", 1)[0] in TYPES:
                        evs.append(ev)
        else:
            with open(path, encoding="utf-8") as fh:
                evs = [r["suppa_event"] for r in csv.DictReader(fh, delimiter="\t")]
        c = collections.defaultdict(collections.Counter)
        for ev in evs:
            etype, det, both = event_status(ev, observed)
            c[etype]["n"] += 1
            c[etype]["detected"] += det
            c[etype]["both"] += bool(both)
            if path is not None:
                rows.append(dict(set=label, suppa_event=ev, event_type=etype, detected="yes" if det else "no",
                                 both_forms_observed="" if both is None else ("yes" if both else "no")))
        tot = collections.Counter()
        for t in TYPES:
            if t != "RI":
                tot.update(c[t])
        lines.append(f"{label}: SE, A3, A5 and MX events {tot['n']}, detected {share(tot['detected'], tot['n'])}, "
                     f"both forms observed among detected {share(tot['both'], tot['detected'])}")
        for t in TYPES:
            if c[t]["n"]:
                extra = "" if t == "RI" else f", both forms among detected {share(c[t]['both'], c[t]['detected'])}"
                lines.append(f"  {t}: {c[t]['n']} events, detected {share(c[t]['detected'], c[t]['n'])}{extra}")
    if rows:
        with open(a.out_prefix + "_events.tsv", "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    with open(a.out_prefix + "_summary.txt", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# minimum anchor on both sides of an intron: {a.min_anchor} aligned reference bases\n")
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
