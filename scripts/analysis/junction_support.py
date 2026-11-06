"""Illumina support of the junctions of the expressed transcript set that are neither annotated nor pooled.

The pooled junction table keeps novel junctions with at least --min-reads uniquely mapped reads in at least
2 percent of the merged runs. Transcripts of the expressed set can still use junctions below that threshold,
because the assembler works per run. For every such junction this script reads the STAR SJ.out.tab of every
aligned run and reports the splice-site motif, the number of runs with at least one and with at least
--min-reads uniquely mapped reads, and the total number of uniquely mapped reads.

STAR motif codes: 0 non-canonical, 1 and 2 GT-AG, 3 and 4 GC-AG, 5 and 6 AT-AC (odd codes plus strand).

Usage
  python junction_support.py --gtf filtered/final.gtf --ref-junctions merged/ref_junctions.tsv \
         --pooled-junctions merged/pooled_junctions.tsv --runres-dir per_run --sheet metadata/samples_all.tsv \
         --out-prefix junction_evidence/support
"""
import argparse
import collections
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from longread_junctions import gtf_junctions, read_junction_table  # noqa: E402

MOTIF = {"0": "non-canonical", "1": "GT-AG", "2": "GT-AG", "3": "GC-AG", "4": "GC-AG", "5": "AT-AC", "6": "AT-AC"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--ref-junctions", required=True)
    ap.add_argument("--pooled-junctions", required=True)
    ap.add_argument("--runres-dir", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-reads", type=int, default=3)
    a = ap.parse_args()

    ref = read_junction_table(a.ref_junctions)
    pooled = read_junction_table(a.pooled_junctions)
    target = {j for j in gtf_junctions(a.gtf) if j not in ref and j not in pooled}
    with open(a.sheet, encoding="utf-8") as fh:
        runs = [r["run"] for r in csv.DictReader(fh, delimiter="\t") if r.get("use_das") != "excluded"]
    motif = {}
    runs_any, runs_min, reads = collections.Counter(), collections.Counter(), collections.Counter()
    missing = []
    for run in runs:
        path = os.path.join(a.runres_dir, run, f"{run}.SJ.out.tab")
        if not os.path.isfile(path):
            missing.append(run)
            continue
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                f = line.split("\t")
                key = (f[0], int(f[1]), int(f[2]))
                if key in target:
                    u = int(f[6])
                    motif[key] = MOTIF.get(f[4], f[4])
                    if u >= 1:
                        runs_any[key] += 1
                    if u >= a.min_reads:
                        runs_min[key] += 1
                    reads[key] += u
    if missing:
        print(f"SJ.out.tab missing for {len(missing)} runs, first: {missing[:3]}")

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    with open(a.out_prefix + "_junctions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"chrom\tintron_start\tintron_end\tmotif\truns_with_reads\truns_with_{a.min_reads}_reads\tunique_reads\n")
        for k in sorted(target):
            fh.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{motif.get(k, 'not in SJ.out.tab')}\t{runs_any[k]}\t{runs_min[k]}\t{reads[k]}\n")
    n = len(target)
    m = collections.Counter(motif.get(k, "not in SJ.out.tab") for k in target)
    lines = [f"runs read {len(runs) - len(missing)}; expressed-set junctions neither annotated nor pooled: {n}",
             "motif: " + ", ".join(f"{k} {v} ({100 * v / n:.2f}%)" for k, v in m.most_common())]
    ra = [runs_any[k] for k in target]
    lines.append(f"runs with at least one unique read: median {statistics.median(ra) if ra else 0}, "
                 f"0 runs {sum(1 for x in ra if x == 0)}, 1 run {sum(1 for x in ra if x == 1)}, "
                 f"2-8 runs {sum(1 for x in ra if 2 <= x <= 8)}, 9 or more runs {sum(1 for x in ra if x >= 9)}")
    rm = [runs_min[k] for k in target]
    lines.append(f"runs with at least {a.min_reads} unique reads: 0 runs {sum(1 for x in rm if x == 0)}, "
                 f"1-8 runs {sum(1 for x in rm if 1 <= x <= 8)}, 9 or more {sum(1 for x in rm if x >= 9)}")
    tr = sorted(reads[k] for k in target)
    if tr:
        lines.append(f"unique reads summed over runs: median {statistics.median(tr)}, 10th percentile {tr[len(tr) // 10]}")
    with open(a.out_prefix + "_summary.txt", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
