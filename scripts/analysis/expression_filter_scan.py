"""How many transcripts would an expression filter keep, for several thresholds?

Reference transcripts (class "=") are always kept by filter_transcripts.py, so only the other
classes are counted here. Two kinds of criteria are evaluated on the transcript TPM matrix:

  runs   TPM >= t in at least n runs of the whole data set
  group  TPM >= t in every run of at least one comparison group (a reference or test group of any
         comparison in the comparison table). This ties the filter to the design: a transcript is
         kept when it is consistently expressed in at least one biological condition that is tested.

Usage
  python expression_filter_scan.py --tpm tpm_all.tsv --map transcript_map.tsv --sheet samples_all.tsv \
         --comparisons comparisons_all.tsv --out expression_filter_scan.tsv
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_splicing_inputs import comparison_groups, read_tsv  # noqa: E402

TPMS = [0.5, 1.0, 2.0, 5.0]
MIN_RUNS = [2, 3, 5, 10, 25, 50]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    code = {}
    with open(a.map, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kept_structural"] == "yes":
                code[r["transcript"]] = r["class_code"]

    sheet = read_tsv(a.sheet)
    groups = []
    for comp in read_tsv(a.comparisons):
        ref, test = comparison_groups(sheet, comp)
        groups += [ref, test]

    with open(a.tpm, encoding="utf-8") as fh:
        samples = fh.readline().rstrip("\n").split("\t")
        col = {s: i for i, s in enumerate(samples)}
        gidx = []
        for g in groups:
            missing = [r for r in g if r not in col]
            if missing:
                sys.exit(f"runs of a comparison group are missing from the TPM matrix: {missing[:5]}")
            gidx.append([col[r] for r in g])

        runs_hist = {t: collections.defaultdict(collections.Counter) for t in TPMS}  # t -> class -> n_runs -> count
        group_pass = {t: collections.Counter() for t in TPMS}                     # t -> class -> count
        total = collections.Counter()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c = code.get(f[0])
            if c is None or c == "=":
                continue
            vals = [float(v) for v in f[1:]]
            total[c] += 1
            for t in TPMS:
                runs_hist[t][c][sum(v >= t for v in vals)] += 1
                if any(min(vals[i] for i in g) >= t for g in gidx):
                    group_pass[t][c] += 1

    classes = sorted(total, key=lambda c: -total[c])
    rows = []
    for t in TPMS:
        for n in MIN_RUNS:
            kept = {c: sum(cnt for k, cnt in runs_hist[t][c].items() if k >= n) for c in classes}
            rows.append(dict(criterion=f"runs>={n}", tpm=t, kept_total=sum(kept.values()), **kept))
        kept = {c: group_pass[t][c] for c in classes}
        rows.append(dict(criterion="all_runs_of_a_group", tpm=t, kept_total=sum(kept.values()), **kept))

    cols = ["criterion", "tpm", "kept_total"] + classes
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# non-reference transcripts after the structural filter: {sum(total.values())} "
                 f"({', '.join(f'{c} {total[c]}' for c in classes)}); comparison groups: {len(groups)}\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"non-reference transcripts: {sum(total.values())}  " + "  ".join(f"{c}={total[c]}" for c in classes))
    print(f"comparison groups: {len(groups)}")
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))


if __name__ == "__main__":
    main()
