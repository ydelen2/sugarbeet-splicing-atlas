"""Filter the merged StringTie transcriptome in two steps.

structural   Uses gffcompare class codes (from gffcmp.annotated.gtf) and exon counts.
             Kept:    =  k  (reference transcripts and their extended versions, any exon number)
                      j  m  n  o  (multi-exon only: novel isoforms of annotated genes)
                      u  (multi-exon only: novel intergenic loci)
             Dropped: single-exon novel transcripts (pre-mRNA or genomic DNA signal cannot be excluded),
                      x (antisense; unreliable because part of the data is unstranded),
                      i, y, p, s, c, e, r and transcripts without a class code.
expression   Keeps every class "=" transcript. Any other transcript must reach --min-tpm
               --criterion runs  (default) in at least --min-runs runs of the Salmon TPM matrix, or
               --criterion group in every run of at least one reference or test group of a comparison
                                 (needs --sheet and --comparisons).
             The TPM matrix is in SUPPA2 format (header = run names only).

Both steps write a GTF and a report. "structural" also writes transcript_map.tsv
(transcript, gene, class code, reference gene and transcript, exon number).

Usage
  python filter_transcripts.py structural --merged merged.gtf --annotated gffcmp.annotated.gtf --out-dir DIR
  python filter_transcripts.py expression --gtf DIR/structural.gtf --map DIR/transcript_map.tsv \
         --tpm tpm_all.tsv --min-tpm 1 --min-runs 2 --out DIR/final.gtf
"""
import argparse
import collections
import os
import re

KEEP_ANY = {"=", "k"}
KEEP_MULTI = {"j", "m", "n", "o", "u"}
ATTR = lambda name: re.compile(name + r' "([^"]*)"')
RE_T, RE_G, RE_C = ATTR("transcript_id"), ATTR("gene_id"), ATTR("class_code")
RE_RG, RE_RT = ATTR("ref_gene_id"), ATTR("cmp_ref")


def get(rx, s):
    m = rx.search(s)
    return m.group(1) if m else ""


def cmd_structural(a):
    info = {}                      # transcript -> dict
    exons = collections.Counter()
    with open(a.annotated, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = get(RE_T, f[8])
            if f[2] == "transcript":
                info[t] = dict(gene=get(RE_G, f[8]), code=get(RE_C, f[8]) or "none", chrom=f[0],
                               ref_gene=get(RE_RG, f[8]), ref_tx=get(RE_RT, f[8]))
            elif f[2] == "exon":
                exons[t] += 1

    keep = set()
    report = collections.defaultdict(lambda: [0, 0])   # code -> [kept, dropped]
    n_organelle = 0
    for t, d in info.items():
        multi = exons[t] > 1
        ok = d["code"] in KEEP_ANY or (d["code"] in KEEP_MULTI and multi)
        if a.drop_seq_prefix and d["chrom"].startswith(a.drop_seq_prefix):
            ok = False
            n_organelle += 1
        report[d["code"]][0 if ok else 1] += 1
        if ok:
            keep.add(t)

    os.makedirs(a.out_dir, exist_ok=True)
    n_lines = 0
    with open(a.merged, encoding="utf-8") as fin, \
         open(os.path.join(a.out_dir, "structural.gtf"), "w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) >= 9 and get(RE_T, f[8]) in keep:
                fout.write(line)
                n_lines += 1

    with open(os.path.join(a.out_dir, "transcript_map.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("transcript\tgene\tclass_code\tref_gene\tref_transcript\tn_exons\tkept_structural\n")
        for t in sorted(info):
            d = info[t]
            fh.write(f"{t}\t{d['gene']}\t{d['code']}\t{d['ref_gene']}\t{d['ref_tx']}\t{exons[t]}\t"
                     f"{'yes' if t in keep else 'no'}\n")

    with open(os.path.join(a.out_dir, "structural_report.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("class_code\tkept\tdropped\n")
        for code in sorted(report, key=lambda c: -sum(report[c])):
            fh.write(f"{code}\t{report[code][0]}\t{report[code][1]}\n")
            print(f"{code}\tkept {report[code][0]}\tdropped {report[code][1]}")
    with_ref = sum(1 for t in keep if info[t]["ref_gene"])
    print(f"transcripts in: {len(info)}  kept: {len(keep)}  GTF lines written: {n_lines}")
    print(f"kept transcripts with a reference gene id: {with_ref}")
    if a.drop_seq_prefix:
        print(f"transcripts dropped because their sequence name starts with '{a.drop_seq_prefix}': {n_organelle}")


def cmd_expression(a):
    code = {}
    with open(a.map, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            code[f[0]] = f[2]

    passed = set()
    with open(a.tpm, encoding="utf-8") as fh:
        samples = fh.readline().rstrip("\n").split("\t")
        n_runs = len(samples)
        if a.criterion == "group":
            # every run of at least one reference or test group of the comparison table
            import sys
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from prepare_splicing_inputs import comparison_groups, read_tsv
            if not (a.sheet and a.comparisons):
                raise SystemExit("--criterion group needs --sheet and --comparisons")
            col = {s: i for i, s in enumerate(samples)}
            sheet = read_tsv(a.sheet)
            gidx = []
            for comp in read_tsv(a.comparisons):
                for g in comparison_groups(sheet, comp):
                    missing = [r for r in g if r not in col]
                    if missing:
                        raise SystemExit(f"runs of comparison {comp['comparison_id']} missing from {a.tpm}: {missing[:5]}")
                    gidx.append([col[r] for r in g])
            rule = f"TPM >= {a.min_tpm} in every run of at least one of {len(gidx)} comparison groups"
        else:
            rule = f"TPM >= {a.min_tpm} in >= {a.min_runs} of {n_runs} runs"
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals = f[1:]
            if len(vals) != n_runs:
                raise SystemExit(f"row {f[0]}: {len(vals)} values, expected {n_runs}")
            if a.criterion == "group":
                v = [float(x) for x in vals]
                if any(min(v[i] for i in g) >= a.min_tpm for g in gidx):
                    passed.add(f[0])
            elif sum(float(x) >= a.min_tpm for x in vals) >= a.min_runs:
                passed.add(f[0])

    kept, stats = set(), collections.defaultdict(lambda: [0, 0])
    seen = set()
    with open(a.gtf, encoding="utf-8") as fin:
        for line in fin:
            f = line.split("\t")
            if len(f) >= 9:
                seen.add(get(RE_T, f[8]))
    for t in seen:
        c = code.get(t, "none")
        ok = c == "=" or t in passed
        stats[c][0 if ok else 1] += 1
        if ok:
            kept.add(t)

    with open(a.gtf, encoding="utf-8") as fin, open(a.out, "w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            f = line.split("\t")
            if len(f) >= 9 and get(RE_T, f[8]) in kept:
                fout.write(line)

    rep = os.path.splitext(a.out)[0] + "_report.tsv"
    with open(rep, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# class '=' always kept; others: {rule}\n")
        fh.write("class_code\tkept\tdropped\n")
        for c in sorted(stats, key=lambda c: -sum(stats[c])):
            fh.write(f"{c}\t{stats[c][0]}\t{stats[c][1]}\n")
            print(f"{c}\tkept {stats[c][0]}\tdropped {stats[c][1]}")
    print(f"transcripts in: {len(seen)}  kept: {len(kept)}  (class '=' always kept; others: {rule})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("structural")
    p.add_argument("--merged", required=True); p.add_argument("--annotated", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--drop-seq-prefix", default="", help="drop transcripts on sequences whose name starts with this")
    p.set_defaults(func=cmd_structural)
    p = sub.add_parser("expression")
    p.add_argument("--gtf", required=True); p.add_argument("--map", required=True); p.add_argument("--tpm", required=True)
    p.add_argument("--min-tpm", type=float, default=1.0); p.add_argument("--min-runs", type=int, default=2)
    p.add_argument("--criterion", choices=["runs", "group"], default="runs")
    p.add_argument("--sheet", help="sample sheet (criterion group)")
    p.add_argument("--comparisons", help="comparison table (criterion group)")
    p.add_argument("--out", required=True); p.set_defaults(func=cmd_expression)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
