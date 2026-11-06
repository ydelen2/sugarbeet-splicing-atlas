"""Prepare input files for SUPPA2 and rMATS from Salmon and STAR outputs.

Subcommands
  tpm     Merge Salmon quant.sf files into one SUPPA2 expression matrix (transcript TPM).
  split   Split a SUPPA2 matrix (TPM or PSI) into one file per comparison group.
  rmats   Write rMATS BAM lists and run parameters (read length, library type) per comparison.

SUPPA2 matrix format: the first line holds the sample names only; every following line holds
an identifier followed by one value per sample, tab separated.

Comparisons are read from a table with the columns
  comparison_id, study, subset (column=value[;column=value] or -), group_column, reference, test

Examples
  python prepare_splicing_inputs.py tpm   --sheet samples_pilot.tsv --salmon-dir results/salmon --out tpm.tsv
  python prepare_splicing_inputs.py split --sheet samples_pilot.tsv --comparisons comparisons_pilot.tsv \
         --matrix events.psi --suffix psi --out-dir results/suppa/groups
  python prepare_splicing_inputs.py rmats --sheet samples_pilot.tsv --comparisons comparisons_pilot.tsv \
         --bam-dir bulk/bam --runres-dir results/per_run --out-dir results/rmats/inputs
"""
import argparse
import csv
import os
import sys


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def comparison_groups(sheet_rows, comp):
    """Return (reference_runs, test_runs) for one comparison row."""
    # only runs cleared for differential analysis take part in comparisons
    rows = [r for r in sheet_rows if r["study"] == comp["study"] and r.get("use_das", "yes") == "yes"]
    subset = comp["subset"].strip()
    if subset and subset != "-":
        for cond in subset.split(";"):
            col, val = cond.split("=", 1)
            rows = [r for r in rows if r[col.strip()] == val.strip()]
    col = comp["group_column"]
    ref = sorted(r["run"] for r in rows if r[col] == comp["reference"])
    test = sorted(r["run"] for r in rows if r[col] == comp["test"])
    if len(ref) < 2 or len(test) < 2:
        sys.exit(f"comparison {comp['comparison_id']}: needs at least 2 runs per group, got {len(ref)} and {len(test)}")
    return ref, test


def cmd_tpm(args):
    # runs marked use_das=excluded were never quantified
    runs = [r["run"] for r in read_tsv(args.sheet) if r.get("use_das") != "excluded"]
    tpm, missing = {}, []
    ids = None
    for run in runs:
        qf = os.path.join(args.salmon_dir, run, "quant.sf")
        if not os.path.isfile(qf):
            missing.append(run)
            continue
        vals = {}
        with open(qf, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                vals[r["Name"]] = r["TPM"]
        if ids is None:
            ids = list(vals)
        elif set(vals) != set(ids):
            sys.exit(f"{run}: transcript set differs from the first sample; all runs must use the same index")
        tpm[run] = vals
    if missing:
        sys.exit(f"missing quant.sf for {len(missing)} runs, first: {missing[:5]}")
    samples = list(tpm)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(samples) + "\n")
        for t in ids:
            fh.write(t + "\t" + "\t".join(tpm[s][t] for s in samples) + "\n")
    print(f"{args.out}: {len(ids)} transcripts x {len(samples)} samples")


def read_matrix(path):
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        body = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    # SUPPA writes the header without a leading field; tolerate one if present
    if body and len(body[0]) == len(header):
        header = header[1:]
    return header, body


def cmd_split(args):
    sheet = read_tsv(args.sheet)
    header, body = read_matrix(args.matrix)
    col_of = {s: i + 1 for i, s in enumerate(header)}
    for comp in read_tsv(args.comparisons):
        ref, test = comparison_groups(sheet, comp)
        outdir = os.path.join(args.out_dir, comp["comparison_id"])
        os.makedirs(outdir, exist_ok=True)
        for name, runs in ((comp["reference"], ref), (comp["test"], test)):
            absent = [r for r in runs if r not in col_of]
            if absent:
                sys.exit(f"{comp['comparison_id']}: runs not in {args.matrix}: {absent}")
            out = os.path.join(outdir, f"{name}.{args.suffix}")
            with open(out, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\t".join(runs) + "\n")
                for row in body:
                    fh.write(row[0] + "\t" + "\t".join(row[col_of[r]] for r in runs) + "\n")
        print(f"{comp['comparison_id']}: {comp['reference']} {len(ref)} runs, {comp['test']} {len(test)} runs")


LIBTYPE = {"unstranded": "fr-unstranded", "reverse": "fr-firststrand", "forward": "fr-secondstrand"}


def cmd_rmats(args):
    sheet = read_tsv(args.sheet)
    by_run = {r["run"]: r for r in sheet}
    comps = read_tsv(args.comparisons)
    if args.only:
        comps = [c for c in comps if c["comparison_id"] == args.only]
        if not comps:
            sys.exit(f"comparison {args.only} not found in {args.comparisons}")
    for comp in comps:
        ref, test = comparison_groups(sheet, comp)
        runs = ref + test
        strands = set()
        for run in runs:
            sf = os.path.join(args.runres_dir, run, "strandedness.tsv")
            if not os.path.isfile(sf):
                sys.exit(f"{run}: {sf} not found (run 04a_stringtie first)")
            with open(sf, encoding="utf-8") as fh:
                strands.add(fh.read().split("\t")[1])
        if len(strands) != 1:
            sys.exit(f"{comp['comparison_id']}: mixed library strandedness {strands}")
        bam = lambda r: os.path.join(args.bam_dir, r, f"{r}.Aligned.sortedByCoord.out.bam")
        outdir = os.path.join(args.out_dir, comp["comparison_id"])
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "b1.txt"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(",".join(bam(r) for r in ref) + "\n")
        with open(os.path.join(outdir, "b2.txt"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(",".join(bam(r) for r in test) + "\n")
        read_len = max(int(by_run[r]["read_len"]) for r in runs)
        with open(os.path.join(outdir, "params.txt"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"READ_LEN={read_len}\nLIBTYPE={LIBTYPE[strands.pop()]}\n"
                     f"B1_LABEL={comp['reference']}\nB2_LABEL={comp['test']}\n")
        print(f"{comp['comparison_id']}: b1={len(ref)} b2={len(test)} read_len={read_len}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("tpm"); p.add_argument("--sheet", required=True); p.add_argument("--salmon-dir", required=True)
    p.add_argument("--out", required=True); p.set_defaults(func=cmd_tpm)
    p = sub.add_parser("split"); p.add_argument("--sheet", required=True); p.add_argument("--comparisons", required=True)
    p.add_argument("--matrix", required=True); p.add_argument("--suffix", required=True, choices=["psi", "tpm"])
    p.add_argument("--out-dir", required=True); p.set_defaults(func=cmd_split)
    p = sub.add_parser("rmats"); p.add_argument("--sheet", required=True); p.add_argument("--comparisons", required=True)
    p.add_argument("--bam-dir", required=True); p.add_argument("--runres-dir", required=True)
    p.add_argument("--only", help="prepare this comparison only, so a problem in another one cannot stop the job")
    p.add_argument("--out-dir", required=True); p.set_defaults(func=cmd_rmats)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
