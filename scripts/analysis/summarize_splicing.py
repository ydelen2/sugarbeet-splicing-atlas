"""Summarize differential splicing results from rMATS and SUPPA2 for every comparison.

rMATS (JCEC tables)
  An event is "covered" when every replicate of both groups has at least --min-reads reads
  supporting the event (inclusion + skipping counts). Significant: FDR < --fdr and
  |IncLevelDifference| >= --dpsi, counted among covered events only.
  rMATS reports IncLevelDifference as group 1 minus group 2 (here reference minus test); the
  summary reports the direction as test minus reference so it can be compared with SUPPA2.

SUPPA2 (events.dpsi)
  Significant: p < --pval and |dPSI| >= --dpsi. The direction of dPSI is checked against the group
  PSI files (mean test minus mean reference) and the share of agreeing events is printed.

Usage
  python summarize_splicing.py --comparisons comparisons_pilot.tsv --rmats-dir results/rmats \
         --suppa-dir results/suppa --out results/summary/splicing_summary.tsv
"""
import argparse
import csv
import math
import os

RMATS_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (TypeError, ValueError):
        return None


def rmats_summary(path, min_reads, fdr_cut, dpsi_cut):
    total = covered = sig = up = down = 0
    genes = set()
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            total += 1
            counts = []
            for grp in ("1", "2"):
                inc = r[f"IJC_SAMPLE_{grp}"].split(",")
                skp = r[f"SJC_SAMPLE_{grp}"].split(",")
                counts += [(fnum(a) or 0) + (fnum(b) or 0) for a, b in zip(inc, skp)]
            if not counts or min(counts) < min_reads:
                continue
            covered += 1
            fdr, d = fnum(r["FDR"]), fnum(r["IncLevelDifference"])
            if fdr is None or d is None or fdr >= fdr_cut or abs(d) < dpsi_cut:
                continue
            sig += 1
            genes.add(r["GeneID"].strip('"'))
            if -d > 0:          # test minus reference
                up += 1
            else:
                down += 1
    return dict(total=total, covered=covered, significant=sig, higher_in_test=up, lower_in_test=down,
                genes=genes)


def read_psi_means(path):
    means = {}
    with open(path, encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals = [v for v in map(fnum, f[1:]) if v is not None]
            if vals:
                means[f[0]] = sum(vals) / len(vals)
    return means


def suppa_summary(suppa_dir, cid, ref, test, p_cut, dpsi_cut):
    dpsi_file = os.path.join(suppa_dir, "diff", cid, "events.dpsi")
    ref_mean = read_psi_means(os.path.join(suppa_dir, "groups", "events", cid, f"{ref}.psi"))
    test_mean = read_psi_means(os.path.join(suppa_dir, "groups", "events", cid, f"{test}.psi"))
    out = {}
    agree = checked = 0
    with open(dpsi_file, encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            ev, d, p = line.rstrip("\n").split("\t")[:3]
            etype = ev.split(";")[1].split(":")[0]
            s = out.setdefault(etype, dict(total=0, covered=0, significant=0, higher_in_test=0,
                                           lower_in_test=0, genes=set()))
            s["total"] += 1
            d, p = fnum(d), fnum(p)
            if d is None or p is None:
                continue
            s["covered"] += 1          # testable: dPSI and p-value available
            if p >= p_cut or abs(d) < dpsi_cut:
                continue
            s["significant"] += 1
            s["genes"].add(ev.split(";")[0])
            if ev in ref_mean and ev in test_mean and abs(test_mean[ev] - ref_mean[ev]) > 1e-6:
                checked += 1
                agree += (test_mean[ev] - ref_mean[ev] > 0) == (d > 0)
            if d > 0:
                s["higher_in_test"] += 1
            else:
                s["lower_in_test"] += 1
    direction = f"{agree}/{checked} significant events have dPSI with the sign of mean(test) minus mean(reference)"
    return out, direction


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--rmats-dir", required=True)
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-reads", type=int, default=20)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--dpsi", type=float, default=0.1)
    a = ap.parse_args()

    rows, overlap = [], []
    with open(a.comparisons, encoding="utf-8") as fh:
        comps = list(csv.DictReader(fh, delimiter="\t"))
    for c in comps:
        cid = c["comparison_id"]
        rmats_genes, suppa_genes = set(), set()
        for et in RMATS_TYPES:
            p = os.path.join(a.rmats_dir, cid, f"{et}.MATS.JCEC.txt")
            if os.path.isfile(p):
                s = rmats_summary(p, a.min_reads, a.fdr, a.dpsi)
                rmats_genes |= s["genes"]
                rows.append(dict(comparison=cid, tool="rMATS", event_type=et, **s))
        if os.path.isfile(os.path.join(a.suppa_dir, "diff", cid, "events.dpsi")):
            res, direction = suppa_summary(a.suppa_dir, cid, c["reference"], c["test"], a.pval, a.dpsi)
            print(f"{cid}: SUPPA2 direction check: {direction}")
            for et in sorted(res):
                suppa_genes |= res[et]["genes"]
                rows.append(dict(comparison=cid, tool="SUPPA2", event_type=et, **res[et]))
        # SUPPA2 joins the ids of overlapping genes with "_and_"; split them for the comparison
        suppa_split = set(g for x in suppa_genes for g in x.split("_and_"))
        both = rmats_genes & suppa_split
        overlap.append(dict(comparison=cid, rmats_genes=len(rmats_genes), suppa_genes=len(suppa_split),
                            both=len(both),
                            pct_of_suppa=round(100 * len(both) / len(suppa_split), 1) if suppa_split else 0,
                            pct_of_rmats=round(100 * len(both) / len(rmats_genes), 1) if rmats_genes else 0))
    for r in rows:
        r["genes"] = len(r["genes"])

    cols = ["comparison", "tool", "event_type", "total", "covered", "significant", "higher_in_test",
            "lower_in_test", "genes"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# rMATS: covered = every replicate >= {a.min_reads} reads, significant = FDR < {a.fdr} and "
                 f"|dPSI| >= {a.dpsi}. SUPPA2: covered = testable, significant = p < {a.pval} and |dPSI| >= {a.dpsi}\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))

    ocols = ["comparison", "rmats_genes", "suppa_genes", "both", "pct_of_suppa", "pct_of_rmats"]
    opath = os.path.splitext(a.out)[0] + "_gene_overlap.tsv"
    with open(opath, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# genes with at least one significant event of any type, per tool, and their overlap\n")
        w = csv.DictWriter(fh, fieldnames=ocols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(overlap)
    print("\n" + "\t".join(ocols))
    for r in overlap:
        print("\t".join(str(r[c]) for c in ocols))


if __name__ == "__main__":
    main()
