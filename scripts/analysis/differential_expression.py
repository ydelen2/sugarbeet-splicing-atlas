"""Differential gene expression per comparison (PyDESeq2 on STAR gene counts) and its overlap with the
high-confidence differential splicing events.

Counts: <runres-dir>/<run>/<run>.ReadsPerGene.out.tab written by STAR with --quantMode GeneCounts against
the BeetSet-3 annotation. The count column follows the strandedness call of 04a_stringtie
(<runres-dir>/<run>/strandedness.tsv): unstranded column 2, forward column 3, reverse column 4.
Comparison groups are resolved exactly as for SUPPA2 and rMATS (comparison_groups of
prepare_splicing_inputs.py), so the same runs enter all three analyses.

Model: ~ group (reference, test) per comparison; genes with fewer than --min-count reads in total over the
runs of the comparison are removed before fitting. Differentially expressed (DE): BH-adjusted p < --fdr
and |log2 fold change| >= --lfc (test over reference, unshrunk estimates).

Overlap with differential splicing (DAS): the DAS genes of a comparison are the reference genes of its
high-confidence events. The population for the test is the set of genes that both DESeq2 tested (adjusted
p not NA) and SUPPA2 could test in that comparison (at least one event with a dPSI and a p-value, mapped
through transcript_map.tsv); DE and DAS gene sets are intersected with it, the overlap expected by chance
is DE x DAS / population and the p-value is the hypergeometric upper tail.

Outputs (<out-prefix>_...)
  de_summary.tsv        per comparison: runs, genes tested, DE genes (up, down), DAS genes, overlap and test
  de_das_genes.tsv      the DAS genes of every comparison with their expression result and event types
  de_regulators.tsv     expression result of the splicing-regulator family genes (--families) per comparison
  de_results/<cid>.tsv  the full DESeq2 result table of every comparison

Usage
  python differential_expression.py --sheet metadata/samples_all.tsv --comparisons metadata/comparisons_all.tsv \
         --runres-dir results_full/per_run --events results_full/summary/event_concordance_high_confidence.tsv \
         --suppa-dir results_full/suppa --map results_full/filtered/transcript_map.tsv \
         --families results_full/domains/splicing_regulators_families.tsv \
         --out-prefix results_full/expression/de --threads 8
"""
import argparse
import collections
import csv
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_splicing_inputs import comparison_groups, read_tsv  # noqa: E402
from stats_util import hypergeom_upper  # noqa: E402

COUNT_COLUMN = {"unstranded": 1, "forward": 2, "reverse": 3}


def read_counts(runres_dir, run):
    """(pandas Series of gene counts, strandedness) of one run."""
    with open(os.path.join(runres_dir, run, "strandedness.tsv"), encoding="utf-8") as fh:
        strand = fh.read().split("\t")[1].strip()
    col = COUNT_COLUMN[strand]
    counts = {}
    with open(os.path.join(runres_dir, run, f"{run}.ReadsPerGene.out.tab"), encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("N_"):
                continue
            counts[f[0]] = int(f[col])
    return pd.Series(counts, name=run), strand


def testable_ref_genes(suppa_dir, cid, ref_of):
    """Reference genes with at least one testable SUPPA2 event in one comparison."""
    genes = set()
    with open(os.path.join(suppa_dir, "diff", cid, "events.dpsi"), encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            ev, d, p = line.rstrip("\n").split("\t")[:3]
            if d != "nan" and p != "nan":
                for g in ev.split(";", 1)[0].split("_and_"):
                    genes.update(ref_of.get(g, ()))
    return genes


def run_deseq(counts, groups, threads):
    """DESeq2 result table (index gene) for a counts table (index run) and a run -> level mapping."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    meta = pd.DataFrame({"group": [groups[r] for r in counts.index]}, index=counts.index)
    inference = DefaultInference(n_cpus=threads)
    try:
        dds = DeseqDataSet(counts=counts, metadata=meta, design="~group", inference=inference, quiet=True)
    except TypeError:                                    # PyDESeq2 releases before 0.5 use design_factors
        dds = DeseqDataSet(counts=counts, metadata=meta, design_factors="group", inference=inference, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "test", "ref"], inference=inference, quiet=True)
    st.summary()
    return st.results_df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--runres-dir", required=True)
    ap.add_argument("--events", required=True, help="event_concordance_high_confidence.tsv")
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--map", required=True, help="transcript_map.tsv")
    ap.add_argument("--families", help="splicing_regulators_families.tsv")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--min-count", type=int, default=10)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    sheet = read_tsv(a.sheet)
    comps = read_tsv(a.comparisons)
    ref_of = collections.defaultdict(set)
    with open(a.map, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ref_gene"]:
                ref_of[r["gene"]].add(r["ref_gene"])

    das = collections.defaultdict(lambda: collections.defaultdict(list))   # cid -> gene -> [event types]
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for g in r["ref_genes"].split(","):
                if g:
                    das[r["comparison"]][g].append(r["event_type"])

    fam_of = {}
    if a.families:
        with open(a.families, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                fam_of[r["gene"]] = (r["family"], r["subfamily"], r["ath_gene"], r["best_ath"])

    res_dir = a.out_prefix + "_results"
    os.makedirs(res_dir, exist_ok=True)
    cache = {}
    summary, gene_rows, reg_rows = [], [], []
    for c in comps:
        cid = c["comparison_id"]
        ref, test = comparison_groups(sheet, c)
        strands = set()
        series = []
        for run in ref + test:
            if run not in cache:
                cache[run] = read_counts(a.runres_dir, run)
            s, strand = cache[run]
            series.append(s)
            strands.add(strand)
        counts = pd.concat(series, axis=1).T.fillna(0).astype(int)
        counts = counts.loc[:, counts.sum(axis=0) >= a.min_count]
        groups = {r: "ref" for r in ref}
        groups.update({r: "test" for r in test})
        res = run_deseq(counts, groups, a.threads)
        res.index.name = "gene"
        res.to_csv(os.path.join(res_dir, f"{cid}.tsv"), sep="\t", float_format="%.6g")

        tested = set(res.index[res["padj"].notna()])
        de_up = set(res.index[(res["padj"] < a.fdr) & (res["log2FoldChange"] >= a.lfc)])
        de_down = set(res.index[(res["padj"] < a.fdr) & (res["log2FoldChange"] <= -a.lfc)])
        de = de_up | de_down
        das_genes = set(das.get(cid, {}))
        pop = tested & testable_ref_genes(a.suppa_dir, cid, ref_of)
        de_p, das_p = de & pop, das_genes & pop
        both = de_p & das_p
        expected = len(de_p) * len(das_p) / len(pop) if pop else 0.0
        pval = hypergeom_upper(len(both), len(das_p), len(de_p), len(pop))
        summary.append(dict(
            comparison=cid, n_reference=len(ref), n_test=len(test), strandedness=";".join(sorted(strands)),
            genes_in_model=counts.shape[1], genes_tested=len(tested), de_genes=len(de), de_up=len(de_up),
            de_down=len(de_down), das_genes=len(das_genes), population=len(pop), de_in_population=len(de_p),
            das_in_population=len(das_p), de_and_das=len(both), de_only=len(de_p) - len(both),
            das_only=len(das_p) - len(both), expected_overlap=round(expected, 2),
            enrichment=("" if expected == 0 else round(len(both) / expected, 2)), p_hypergeometric=f"{pval:.3g}"))
        for g in sorted(das_genes):
            row = res.loc[g] if g in res.index else None
            gene_rows.append(dict(
                comparison=cid, gene=g, event_types=",".join(sorted(set(das[cid][g]))),
                base_mean=("" if row is None else f"{row['baseMean']:.2f}"),
                log2fc=("" if row is None else f"{row['log2FoldChange']:.3f}"),
                padj=("" if row is None or pd.isna(row["padj"]) else f"{row['padj']:.3g}"),
                de=("yes" if g in de else ("no" if g in tested else "not tested"))))
        for g, (fam, sub, ath, ent) in sorted(fam_of.items()):
            if g not in res.index:
                continue
            row = res.loc[g]
            reg_rows.append(dict(
                comparison=cid, gene=g, family=fam, subfamily=sub, ath_gene=ath, best_ath=ent,
                base_mean=f"{row['baseMean']:.2f}", log2fc=f"{row['log2FoldChange']:.3f}",
                padj=("" if pd.isna(row["padj"]) else f"{row['padj']:.3g}"),
                de=("yes" if g in de else ("no" if g in tested else "not tested")),
                das=("yes" if g in das_genes else "no")))
        print(f"{cid}: {len(ref)} vs {len(test)} runs, {counts.shape[1]} genes in the model, {len(tested)} tested, "
              f"DE {len(de)} (up {len(de_up)}, down {len(de_down)}); DAS {len(das_genes)}, population {len(pop)}, "
              f"DE and DAS {len(both)} (expected {expected:.1f}, p {pval:.3g})")

    def write(path, rows):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"table: {path}")

    with open(a.out_prefix + "_summary.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# DE: BH-adjusted p < {a.fdr} and |log2 fold change| >= {a.lfc} (test over reference); genes with "
                 f"fewer than {a.min_count} reads over the runs of a comparison were removed before the fit; "
                 "population = genes tested by DESeq2 and testable by SUPPA2 in the comparison; expected_overlap = "
                 "de_in_population x das_in_population / population; p_hypergeometric = upper tail\n")
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    print(f"table: {a.out_prefix}_summary.tsv")
    write(a.out_prefix + "_das_genes.tsv", gene_rows)
    if reg_rows:
        write(a.out_prefix + "_regulators.tsv", reg_rows)


if __name__ == "__main__":
    main()
