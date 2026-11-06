"""Expression of the nonsense-mediated decay factors across the comparisons.

If the premature-stop forms accumulate under a condition because the decay pathway is weakened rather than
because splicing shifts, the core factors of the pathway should be down-regulated in that condition. This script
pulls the DESeq2 results (differential_expression.py, <out-prefix>_results/<comparison>.tsv) of the sugar beet genes
whose best Arabidopsis hit is one of the core NMD factors UPF1 (AT5G47010), UPF2 (AT2G39260), UPF3 (AT1G33980)
and SMG7 (AT5G19400), identified through the UniProt id mapping of the Arabidopsis reference proteome.

Outputs (<out-prefix>_...)
  genes.tsv     the sugar beet genes assigned to each factor with their best hit, identity and coverage
  table.tsv     log2 fold change (test over reference), adjusted p and the DE call of every gene in every comparison
  summary.txt   per condition group, how many of the factor genes are significantly down- or up-regulated

Usage
  python nmd_factor_expression.py --beet-ath results_full/domains/beet_vs_ath.tsv \
         --idmapping ref/proteomes/UP000006548_3702.idmapping.gz --de-dir results_full/expression/de_results \
         --comparisons metadata/comparisons_all.tsv --out-prefix results_full/summary/nmd_factors
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_species_cold import read_idmapping  # noqa: E402
from splicing_regulators import read_diamond  # noqa: E402

FACTORS = {"AT5G47010": "UPF1", "AT2G39260": "UPF2", "AT1G33980": "UPF3", "AT5G19400": "SMG7"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beet-ath", required=True)
    ap.add_argument("--idmapping", required=True)
    ap.add_argument("--de-dir", required=True, help="directory with one DESeq2 result table per comparison")
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    a = ap.parse_args()

    agi_of = read_idmapping(a.idmapping)
    best, _ = read_diamond(a.beet_ath)
    genes = {}                                         # beet gene -> (factor, hit row)
    for prot, b in best.items():
        acc = b["subject"].split("|")[1] if b["subject"].count("|") == 2 else b["subject"]
        codes = agi_of.get(acc, set()) & set(FACTORS)
        if not codes:
            continue
        gene = prot.rsplit(".", 1)[0]
        if gene not in genes or b["bitscore"] > genes[gene][1]["bitscore"]:
            genes[gene] = (FACTORS[sorted(codes)[0]], b)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    with open(a.out_prefix + "_genes.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("gene\tfactor\tbest_ath\tpident\tqcov\tscov\tbitscore\n")
        for g, (fac, b) in sorted(genes.items(), key=lambda x: (x[1][0], x[0])):
            fh.write(f"{g}\t{fac}\t{b['subject']}\t{b['pident']}\t{b['qcov']}\t{b['scov']}\t{b['bitscore']}\n")
    print(f"factor genes: {len(genes)} ({collections.Counter(f for f, _ in genes.values())})")

    with open(a.comparisons, encoding="utf-8") as fh:
        comps = [r["comparison_id"] for r in csv.DictReader(fh, delimiter="\t")]
    rows, calls = [], collections.defaultdict(collections.Counter)
    for cid in comps:
        path = os.path.join(a.de_dir, cid + ".tsv")
        if not os.path.exists(path):
            print(f"missing DE table for {cid}", file=sys.stderr)
            continue
        with open(path, encoding="utf-8") as fh:
            res = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
        group = cid.split("_", 1)[0]
        for g, (fac, b) in sorted(genes.items()):
            r = res.get(g)
            lfc = float(r["log2FoldChange"]) if r and r["log2FoldChange"] not in ("", "nan") else float("nan")
            padj = float(r["padj"]) if r and r["padj"] not in ("", "nan", "NA") else float("nan")
            call = "down" if padj < a.fdr and lfc <= -a.lfc else "up" if padj < a.fdr and lfc >= a.lfc else "no"
            rows.append(dict(comparison=cid, group=group, gene=g, factor=fac,
                             base_mean=(r["baseMean"] if r else ""), log2fc=(f"{lfc:.3f}" if lfc == lfc else ""),
                             padj=(f"{padj:.3g}" if padj == padj else ""), de=call))
            calls[group][call] += 1
    with open(a.out_prefix + "_table.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(a.out_prefix + "_summary.txt", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# NMD factor genes ({len(genes)}) over {len(comps)} comparisons; DE at adjusted p < {a.fdr} and "
                 f"|log2 fold change| >= {a.lfc}\n")
        fh.write("group\tgene_comparison_pairs\tdown\tup\tnot_de\n")
        for group in sorted(calls):
            c = calls[group]
            fh.write(f"{group}\t{sum(c.values())}\t{c['down']}\t{c['up']}\t{c['no']}\n")
        fh.write("\n# per gene: comparisons in which the gene is down / up\n")
        for g, (fac, _) in sorted(genes.items(), key=lambda x: (x[1][0], x[0])):
            down = [r["comparison"] for r in rows if r["gene"] == g and r["de"] == "down"]
            up = [r["comparison"] for r in rows if r["gene"] == g and r["de"] == "up"]
            fh.write(f"{fac}\t{g}\tdown: {', '.join(down) or '-'}\tup: {', '.join(up) or '-'}\n")
    print(open(a.out_prefix + "_summary.txt", encoding="utf-8").read())


if __name__ == "__main__":
    main()
