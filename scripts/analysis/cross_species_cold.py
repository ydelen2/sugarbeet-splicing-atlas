"""Cold-responsive differential splicing in sugar beet against the Arabidopsis cold DAS gene list of
Calixto et al. 2018 (Plant Cell 30:1424, doi 10.1105/tpc.18.00177).

Sugar beet genes are linked to Arabidopsis genes through the DIAMOND search of 10_domains.sbatch
(beet_vs_ath.tsv, best hit by bit score per protein, one line per beet gene through its best-scoring
protein) and the UniProt id mapping of the Arabidopsis reference proteome (accession to AGI locus code,
any AGI code listed for the accession). The Arabidopsis list is read from any text file: every AGI code
(AT[1-5CM]G followed by five digits) found in the file is taken, so the supplemental table can be given as
exported text or CSV without editing; make sure the file holds the DAS gene list and nothing else.

Sets compared (beet genes with an Arabidopsis link)
  background   genes with a testable SUPPA2 event in the comparisons whose identifier starts with
               --background-prefix (default cold, so that the background is what could be detected in the cold);
               with --ath-expressed, only genes whose Arabidopsis counterpart was expressed in the Arabidopsis
               experiment (and so could have been called there) are kept in every set
  cold DAS     genes with a high-confidence event in at least one cold comparison
  cold core    genes with a high-confidence event in all --n-cold cold comparisons
For each set the share of genes whose Arabidopsis counterpart is in the list is reported with the
hypergeometric upper-tail p-value against the background, and the number of distinct Arabidopsis counterparts
is reported as well because several beet genes can share one. The overlapping genes are written with their
event types and directions.

Usage
  python cross_species_cold.py --events results_full/summary/events_high_confidence_annotated.tsv \
         --suppa-dir results_full/suppa --comparisons metadata/comparisons_all.tsv \
         --map results_full/filtered/transcript_map.tsv \
         --beet-ath results_full/domains/beet_vs_ath.tsv --idmapping ref/proteomes/UP000006548_3702.idmapping.gz \
         --ath-list ref/calixto2018_cold_das_genes.txt --ath-expressed ref/calixto2018_expressed_genes.txt \
         --out-prefix results_full/summary/cross_species_cold
"""
import argparse
import collections
import csv
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splicing_regulators import read_diamond, testable_genes  # noqa: E402
from stats_util import hypergeom_upper  # noqa: E402

AGI = re.compile(r"\bAT[1-5CM]G\d{5}\b", re.I)


def read_idmapping(path):
    """{uniprot accession: {AGI codes}} from the reference proteome idmapping file."""
    agi_of = collections.defaultdict(set)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            for m in AGI.findall(f[2]):
                agi_of[f[0]].add(m.upper())
    return agi_of


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True, help="events_high_confidence_annotated.tsv")
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--comparisons", help="comparison table that defines the testable-gene background")
    ap.add_argument("--beet-ath", required=True)
    ap.add_argument("--idmapping", required=True)
    ap.add_argument("--ath-list", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-cold", type=int, default=6, help="number of cold comparisons that define the core")
    ap.add_argument("--background-prefix", default="cold", help="comparison identifiers that define the background")
    ap.add_argument("--ath-expressed", help="file with the AGI codes expressed in the Arabidopsis experiment")
    a = ap.parse_args()

    with open(a.ath_list, encoding="utf-8", errors="replace") as fh:
        ath_list = {m.upper() for m in AGI.findall(fh.read())}
    print(f"Arabidopsis list: {len(ath_list)} distinct AGI codes")
    agi_of = read_idmapping(a.idmapping)
    best, _ = read_diamond(a.beet_ath)
    link = {}                                        # beet gene -> (AGI codes, uniprot entry, bitscore)
    for prot, b in best.items():
        gene = prot.rsplit(".", 1)[0]
        acc = b["subject"].split("|")[1] if b["subject"].count("|") == 2 else b["subject"]
        codes = agi_of.get(acc, set())
        if codes and (gene not in link or b["bitscore"] > link[gene][2]):
            link[gene] = (codes, b["subject"], b["bitscore"])
    print(f"beet genes with a best Arabidopsis hit that has an AGI code: {len(link)}")

    comp_file = a.comparisons
    if a.comparisons and a.background_prefix:
        with open(a.comparisons, encoding="utf-8") as fh:
            lines = fh.readlines()
        keep = [lines[0]] + [l for l in lines[1:] if l.split("\t", 1)[0].startswith(a.background_prefix)]
        comp_file = a.out_prefix + "_background_comparisons.tsv"
        os.makedirs(os.path.dirname(os.path.abspath(comp_file)), exist_ok=True)
        with open(comp_file, "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(keep)
        print(f"background comparisons: {len(keep) - 1} with prefix {a.background_prefix}")
    if a.ath_expressed:
        with open(a.ath_expressed, encoding="utf-8", errors="replace") as fh:
            expressed = {m.upper() for m in AGI.findall(fh.read())}
        link = {g: v for g, v in link.items() if v[0] & expressed}
        print(f"Arabidopsis expressed genes: {len(expressed)}; beet genes whose counterpart is expressed: {len(link)}")
    bg = testable_genes(a.suppa_dir, a.map, comp_file) & set(link)
    cold_das, core = set(), set()
    ev_of = collections.defaultdict(list)
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cold_comps = [c for c in r["comparisons"].split(",") if c.startswith("cold")]
            if not cold_comps:
                continue
            for g in r["ref_gene"].split(","):
                if g:
                    cold_das.add(g)
                    ev_of[g].append(f"{r['event_type']}:{r['direction']}:{len(cold_comps)}")
                    if len(cold_comps) >= a.n_cold:
                        core.add(g)
    in_list = lambda g: bool(link[g][0] & ath_list)
    bg_hit = {g for g in bg if in_list(g)}
    rows = []
    for name, genes in (("background", bg), ("cold_das", cold_das & bg), ("cold_core", core & bg)):
        k = sum(1 for g in genes if in_list(g))
        expected = len(genes) * len(bg_hit) / len(bg) if bg else 0.0
        agi_in = {c for g in genes for c in link[g][0] & ath_list}
        rows.append(dict(set=name, genes=len(genes), counterpart_in_list=k, distinct_counterparts_in_list=len(agi_in),
                         share=round(k / len(genes), 4) if genes else "", expected=round(expected, 2),
                         enrichment=("" if expected == 0 or name == "background" else round(k / expected, 2)),
                         p_hypergeometric=("" if name == "background" else f"{hypergeom_upper(k, len(genes), len(bg_hit), len(bg)):.3g}")))
    with open(a.out_prefix + "_summary.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    genes_out = []
    for g in sorted(cold_das & bg):
        codes, subj, score = link[g]
        genes_out.append(dict(gene=g, ath_best_hit=subj, agi=";".join(sorted(codes)), bitscore=score,
                              in_arabidopsis_list=("yes" if in_list(g) else "no"),
                              cold_core=("yes" if g in core else "no"), events=";".join(ev_of[g])))
    with open(a.out_prefix + "_genes.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(genes_out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(genes_out)
    for r in rows:
        print(f"  {r['set']}: {r['counterpart_in_list']} of {r['genes']} genes have a counterpart in the list "
              f"(expected {r['expected']}, p {r['p_hypergeometric'] or 'reference'})")
    print(f"tables: {a.out_prefix}_summary.tsv  {a.out_prefix}_genes.tsv")


if __name__ == "__main__":
    main()
