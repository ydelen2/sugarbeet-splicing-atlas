"""SR protein candidates of one proteome, with the same rule as splicing_regulators.py.

A protein is an SR candidate when it has one or two RRM_1 (PF00076) domains and a downstream RS domain
by the plant definition of Barta et al. (2010): a 50-residue window with at least 20 percent RS or SR
dipeptides (five dipeptides in a 50-residue window). The number of zinc knuckles (zf-CCHC) is reported for the RSZ and RS2Z
subfamilies. Its subfamily label is the best
DIAMOND hit in the Arabidopsis UniProt proteome when that hit is one of the SR entries of the family
table; otherwise the protein is kept as "SR-related" (rule only). One protein per gene is kept, the
longest; the protein to gene map comes from the NCBI GFF3 (CDS features carry protein_id and gene or locus_tag) or,
for BeetSet-3, from the protein id itself.

Usage
  python sr_candidates.py --species quinoa --proteins protein.faa --domtbl rrm.domtblout \
      --diamond vs_ath.tsv --families splicing_factor_families.tsv [--gff genomic.gff] --out quinoa_sr.tsv
The output table has: species, protein, gene, length, n_rrm, rs_dipeptides_max, rs_fraction_max, n_zinc_knuckles,
best_ath, subfamily, ath_gene, pident, and the sequence, so the trees can be built from the tables alone.
"""
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splicing_regulators import read_fasta, read_domtbl, read_diamond, entry_name, rs_domain, is_sr, RRM, ZNF  # noqa: E402


def gene_map_from_gff(path):
    """{protein_id: gene} from NCBI GFF3 CDS features."""
    m = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            # RefSeq CDS features carry gene=, GenBank submissions carry locus_tag=; the mRNA Parent is the last resort
            pid, gene = attrs.get("protein_id"), attrs.get("gene") or attrs.get("locus_tag") or attrs.get("Parent", "")
            if pid and pid not in m:
                m[pid] = gene
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--diamond", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    seqs = read_fasta(a.proteins)
    doms = read_domtbl(a.domtbl)
    best, _ = read_diamond(a.diamond)
    sr_entries = {}
    with open(a.families, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"] in ("SR", "SR-like"):
                sr_entries[r["uniprot_entry"]] = r
    gene_of = gene_map_from_gff(a.gff) if a.gff else {}

    rows = []
    for prot, seq in seqs.items():
        dlist = doms.get(prot, [])
        rrms = [d for d in dlist if d[1] == RRM]
        if not rrms or len(rrms) > 2:
            continue
        dip, wstart, frac = rs_domain(seq, rrms[-1][3])
        if not is_sr(rrms, dip):
            continue
        n_znf = sum(1 for d in dlist if d[1] == ZNF)
        b = best.get(prot)
        ent = entry_name(b["subject"]) if b else ""
        fam = sr_entries.get(ent)
        gene = gene_of.get(prot) or (prot.rsplit(".", 1)[0] if re.match(r"^Bv", prot) else prot)
        rows.append(dict(species=a.species, protein=prot, gene=gene, length=len(seq), n_rrm=len(rrms),
                         rs_dipeptides_max=dip, rs_fraction_max=round(frac, 3), n_zinc_knuckles=n_znf, best_ath=ent,
                         subfamily=fam["subfamily"] if fam else "SR-related", ath_gene=fam["gene"] if fam else "",
                         pident=b["pident"] if b else "", sequence=seq))
    # longest protein per gene
    by_gene = {}
    for r in rows:
        if r["gene"] not in by_gene or r["length"] > by_gene[r["gene"]]["length"]:
            by_gene[r["gene"]] = r
    out = sorted(by_gene.values(), key=lambda r: (r["subfamily"], r["gene"]))
    cols = list(out[0].keys()) if out else ["species", "protein", "gene"]
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n_sr = sum(1 for r in out if r["subfamily"] != "SR-related")
    print(f"{a.species}: proteins with 1-2 RRMs and an RS domain: {len(rows)}; genes: {len(out)}; "
          f"with an Arabidopsis SR best hit: {n_sr}; SR-related (rule only): {len(out) - n_sr}")


if __name__ == "__main__":
    main()
