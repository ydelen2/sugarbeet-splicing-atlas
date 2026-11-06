"""Pfam domains affected by the coding high-confidence events.

For every event inside a CDS the annotated mRNA that carries one form of the event is identified
(event_forms.py) and the alternative region is placed on its protein:
  region in the annotated form   the region is exonic in the mRNA, so its bases map to a protein interval
                                 [aa_start, aa_end]; the other form removes these residues
  region absent                  the mRNA lacks the region, so the other form inserts residues at the
                                 protein position of the junction (aa_start = aa_end = that position)
Pfam domains of the same protein (hmmsearch domtblout of 10_domains.sbatch, envelope coordinates) that
overlap the interval or contain the insertion point are listed with the number of their residues inside the
region. For MX events the region of the form carried by the mRNA is used and the other exon is the
replacement.

Only frame-preserving events (region length divisible by 3, MX: length difference divisible by 3) yield an
in-frame insertion or deletion; frame-shifting events are kept in the table with their frame so that the
domain at the shift point is known, but the summary counts the frame-preserving ones separately.

Usage
  python domain_impact.py --consequences results_full/summary/events_high_confidence_consequences.tsv \
         --gff ref/BeetSet-3.genes.2604.gff --domtbl results_full/domains/pfam.domtblout \
         --out results_full/summary/events_high_confidence_domains.tsv
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from event_forms import parse_event, read_models, which_form  # noqa: E402
from splicing_regulators import read_domtbl  # noqa: E402


def cds_offset(cds_list, strand, gpos):
    """0-based coding offset of a genomic position inside the annotated CDS, or None."""
    off = 0
    blocks = cds_list if strand == "+" else list(reversed(cds_list))
    for s, e, _ in blocks:
        if s <= gpos <= e:
            return off + (gpos - s if strand == "+" else e - gpos)
        off += e - s + 1
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--consequences", required=True, help="events_high_confidence_consequences.tsv")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    mrnas_of, info, exons, cds = read_models(a.gff)
    doms = read_domtbl(a.domtbl)
    shared = len(set(doms) & set(info))
    print(f"proteins with domains {len(doms)}, of them with an mRNA id in the GFF {shared}")
    if shared == 0:
        sys.exit("protein ids of the domain table do not match the mRNA ids of the GFF")

    out = []
    with open(a.consequences, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["cds_context"] == "CDS"]
    for r in rows:
        chrom, strand, etype, f1, f2, regions = parse_event(r["suppa_event"])
        hit = None
        for gene in [g for g in r["ref_gene"].split(",") if g]:
            for m in mrnas_of.get(gene, []):
                if not cds.get(m) or info[m][0] != chrom or info[m][1] != strand:
                    continue
                form = which_form(exons[m], etype, f1, f2, regions)
                if form is None:
                    continue
                region = regions[0] if etype != "MX" else regions[form - 1]
                other = None if etype != "MX" else regions[2 - form]
                s, e = region
                inside = any(es <= s and e <= ee for es, ee in exons[m])
                if inside:
                    if strand == "+":
                        o1, o2 = cds_offset(cds[m], strand, s), cds_offset(cds[m], strand, e)
                    else:
                        o1, o2 = cds_offset(cds[m], strand, e), cds_offset(cds[m], strand, s)
                    if o1 is None or o2 is None:
                        continue
                    aa1, aa2 = o1 // 3 + 1, o2 // 3 + 1
                    placement = "region in annotated form"
                else:
                    # insertion point: the exon end just before the region in transcription order
                    if strand == "+":
                        anchor = max((ee for es, ee in exons[m] if ee < s), default=None)
                    else:
                        anchor = min((es for es, ee in exons[m] if es > e), default=None)
                    o = None if anchor is None else cds_offset(cds[m], strand, anchor)
                    if o is None:
                        continue
                    aa1 = aa2 = o // 3 + 1
                    placement = "region absent, insertion point"
                hit = (gene, m, form, aa1, aa2, placement, region, other)
                break
            if hit:
                break
        if not hit:
            out.append(dict(suppa_event=r["suppa_event"], event_type=etype, gene=r["ref_gene"], frame=r["frame"], mrna="",
                            annotated_form="", placement="unresolved", aa_start="", aa_end="", region_aa="",
                            protein_length="", domains_affected="", n_domains=0, n_comparisons=r["n_comparisons"],
                            groups=r["groups"], direction=r["direction"], description=r["description"]))
            continue
        gene, m, form, aa1, aa2, placement, region, other = hit
        region_aa = (region[1] - region[0] + 1) / 3
        plen = sum(e - s + 1 for s, e, _ in cds[m]) // 3
        affected = []
        for name, acc, d_from, d_to, evalue in doms.get(m, []):
            if d_from <= aa2 and d_to >= aa1:
                inside = min(aa2, d_to) - max(aa1, d_from) + 1
                affected.append(f"{name}({d_from}-{d_to}):{inside}")
        out.append(dict(suppa_event=r["suppa_event"], event_type=etype, gene=gene, frame=r["frame"], mrna=m,
                        annotated_form=form, placement=placement, aa_start=aa1, aa_end=aa2,
                        region_aa=f"{region_aa:.1f}", protein_length=plen, domains_affected=";".join(affected),
                        n_domains=len(affected), n_comparisons=r["n_comparisons"], groups=r["groups"],
                        direction=r["direction"], description=r["description"]))

    cols = list(out[0].keys())
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    print(f"CDS events {len(out)}, placed on a protein {sum(1 for d in out if d['mrna'])}")
    for label, keep in (("frame-preserving", lambda d: d["frame"] == "0"), ("frame-shifting", lambda d: d["frame"] != "0")):
        rows_k = [d for d in out if keep(d) and d["mrna"]]
        with_dom = [d for d in rows_k if d["n_domains"]]
        by_type = collections.Counter(d["event_type"] for d in with_dom)
        print(f"  {label}: {len(rows_k)} placed events, {len(with_dom)} overlap a Pfam domain ("
              + ", ".join(f"{t} {by_type[t]}" for t in ("SE", "A3", "A5", "RI", "MX") if by_type[t]) + ")")
    # one count per event and domain name (a region can overlap several copies of a repeated domain)
    dom_names = collections.Counter(n for d in out if d["frame"] == "0"
                                    for n in {x.split("(")[0] for x in d["domains_affected"].split(";") if x})
    print("  events per domain, frame-preserving events (stop codons not considered here, see nmd_features.py): "
          + ", ".join(f"{k} {v}" for k, v in dom_names.most_common(10)))
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
