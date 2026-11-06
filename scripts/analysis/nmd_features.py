"""Premature termination codons and nonsense-mediated decay features of the high-confidence events.

For every event the two transcript forms are built on an annotated mRNA of the reference gene whose exon
chain contains the junctions of one form (event_forms.py). Both forms are translated from the annotated
start codon of that mRNA to the first in-frame stop codon, and for each form the script reports
  stop            genomic position of the first base of the stop codon (or none when no stop is reached)
  ptc             yes when the stop lies upstream of the annotated stop codon of the mRNA in transcription
                  order (a premature termination codon), no when it is the annotated stop, downstream when
                  the annotated stop is skipped and a later one is used
  junction_dist   bases between the end of the stop codon and the last exon junction of the form; negative
                  when the stop is in the last exon
  utr3            length of the transcript after the stop codon
  nmd             features that mark plant NMD targets: "junction" when the stop is at least --min-junction
                  bases upstream of the last exon junction (50-nucleotide rule), "long_utr3" when the
                  3' UTR is longer than --min-utr3 bases; empty when neither applies
The annotated form of the mRNA is checked first: its ORF must end exactly at the annotated stop codon,
otherwise the mRNA is not used (counted as inconsistent). Events whose junctions are absent from every
annotated mRNA of the gene, or whose start codon is not exonic in one form, are reported as unresolved.

ptc_form summarizes the event: 1 or 2 (which form carries the premature stop), both, none, or unresolved.
Together with the direction of form 1 in the test group (from the annotated events table) this tells
whether the PTC-containing form rises or falls under the condition.

Usage
  python nmd_features.py --events results_full/summary/events_high_confidence_annotated.tsv \
         --consequences results_full/summary/events_high_confidence_consequences.tsv \
         --gff ref/BeetSet-3.genes.2604.gff --genome ref/RefBeet-3.0.fa \
         --out results_full/summary/events_high_confidence_nmd.tsv
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from event_consequences import Fasta  # noqa: E402
from event_forms import Form, other_form, parse_event, read_models, start_and_stop, which_form  # noqa: E402


def describe(form, start_g, stop_g_annot, strand, min_junction, min_utr3):
    """Dictionary of the ORF features of one form (empty strings when the start codon is not exonic)."""
    i, j = form.orf(start_g)
    if i is None:
        return None
    if j is None:
        return dict(stop="none", ptc="no_stop", junction_dist="", utr3="", nmd="")
    stop_g = form.g_pos(j)
    annot_idx = form.tx_pos(stop_g_annot)
    if stop_g == stop_g_annot:
        ptc = "no"
    elif annot_idx is None or j < annot_idx:
        ptc = "yes"
    else:
        ptc = "downstream"
    dist = form.distance_to_last_junction(j)
    utr3 = form.length - (j + 3)
    feats = []
    if dist is not None and dist >= min_junction:
        feats.append("junction")
    if utr3 > min_utr3:
        feats.append("long_utr3")
    return dict(stop=stop_g, ptc=ptc, junction_dist=("" if dist is None else dist), utr3=utr3, nmd=";".join(feats))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True, help="events_high_confidence_annotated.tsv")
    ap.add_argument("--consequences", help="events_high_confidence_consequences.tsv, for the cross-check of ri_ptc")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-junction", type=int, default=50)
    ap.add_argument("--min-utr3", type=int, default=350)
    a = ap.parse_args()

    mrnas_of, info, exons, cds = read_models(a.gff)
    fasta = Fasta(a.genome)
    ri_ptc = {}
    if a.consequences:
        with open(a.consequences, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["event_type"] == "RI" and r["ri_ptc"]:
                    ri_ptc[r["suppa_event"]] = r["ri_ptc"]

    out = []
    inconsistent = set()
    checked = {}
    with open(a.events, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        chrom, strand, etype, f1, f2, regions = parse_event(r["suppa_event"])
        result, note = None, ""
        for gene in [g for g in r["ref_gene"].split(",") if g]:
            for m in mrnas_of.get(gene, []):
                if not cds.get(m) or info[m][0] != chrom or info[m][1] != strand:
                    continue
                form = which_form(exons[m], etype, f1, f2, regions)
                if form is None:
                    continue
                start_g, stop_g = start_and_stop(cds[m], strand)
                if m not in checked:
                    own = Form(chrom, strand, exons[m], fasta)
                    i, j = own.orf(start_g)
                    first_phase = cds[m][0][2] if strand == "+" else cds[m][-1][2]
                    checked[m] = first_phase == 0 and j is not None and own.g_pos(j) == stop_g
                if not checked[m]:
                    inconsistent.add(m)
                    note = "annotated ORF inconsistent"
                    continue
                other = other_form(exons[m], etype, form, f1, f2, regions)
                if other is None:
                    note = "other form not constructible"
                    continue
                chains = {form: exons[m], 3 - form: other}
                d1 = describe(Form(chrom, strand, chains[1], fasta), start_g, stop_g, strand, a.min_junction, a.min_utr3)
                d2 = describe(Form(chrom, strand, chains[2], fasta), start_g, stop_g, strand, a.min_junction, a.min_utr3)
                if d1 is None or d2 is None:
                    note = "start codon not exonic in one form"
                    continue
                result = (gene, m, form, d1, d2)
                break
            if result:
                break
        if result:
            gene, m, form, d1, d2 = result
            p1, p2 = d1["ptc"] == "yes", d2["ptc"] == "yes"
            ptc_form = "both" if p1 and p2 else "1" if p1 else "2" if p2 else "none"
            nmd_form = ";".join(k for k, d, p in (("1", d1, p1), ("2", d2, p2)) if p and d["nmd"])
            direction = r["direction"]
            ptc_dir = ""
            if ptc_form in ("1", "2") and direction in ("higher", "lower"):
                ptc_dir = direction if ptc_form == "1" else ("lower" if direction == "higher" else "higher")
            out.append(dict(
                suppa_event=r["suppa_event"], event_type=etype, gene=gene, mrna=m, annotated_form=form,
                form1_stop=d1["stop"], form1_ptc=d1["ptc"], form1_junction_dist=d1["junction_dist"], form1_utr3=d1["utr3"],
                form1_nmd=d1["nmd"], form2_stop=d2["stop"], form2_ptc=d2["ptc"], form2_junction_dist=d2["junction_dist"],
                form2_utr3=d2["utr3"], form2_nmd=d2["nmd"], ptc_form=ptc_form, nmd_form=nmd_form,
                form1_in_test=direction, ptc_form_in_test=ptc_dir, ri_ptc_consequences=ri_ptc.get(r["suppa_event"], ""),
                n_comparisons=r["n_comparisons"], groups=r["groups"], description=r["description"], note=""))
        else:
            out.append(dict(
                suppa_event=r["suppa_event"], event_type=etype, gene=r["ref_gene"], mrna="", annotated_form="",
                form1_stop="", form1_ptc="", form1_junction_dist="", form1_utr3="", form1_nmd="", form2_stop="",
                form2_ptc="", form2_junction_dist="", form2_utr3="", form2_nmd="", ptc_form="unresolved", nmd_form="",
                form1_in_test=r["direction"], ptc_form_in_test="", ri_ptc_consequences=ri_ptc.get(r["suppa_event"], ""),
                n_comparisons=r["n_comparisons"], groups=r["groups"], description=r["description"],
                note=note or ("no reference gene" if not r["ref_gene"] else "no annotated mRNA has the event junctions")))

    cols = list(out[0].keys())
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    print(f"events {len(out)}; mRNAs with an inconsistent annotated ORF {len(inconsistent)}")
    print("resolution and PTC form by event type:")
    for t in ("SE", "A3", "A5", "RI", "MX"):
        rows_t = [d for d in out if d["event_type"] == t]
        if not rows_t:
            continue
        c = collections.Counter(d["ptc_form"] for d in rows_t)
        n = collections.Counter(d["nmd_form"] for d in rows_t if d["ptc_form"] in ("1", "2"))
        print(f"  {t}: {len(rows_t)} events; ptc_form " + ", ".join(f"{k} {c[k]}" for k in ("1", "2", "both", "none", "unresolved") if c[k])
              + f"; PTC form with NMD features {sum(v for k, v in n.items() if k)} of {sum(n.values())}")
    res = [d for d in out if d["ptc_form"] in ("1", "2")]
    dirs = collections.Counter(d["ptc_form_in_test"] for d in res)
    print(f"events with a PTC in exactly one form: {len(res)}; that form higher in test {dirs['higher']}, "
          f"lower {dirs['lower']}, mixed or undefined {len(res) - dirs['higher'] - dirs['lower']}")
    if ri_ptc:
        ri = [d for d in out if d["event_type"] == "RI" and d["ri_ptc_consequences"] and d["ptc_form"] != "unresolved"]
        agree = sum(1 for d in ri if (d["form1_ptc"] == "yes") == (d["ri_ptc_consequences"] == "yes"))
        print(f"RI cross-check with event_consequences.py: retained form PTC agrees for {agree} of {len(ri)} events")
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
