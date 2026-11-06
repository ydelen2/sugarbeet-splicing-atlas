"""Coding consequences of the high-confidence splicing events.

For every event the alternative region is derived from the SUPPA2 event id:
  SE   the skipped exon
  A3   the segment between the two 3' splice sites (part of the long exon form)
  A5   the segment between the two 5' splice sites
  RI   the retained intron
  MX   the two mutually exclusive exons (both regions)
In a SUPPA2 id every "x-y" pair is a junction written in genomic order, so x is an exon end and y an exon
start on either strand; the two non-shared coordinates of an A3 or A5 event are therefore both exon ends
(region x1+1..x2) or both exon starts (region y1..y2-1).

The region is placed on the reference gene models (GFF3, mRNA with CDS children):
  cds_context   CDS        the region lies inside the CDS span of at least one mRNA of the gene
                UTR        inside the span of an mRNA that has a CDS, but outside every CDS span
                noncoding  inside an mRNA span, but no mRNA of the gene has a CDS
                none       no reference gene, or the region is not inside any mRNA span of the gene
  frame         region length modulo 3 (MX: difference of the two exon lengths modulo 3); 0 means the
                inclusion or exclusion of the region preserves the reading frame
  ri_ptc        RI events inside a CDS only: yes when the retained intron, read in the frame of the
                flanking CDS, contains a stop codon (TAA, TAG or TGA); the codon that spans the exon-intron
                boundary is included, the codon that spans the intron-exon boundary is not
The frame at the intron boundary comes from the CDS feature that ends exactly at the intron (plus strand)
or starts exactly after it (minus strand) in one mRNA of the gene; if no CDS feature borders the intron,
ri_ptc is left empty and ri_frame_source says why.

Genome access reads the FASTA through its .fai index (samtools faidx layout), no extra library needed.

Usage
  python event_consequences.py --events results_full/summary/events_high_confidence_annotated.tsv \
         --gff ref/BeetSet-3.genes.2604.gff --genome ref/RefBeet-3.0.fa \
         --out results_full/summary/events_high_confidence_consequences.tsv
"""
import argparse
import collections
import csv
import re

STOPS = {"TAA", "TAG", "TGA"}
RANK = {"none": 0, "noncoding": 1, "UTR": 2, "CDS": 3}   # the most informative context wins for multi-gene loci
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Fasta:
    """Random access to a FASTA file through its .fai index; 1-based inclusive coordinates."""

    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai", encoding="utf-8") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")[:5]
                self.idx[name] = (int(length), int(offset), int(linebases), int(linewidth))

    def fetch(self, name, start, end):
        length, offset, linebases, linewidth = self.idx[name]
        start, end = max(1, start), min(end, length)
        if end < start:
            return ""
        first = offset + (start - 1) // linebases * linewidth + (start - 1) % linebases
        last = offset + (end - 1) // linebases * linewidth + (end - 1) % linebases
        self.fh.seek(first)
        raw = self.fh.read(last - first + 1)
        return raw.decode("ascii").replace("\n", "").replace("\r", "").upper()


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_gff(path):
    """Gene models: {gene: [mrna]}, {mrna: (chrom, strand, start, end)}, {mrna: [(cds_start, cds_end, phase)]}."""
    mrnas_of = collections.defaultdict(list)
    mrna = {}
    cds = collections.defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "mRNA":
                mrna[attrs["ID"]] = (f[0], f[6], int(f[3]), int(f[4]))
                mrnas_of[attrs.get("Parent", "")].append(attrs["ID"])
            elif f[2] == "CDS":
                for parent in attrs.get("Parent", "").split(","):
                    cds[parent].append((int(f[3]), int(f[4]), int(f[7]) if f[7] in "012" else 0))
    return mrnas_of, mrna, cds


def regions(ev):
    """(chrom, strand, event type, [(start, end), ...]) of the alternative region(s)."""
    body = ev.split(";", 1)[1]
    parts = body.split(":")
    etype, chrom, strand, fields = parts[0], parts[1], parts[-1], parts[2:-1]
    if etype == "SE":
        e1, s2 = map(int, fields[0].split("-"))
        e2, s3 = map(int, fields[1].split("-"))
        return chrom, strand, etype, [(s2, e2)]
    if etype == "RI":
        e1, s2 = map(int, fields[1].split("-"))
        return chrom, strand, etype, [(e1 + 1, s2 - 1)]
    if etype in ("A3", "A5"):
        (x1, y1), (x2, y2) = (tuple(map(int, f.split("-"))) for f in fields)
        if y1 == y2:                       # shared exon start, the alternative sites are exon ends
            a, b = sorted((x1, x2))
            return chrom, strand, etype, [(a + 1, b)]
        a, b = sorted((y1, y2))            # shared exon end, the alternative sites are exon starts
        return chrom, strand, etype, [(a, b - 1)]
    if etype == "MX":
        e1, s2 = map(int, fields[0].split("-"))
        e2, s4 = map(int, fields[1].split("-"))
        _, s3 = map(int, fields[2].split("-"))
        e3, _ = map(int, fields[3].split("-"))
        return chrom, strand, etype, [(s2, e2), (s3, e3)]
    return chrom, strand, etype, []


def cds_context(gene, regs, mrnas_of, mrna, cds):
    inside_coding, inside_noncoding, inside_cds = False, False, False
    for m in mrnas_of.get(gene, []):
        _, _, ms, me = mrna[m]
        if all(ms <= a and b <= me for a, b in regs):
            if cds.get(m):
                inside_coding = True
                cs, ce = min(c[0] for c in cds[m]), max(c[1] for c in cds[m])
                if all(cs <= a and b <= ce for a, b in regs):
                    inside_cds = True
            else:
                inside_noncoding = True
    if inside_cds:
        return "CDS"
    if inside_coding:
        return "UTR"
    return "noncoding" if inside_noncoding else "none"


def ri_ptc(gene, chrom, strand, intron, mrnas_of, cds, fasta):
    """(ptc yes/no/'', frame source) for a retained intron inside a CDS."""
    a, b = intron
    for m in mrnas_of.get(gene, []):
        for cs, ce, ph in cds.get(m, []):
            if strand == "+" and ce == a - 1:
                k = (a - cs - ph) % 3               # exon bases already in the codon at the boundary
                seq = fasta.fetch(chrom, a - k, b)
                break
            if strand == "-" and cs == b + 1:
                k = (ce - ph - b) % 3
                seq = revcomp(fasta.fetch(chrom, a, b + k))
                break
        else:
            continue
        codons = {seq[i:i + 3] for i in range(0, len(seq) - 2, 3)}
        return ("yes" if codons & STOPS else "no"), m
    return "", "no CDS feature borders the intron"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True, help="events_high_confidence_annotated.tsv")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True, help="FASTA with a .fai index next to it")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    mrnas_of, mrna, cds = read_gff(a.gff)
    fasta = Fasta(a.genome)
    print(f"gene models: {len(mrnas_of)} genes, {len(mrna)} mRNAs, {sum(len(v) for v in cds.values())} CDS features")

    out = []
    with open(a.events, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            chrom, strand, etype, regs = regions(r["suppa_event"])
            genes = [g for g in r["ref_gene"].split(",") if g]
            lengths = [b - a + 1 for a, b in regs]
            frame = (abs(lengths[0] - lengths[1]) if etype == "MX" else lengths[0]) % 3 if regs else ""
            ctx = ""
            for g in genes:
                c = cds_context(g, regs, mrnas_of, mrna, cds)
                if RANK[c] > RANK.get(ctx, -1):
                    ctx = c
            ptc, src = "", ""
            if etype == "RI" and ctx == "CDS":
                for g in genes:
                    ptc, src = ri_ptc(g, chrom, strand, regs[0], mrnas_of, cds, fasta)
                    if ptc:
                        break
            out.append(dict(
                suppa_event=r["suppa_event"], event_type=etype, ref_gene=r["ref_gene"],
                region=";".join(f"{chrom}:{s}-{e}" for s, e in regs), region_length=";".join(map(str, lengths)),
                cds_context=ctx if genes else "none", frame=frame, ri_ptc=ptc, ri_frame_source=src,
                n_comparisons=r["n_comparisons"], groups=r["groups"], direction=r["direction"],
                description=r["description"]))

    cols = list(out[0].keys()) if out else []
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    print("event type, CDS context and frame (frame 0 = length divisible by 3):")
    tab = collections.Counter((d["event_type"], d["cds_context"], d["frame"]) for d in out)
    for k in sorted(tab):
        print("  " + "\t".join(map(str, k)) + f"\t{tab[k]}")
    ri = collections.Counter(d["ri_ptc"] for d in out if d["event_type"] == "RI" and d["cds_context"] == "CDS")
    print(f"RI inside CDS: stop codon in the retained intron yes {ri['yes']}, no {ri['no']}, undetermined {ri['']}")
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
