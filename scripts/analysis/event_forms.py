"""Transcript forms of SUPPA2 events on the reference gene models, shared by nmd_features.py and
domain_impact.py.

A SUPPA2 event id lists junctions as "x-y" pairs in genomic order (x an exon end, y an exon start):
  SE  e1-s2:e2-s3            form 1 includes the exon s2..e2, form 2 joins e1 to s3
  RI  s1:e1-s2:e2            form 1 retains the intron e1+1..s2-1, form 2 splices it (junction e1-s2)
  A3, A5  x1-y1:x2-y2        form 1 uses the first junction, form 2 the second
  MX  e1-s2:e2-s4:e1-s3:e3-s4   form 1 includes the exon s2..e2, form 2 the exon s3..e3
Form 1 is the form whose inclusion SUPPA2 reports as PSI (the alternative_transcripts of the ioe file), so a
positive dPSI means form 1 is more used in the test group.

An annotated mRNA "is" one form when its exon chain contains every junction of that form (RI form 1: one
exon spans the intron). The other form is then built from the same exon chain by removing or inserting the
skipped exon, merging or splitting the exons around the retained intron, moving one splice site, or swapping
the mutually exclusive exons. Everything downstream is derived from the two exon chains: transcript sequence,
position of the annotated start codon, first in-frame stop codon, exon junctions and 3' UTR.
"""
import collections

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_models(path):
    """({gene: [mrna]}, {mrna: (chrom, strand)}, {mrna: [(exon_start, exon_end)]}, {mrna: [(cds_start, cds_end, phase)]})
    from a GFF3 file; exon and CDS lists are sorted by genomic start."""
    mrnas_of = collections.defaultdict(list)
    info, exons, cds = {}, collections.defaultdict(list), collections.defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "mRNA":
                info[attrs["ID"]] = (f[0], f[6])
                mrnas_of[attrs.get("Parent", "")].append(attrs["ID"])
            elif f[2] in ("exon", "CDS"):
                for parent in attrs.get("Parent", "").split(","):
                    if f[2] == "exon":
                        exons[parent].append((int(f[3]), int(f[4])))
                    else:
                        cds[parent].append((int(f[3]), int(f[4]), int(f[7]) if f[7] in "012" else 0))
    for m in exons:
        exons[m].sort()
    for m in cds:
        cds[m].sort()
    return mrnas_of, info, exons, cds


def parse_event(ev):
    """(chrom, strand, type, form1_junctions, form2_junctions, regions) of a SUPPA2 event id."""
    body = ev.split(";", 1)[1]
    parts = body.split(":")
    etype, chrom, strand, fields = parts[0], parts[1], parts[-1], parts[2:-1]
    pair = lambda f: tuple(map(int, f.split("-")))
    if etype == "SE":
        (e1, s2), (e2, s3) = pair(fields[0]), pair(fields[1])
        return chrom, strand, etype, [(e1, s2), (e2, s3)], [(e1, s3)], [(s2, e2)]
    if etype == "RI":
        e1, s2 = pair(fields[1])
        return chrom, strand, etype, [], [(e1, s2)], [(e1 + 1, s2 - 1)]
    if etype in ("A3", "A5"):
        j1, j2 = pair(fields[0]), pair(fields[1])
        if j1[1] == j2[1]:
            a, b = sorted((j1[0], j2[0]))
            region = (a + 1, b)
        else:
            a, b = sorted((j1[1], j2[1]))
            region = (a, b - 1)
        return chrom, strand, etype, [j1], [j2], [region]
    if etype == "MX":
        (e1, s2), (e2, s4), (_, s3), (e3, _) = (pair(f) for f in fields)
        return chrom, strand, etype, [(e1, s2), (e2, s4)], [(e1, s3), (e3, s4)], [(s2, e2), (s3, e3)]
    raise ValueError(f"unknown event type in {ev}")


def junctions(exons):
    return {(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)}


def valid_chain(exons):
    return all(s <= e for s, e in exons) and all(exons[i][1] < exons[i + 1][0] for i in range(len(exons) - 1))


def which_form(exons, etype, f1, f2, regions):
    """1 or 2 when the exon chain is that form of the event, else None."""
    J = junctions(exons)
    if etype == "RI":
        a, b = regions[0]
        if any(s < a and e > b for s, e in exons):
            return 1
        return 2 if set(f2) <= J else None
    if set(f1) <= J:
        return 1
    if set(f2) <= J:
        return 2
    return None


def other_form(exons, etype, form, f1, f2, regions):
    """Exon chain of the other form of the event, or None when it cannot be built from this chain."""
    ex = list(exons)
    if etype == "SE":
        s2, e2 = regions[0]
        if form == 1:
            if (s2, e2) not in ex:
                return None
            ex.remove((s2, e2))
        else:
            ex.append((s2, e2))
            ex.sort()
    elif etype == "RI":
        a, b = regions[0]
        if form == 1:
            for i, (s, e) in enumerate(ex):
                if s < a and e > b:
                    ex[i:i + 1] = [(s, a - 1), (b + 1, e)]
                    break
        else:
            for i in range(len(ex) - 1):
                if ex[i][1] == a - 1 and ex[i + 1][0] == b + 1:
                    ex[i:i + 2] = [(ex[i][0], ex[i + 1][1])]
                    break
    elif etype in ("A3", "A5"):
        (x_old, y_old), (x_new, y_new) = (f1[0], f2[0]) if form == 1 else (f2[0], f1[0])
        for i, (s, e) in enumerate(ex):
            if x_old == x_new and s == y_old:          # shared exon end, the exon start moves
                ex[i] = (y_new, e)
                break
            if y_old == y_new and e == x_old:          # shared exon start, the exon end moves
                ex[i] = (s, x_new)
                break
        else:
            return None
    elif etype == "MX":
        old, new = (regions[0], regions[1]) if form == 1 else (regions[1], regions[0])
        if old not in ex:
            return None
        ex[ex.index(old)] = new
        ex.sort()
    return ex if valid_chain(ex) else None


class Form:
    """Transcript of one exon chain: sequence, coordinate mapping and the ORF from a given start codon."""

    def __init__(self, chrom, strand, exons, fasta):
        self.strand = strand
        self.exons = list(exons) if strand == "+" else list(reversed(exons))   # transcription order
        self.seq = "".join(fasta.fetch(chrom, s, e) if strand == "+" else revcomp(fasta.fetch(chrom, s, e))
                           for s, e in self.exons)
        self.offsets = []                                     # transcript index of the first base of each exon
        pos = 0
        for s, e in self.exons:
            self.offsets.append(pos)
            pos += e - s + 1
        self.length = pos
        self.junction_positions = self.offsets[1:]            # index of the first base after each junction

    def tx_pos(self, gpos):
        """0-based transcript index of a genomic position, or None when it is not exonic."""
        for (s, e), off in zip(self.exons, self.offsets):
            if s <= gpos <= e:
                return off + (gpos - s if self.strand == "+" else e - gpos)
        return None

    def g_pos(self, idx):
        for (s, e), off in zip(self.exons, self.offsets):
            if off <= idx < off + e - s + 1:
                return s + (idx - off) if self.strand == "+" else e - (idx - off)
        return None

    def orf(self, start_gpos, stops=("TAA", "TAG", "TGA")):
        """(start index, stop index or None) of the ORF starting at the given genomic start codon position."""
        i = self.tx_pos(start_gpos)
        if i is None:
            return None, None
        for j in range(i, self.length - 2, 3):
            if self.seq[j:j + 3] in stops:
                return i, j
        return i, None

    def distance_to_last_junction(self, stop_idx):
        """Bases between the end of the stop codon and the last exon junction; negative inside the last exon."""
        if not self.junction_positions:
            return None
        return self.junction_positions[-1] - (stop_idx + 3)


def start_and_stop(cds_list, strand):
    """Genomic positions of the first base of the start codon and of the stop codon of an annotated CDS."""
    if strand == "+":
        return cds_list[0][0], cds_list[-1][1] - 2
    return cds_list[-1][1], cds_list[0][0] + 2
