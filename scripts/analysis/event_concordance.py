"""Event-level agreement between SUPPA2 and rMATS.

Event types compared: SE (rMATS SE), A3 (A3SS), A5 (A5SS), RI (RI), MX (MXE). SUPPA2 AF and AL
events have no rMATS counterpart.

Matching. The coordinates in a SUPPA2 event id are splice sites in 1-based GTF coordinates. An rMATS
row is converted to the same system (columns ending in Start_0base or ES are 0-based starts and get +1,
columns ending in End or EE are 1-based ends). A SUPPA2 event matches an rMATS row of the same type on
the same sequence and strand when every SUPPA2 coordinate is in the rMATS coordinate set. For RI only
the two intron boundaries of the SUPPA2 id are used, because the outer exon boundaries depend on the
transcript. When several rMATS rows match one SUPPA2 event, the covered row with the lowest FDR is
used and the case is counted in rmats_rows_ambiguous. rMATS writes sequence names with a "chr" prefix
(BvchrU.ctg1145_bp24838 becomes chrBvchrU.ctg1145_bp24838); the prefix is removed when the SUPPA2 file
uses the bare name. For SE events the rule is rMATS exonStart_0base = SUPPA2 s2 - 1,
exonEnd = e2, upstreamEE = e1, downstreamES = s3 - 1, on both strands.

Directions. Both tools are reported as test minus reference; rMATS IncLevelDifference is reference
minus test and is negated. For MX events on the minus strand the two tools report the inclusion of
opposite exons (rMATS always calls the genomically first exon the first one; in the sugar beet data the
per strand correlation of the two dPSI values was 0.72, 0.80 and 0.22 on the plus strand and -0.63,
-0.70 and -0.63 on the minus strand in three comparisons), so the rMATS dPSI of minus-strand MX events
is negated once more. The per type Pearson correlation over the matched testable events (dpsi_pearson)
is the check of these conventions: a negative value would mean that the two tools define inclusion
differently for that type.

Significance. SUPPA2: p < --pval (nominal) and, separately, Benjamini-Hochberg adjusted p < --fdr over
the testable events of the same type in the comparison, both with |dPSI| >= --dpsi. rMATS: FDR < --fdr,
|dPSI| >= --dpsi and every replicate with >= --min-reads supporting reads.

Outputs
  <out-prefix>_summary.tsv          one row per comparison and event type
  <out-prefix>_high_confidence.tsv  events significant in both tools (SUPPA2 nominal p) with the same direction

Usage
  python event_concordance.py --comparisons comparisons_all.tsv --suppa-dir results_full/suppa \
         --rmats-dir results_full/rmats_annot --map results_full/filtered/transcript_map.tsv \
         --out-prefix results_full/summary/event_concordance
"""
import argparse
import collections
import csv
import math
import os

from stats_util import hypergeom_upper

TYPES = {"SE": "SE", "A3": "A3SS", "A5": "A5SS", "RI": "RI", "MX": "MXE"}
RMATS_COORDS = {
    "SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}


def fnum(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values, same order as the input."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvals[i] * n / rank)
        adj[i] = running
    return adj


def suppa_coords(etype, body):
    """Splice-site coordinates of a SUPPA2 event id body (the part after the type)."""
    parts = body.split(":")
    chrom, strand, fields = parts[0], parts[-1], parts[1:-1]
    if etype == "RI":
        fields = [f for f in fields if "-" in f]      # keep e1-s2, drop the outer s1 and e2
    coords = frozenset(int(x) for f in fields for x in f.split("-"))
    return chrom, strand, coords


def read_suppa(path, etype):
    """{event_id: (gene, chrom, strand, coords, dpsi, p)} for testable events of one type."""
    out = {}
    prefix = etype + ":"
    with open(path, encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            ev, d, p = line.rstrip("\n").split("\t")[:3]
            gene, body = ev.split(";", 1)
            if not body.startswith(prefix):
                continue
            d, p = fnum(d), fnum(p)
            if d is None or p is None:
                continue
            chrom, strand, coords = suppa_coords(etype, body[len(prefix):])
            out[ev] = (gene, chrom, strand, coords, d, p)
    return out


def read_rmats(path, rtype, min_reads, seqnames):
    """Rows indexed by (chrom, strand, coordinate): [(coords, dpsi_test_minus_ref, fdr, covered, row_id)]."""
    index = collections.defaultdict(list)
    n = 0
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            chrom = r["chr"]
            if chrom not in seqnames and chrom.startswith("chr") and chrom[3:] in seqnames:
                chrom = chrom[3:]
            coords = set()
            for c in RMATS_COORDS[rtype]:
                v = int(r[c])
                coords.add(v + 1 if c.endswith("Start_0base") or c.endswith("ES") else v)
            counts = []
            for g in ("1", "2"):
                counts += [(fnum(a) or 0) + (fnum(b) or 0) for a, b in
                           zip(r[f"IJC_SAMPLE_{g}"].split(","), r[f"SJC_SAMPLE_{g}"].split(","))]
            d, fdr = fnum(r["IncLevelDifference"]), fnum(r["FDR"])
            covered = bool(counts) and min(counts) >= min_reads and d is not None and fdr is not None
            entry = (frozenset(coords), None if d is None else -d, fdr, covered, n)
            for v in coords:
                index[(chrom, r["strand"], v)].append(entry)
    return index, n


def match(index, chrom, strand, coords):
    """rMATS rows whose coordinate set contains every SUPPA2 coordinate."""
    anchor = min(coords)
    return [e for e in index.get((chrom, strand, anchor), []) if coords <= e[0]]


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return None if sxx == 0 or syy == 0 else sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--suppa-dir", required=True)
    ap.add_argument("--rmats-dir", required=True)
    ap.add_argument("--map", required=True, help="transcript_map.tsv (merged gene to reference gene)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--types", default=",".join(TYPES), help="SUPPA2 event types to compare")
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--dpsi", type=float, default=0.1)
    ap.add_argument("--min-reads", type=int, default=20)
    a = ap.parse_args()

    ref_genes = collections.defaultdict(set)
    with open(a.map, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ref_gene"]:
                ref_genes[r["gene"]].add(r["ref_gene"])

    with open(a.comparisons, encoding="utf-8") as fh:
        comps = list(csv.DictReader(fh, delimiter="\t"))

    summary, hc = [], []
    for c in comps:
        cid = c["comparison_id"]
        sp = os.path.join(a.suppa_dir, "diff", cid, "events.dpsi")
        if not os.path.isfile(sp):
            print(f"{cid}: SUPPA2 file missing, skipped")
            continue
        for etype in a.types.split(","):
            rtype = TYPES[etype]
            rp = os.path.join(a.rmats_dir, cid, f"{rtype}.MATS.JCEC.txt")
            if not os.path.isfile(rp):
                print(f"{cid} {rtype}: rMATS file missing, skipped")
                continue
            s = read_suppa(sp, etype)
            index, n_rows = read_rmats(rp, rtype, a.min_reads, {v[1] for v in s.values()})
            evs = list(s)
            adj = dict(zip(evs, bh([s[e][5] for e in evs])))
            partner, ambiguous = {}, 0
            for e in evs:
                _, chrom, strand, coords, _, _ = s[e]
                rows = match(index, chrom, strand, coords)
                if not rows:
                    continue
                if len(rows) > 1:
                    ambiguous += 1
                    rows.sort(key=lambda x: (not x[3], x[2] if x[2] is not None else 2.0))
                partner[e] = rows[0]
            def rd(e):
                """rMATS dPSI of the partner of SUPPA2 event e, as test minus reference of the SUPPA2 form."""
                return -partner[e][1] if etype == "MX" and s[e][2] == "-" else partner[e][1]
            covered_rows = {r[4] for rs in index.values() for r in rs if r[3]}
            matched = {e for e, r in partner.items() if r[3]}
            s_sig = {e for e in evs if abs(s[e][4]) >= a.dpsi and s[e][5] < a.pval}
            s_sig_bh = {e for e in evs if abs(s[e][4]) >= a.dpsi and adj[e] < a.fdr}
            r_sig_rows = {r[4] for rs in index.values() for r in rs if r[3] and abs(r[1]) >= a.dpsi and r[2] < a.fdr}
            both = {e for e in s_sig & matched if partner[e][4] in r_sig_rows}
            same = {e for e in both if (s[e][4] > 0) == (rd(e) > 0)}
            both_bh = {e for e in s_sig_bh & matched if partner[e][4] in r_sig_rows}
            s_sig_matched = s_sig & matched
            r_sig_matched = {e for e in matched if partner[e][4] in r_sig_rows}
            n_m = len(matched)
            expected = len(s_sig_matched) * len(r_sig_matched) / n_m if n_m else 0.0
            p_overlap = hypergeom_upper(len(both), len(s_sig_matched), len(r_sig_matched), n_m)
            rho = pearson([s[e][4] for e in matched], [rd(e) for e in matched])
            up = sum(1 for e in same if s[e][4] > 0)
            summary.append(dict(
                comparison=cid, event_type=etype, suppa_testable=len(s), rmats_all=n_rows,
                rmats_covered=len(covered_rows),
                suppa_with_rmats_partner_pct=round(100 * len(partner) / len(s), 1) if s else 0,
                rmats_rows_ambiguous=ambiguous, matched_testable=len(matched),
                dpsi_pearson=("" if rho is None else round(rho, 3)),
                suppa_sig_nominal=len(s_sig), suppa_sig_bh=len(s_sig_bh), rmats_sig=len(r_sig_rows),
                sig_both=len(both), sig_both_same_direction=len(same), sig_both_bh=len(both_bh),
                same_direction_higher_in_test=up, same_direction_lower_in_test=len(same) - up,
                suppa_sig_matched=len(s_sig_matched), rmats_sig_matched=len(r_sig_matched),
                sig_both_expected=round(expected, 2),
                sig_both_enrichment=("" if expected == 0 else round(len(both) / expected, 1)),
                sig_both_p=f"{p_overlap:.3g}"))
            for e in sorted(same):
                gene, _, _, _, d, p = s[e]
                hc.append(dict(comparison=cid, event_type=etype, suppa_event=e, gene=gene,
                               ref_genes=",".join(sorted(ref_genes.get(gene, []))),
                               suppa_dpsi=round(d, 4), suppa_p=p, suppa_p_bh=round(adj[e], 6),
                               rmats_dpsi=round(rd(e), 4), rmats_fdr=partner[e][2]))

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    cols = list(summary[0].keys()) if summary else []
    with open(a.out_prefix + "_summary.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# SUPPA2 nominal p < {a.pval} or BH < {a.fdr}; rMATS FDR < {a.fdr}, every replicate >= "
                 f"{a.min_reads} reads; |dPSI| >= {a.dpsi}; directions are test minus reference; the _matched "
                 "columns count significant events among the matched testable events, sig_both_expected is the "
                 "overlap expected by chance for independent calls (product of the two matched counts over "
                 "matched_testable) and sig_both_p its hypergeometric upper-tail probability\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    hcols = ["comparison", "event_type", "suppa_event", "gene", "ref_genes", "suppa_dpsi", "suppa_p", "suppa_p_bh",
             "rmats_dpsi", "rmats_fdr"]
    with open(a.out_prefix + "_high_confidence.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=hcols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(hc)

    print("\t".join(cols))
    for row in summary:
        print("\t".join(str(row[k]) for k in cols))
    by_type = collections.Counter(r["event_type"] for r in hc)
    print("\nhigh-confidence events (both tools, same direction): " + str(len(hc)) + "; by type: "
          + ", ".join(f"{t} {by_type[t]}" for t in TYPES if by_type[t]))
    print(f"tables: {a.out_prefix}_summary.tsv  {a.out_prefix}_high_confidence.tsv")


if __name__ == "__main__":
    main()
