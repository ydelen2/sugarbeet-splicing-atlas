"""Splicing-regulator gene families of sugar beet from the Pfam domain scan and the Arabidopsis searches.

Inputs (written by workflow/10_domains.sbatch)
  --domtbl     hmmsearch --domtblout of Pfam-A against the BeetSet-3 proteins (gathering thresholds)
  --proteins   the BeetSet-3 protein FASTA used for the scan
  --beet-ath   DIAMOND blastp beet versus Arabidopsis, --ath-beet the reverse search; both with the columns
               qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore
  --families   metadata/curation/splicing_factor_families.tsv: family, subfamily, uniprot_entry (entry name
               such as RSZ22_ARATH), gene, source. Membership of a beet protein in a family is given by its
               best Arabidopsis hit being one of the listed entries; reciprocal best hits are marked.
  --events     events_high_confidence_annotated.tsv (optional), to mark the family members with
               high-confidence differential splicing events
  --suppa-dir and --map (optional, with --events): the SUPPA2 diff directories and transcript_map.tsv give
               the background for an enrichment test, the reference genes with at least one testable event
               in any comparison. For every family the number of members in the background with a
               high-confidence event is compared with the background rate by a hypergeometric test
               (upper tail, exact).

SR proteins are defined by sequence with the plant definition of Barta, Kalyna and Reddy (2010, Plant Cell
22:2926): one or two N-terminal RRMs (PF00076) followed by a downstream RS domain of at least 50 residues
with at least 20 percent RS or SR dipeptides. The rule is applied to 50-residue windows downstream of the
last RRM: a window qualifies when at least 20 percent of its residues lie in RS or SR dipeptides, that is,
when it contains at least five dipeptides. This reading was calibrated on the 18 Arabidopsis SR proteins
of Barta et al. (2010), which all satisfy it (minimum six dipeptides in a window), whereas a threshold of
ten dipeptides per window excludes four of them (RS31, RS31a, RS40, RSZ21).
The highest arginine plus serine fraction over the same windows is reported as well, with the flag of the
metazoan definition of Manley and Krainer (2010, more than 40 percent), and the number of zinc knuckles
(Pfam zf-CCHC, PF00098) that distinguishes the RSZ (one) and RS2Z (two) subfamilies. Proteins with more
than two RRMs are not SR proteins by either definition. Whether a protein satisfies the rule is independent
of the Arabidopsis hit, so the two lines of evidence can be compared in the output.

Outputs
  <out-prefix>_proteins.tsv   one row per beet protein with an RRM, a family hit or an SR domain
  <out-prefix>_families.tsv   family members (genes), with the evidence and the DAS events
Only the best scoring protein of each gene is kept in the family table.

Usage
  python splicing_regulators.py --domtbl results_full/domains/pfam.domtblout \
      --proteins results_full/domains/beet_proteins.fa --beet-ath results_full/domains/beet_vs_ath.tsv \
      --ath-beet results_full/domains/ath_vs_beet.tsv --families metadata/curation/splicing_factor_families.tsv \
      --events results_full/summary/events_high_confidence_annotated.tsv \
      --out-prefix results_full/domains/splicing_regulators
"""
import argparse
import collections
import csv
import glob
import math
import os

RRM = "PF00076"
RS_WINDOW = 50
RS_MIN_DIPEPTIDES = 5           # 20 percent of the residues of a 50-residue window in RS or SR dipeptides (Barta et al. 2010)
RS_MIN_FRACTION = 0.40          # arginine plus serine fraction of the metazoan definition (Manley and Krainer 2010)
ZNF = "PF00098"                 # zf-CCHC zinc knuckle


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def read_domtbl(path):
    """{protein: [(pfam_name, pfam_acc, env_from, env_to, i_evalue)]} sorted by position."""
    doms = collections.defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            doms[f[0]].append((f[3], f[4].split(".")[0], int(f[19]), int(f[20]), float(f[12])))
    for p in doms:
        doms[p].sort(key=lambda d: d[2])
    return doms


def read_diamond(path):
    """{query: best hit row} by bitscore, and all hits per query."""
    best, hits = {}, collections.defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = dict(query=f[0], subject=f[1], pident=float(f[2]), length=int(f[3]), qlen=int(f[4]), slen=int(f[5]),
                       qcov=float(f[6]), scov=float(f[7]), evalue=float(f[8]), bitscore=float(f[9]))
            hits[f[0]].append(row)
            if f[0] not in best or row["bitscore"] > best[f[0]]["bitscore"]:
                best[f[0]] = row
    return best, hits


def entry_name(uniprot_id):
    """RSZ22_ARATH from sp|O81126|RSZ22_ARATH."""
    parts = uniprot_id.split("|")
    return parts[2] if len(parts) == 3 else uniprot_id


def rs_domain(seq, start):
    """Windows of 50 residues downstream of position start (0-based).

    Returns (max RS/SR dipeptides in a window, 1-based start of that window, max R+S fraction in a window).
    """
    tail = seq[start:]
    dip_max, dip_start, frac_max = 0, -1, 0.0
    for i in range(0, max(0, len(tail) - RS_WINDOW + 1)):
        w = tail[i:i + RS_WINDOW]
        dip = sum(1 for j in range(len(w) - 1) if w[j:j + 2] in ("RS", "SR"))
        if dip > dip_max:
            dip_max, dip_start = dip, start + i + 1
        frac_max = max(frac_max, sum(1 for c in w if c in "RS") / RS_WINDOW)
    return dip_max, dip_start, frac_max


def is_sr(rrms, dip_max):
    """Plant SR protein rule: one or two RRMs and a 50-residue window with at least five RS or SR dipeptides."""
    return 1 <= len(rrms) <= 2 and dip_max >= RS_MIN_DIPEPTIDES


def comparison_ids(path):
    """Comparison identifiers of a comparison table (first column comparison_id)."""
    with open(path, encoding="utf-8") as fh:
        return [r["comparison_id"] for r in csv.DictReader(fh, delimiter="\t")]


def testable_genes(suppa_dir, map_path, comparisons=None):
    """Reference genes of the merged genes with at least one testable SUPPA2 event in the comparisons of the
    given comparison table (all comparisons with a diff directory when no table is given). The table keeps the
    background of one comparison set independent of other sets analyzed in the same SUPPA2 directory."""
    ref_of = collections.defaultdict(set)
    with open(map_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["ref_gene"]:
                ref_of[r["gene"]].add(r["ref_gene"])
    merged = set()
    if comparisons:
        paths = [os.path.join(suppa_dir, "diff", c, "events.dpsi") for c in comparison_ids(comparisons)]
        missing = [p for p in paths if not os.path.isfile(p)]
        if missing:
            raise SystemExit(f"SUPPA2 results missing for {len(missing)} comparisons, first: {missing[0]}")
    else:
        paths = glob.glob(os.path.join(suppa_dir, "diff", "*", "events.dpsi"))
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                ev, d, p = line.rstrip("\n").split("\t")[:3]
                if d != "nan" and p != "nan":
                    for g in ev.split(";", 1)[0].split("_and_"):
                        merged.add(g)
    return {rg for g in merged for rg in ref_of.get(g, ())}


def hypergeom_upper(k, n, K, N):
    """P(X >= k) for X hypergeometric with n draws, K successes among N."""
    total = math.comb(N, n)
    return sum(math.comb(K, x) * math.comb(N - K, n - x) for x in range(k, min(n, K) + 1)) / total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--beet-ath", required=True)
    ap.add_argument("--ath-beet", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--events")
    ap.add_argument("--suppa-dir")
    ap.add_argument("--map")
    ap.add_argument("--comparisons", help="comparison table that defines the testable-gene background")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    seqs = read_fasta(a.proteins)
    doms = read_domtbl(a.domtbl)
    best_ba, hits_ba = read_diamond(a.beet_ath)
    best_ab, _ = read_diamond(a.ath_beet)

    fam_of = {}
    with open(a.families, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam_of[r["uniprot_entry"]] = r

    das = collections.defaultdict(list)
    if a.events:
        with open(a.events, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for g in r["ref_gene"].split(","):
                    if g:
                        das[g].append(f"{r['event_type']}:{r['direction']}:{r['n_comparisons']}:{r['groups']}")

    rows = []
    for prot, seq in seqs.items():
        d = doms.get(prot, [])
        rrms = [x for x in d if x[1] == RRM]
        arch = ";".join(f"{x[0]}({x[2]}-{x[3]})" for x in d)
        dip, wstart, frac = rs_domain(seq, rrms[-1][3]) if rrms else (0, -1, 0.0)
        sr_rule = is_sr(rrms, dip)
        n_znf = sum(1 for x in d if x[1] == ZNF)
        b = best_ba.get(prot)
        ent = entry_name(b["subject"]) if b else ""
        fam = fam_of.get(ent)
        rbh = bool(b) and best_ab.get(b["subject"], {}).get("subject") == prot
        # any listed entry among the hits, when the best hit is not a listed one
        other = ""
        if not fam:
            for h in hits_ba.get(prot, []):
                e = entry_name(h["subject"])
                if e in fam_of:
                    other = f"{e} (bitscore {h['bitscore']:.0f}, best {b['bitscore']:.0f})"
                    break
        if not (rrms or fam or sr_rule or other):
            continue
        gene = prot.rsplit(".", 1)[0]
        rows.append(dict(
            protein=prot, gene=gene, length=len(seq), n_rrm=len(rrms), domains=arch,
            rs_dipeptides_max=dip, rs_window_start=wstart, sr_by_sequence="yes" if sr_rule else "no",
            rs_fraction_max=round(frac, 3), sr_by_rs_content="yes" if (1 <= len(rrms) <= 2 and frac > RS_MIN_FRACTION) else "no",
            n_zinc_knuckles=n_znf,
            best_ath=ent, best_pident=b["pident"] if b else "", best_qcov=b["qcov"] if b else "", best_scov=b["scov"] if b else "",
            best_bitscore=b["bitscore"] if b else "", reciprocal_best="yes" if rbh else "no",
            family=fam["family"] if fam else "", subfamily=fam["subfamily"] if fam else "",
            ath_gene=fam["gene"] if fam else "", listed_entry_among_hits=other,
            das_events=";".join(das.get(gene, []))))

    rows.sort(key=lambda r: (r["family"] or "~", r["subfamily"], r["gene"]))
    cols = list(rows[0].keys()) if rows else []
    with open(a.out_prefix + "_proteins.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # family table: one row per gene, best scoring protein of the gene
    by_gene = {}
    for r in rows:
        if not r["family"]:
            continue
        g = r["gene"]
        if g not in by_gene or float(r["best_bitscore"]) > float(by_gene[g]["best_bitscore"]):
            by_gene[g] = r
    fam_rows = sorted(by_gene.values(), key=lambda r: (r["family"], r["subfamily"], r["ath_gene"], r["gene"]))
    with open(a.out_prefix + "_families.tsv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(fam_rows)

    print(f"proteins with an RRM: {sum(1 for r in rows if r['n_rrm'])}; genes: {len({r['gene'] for r in rows if r['n_rrm']})}")
    sr = [r for r in rows if r["sr_by_sequence"] == "yes"]
    print(f"SR proteins by the sequence rule: {len(sr)} proteins, {len({r['gene'] for r in sr})} genes")
    sr_fam = {r["gene"] for r in fam_rows if r["family"] == "SR"}
    sr_seq = {r["gene"] for r in sr}
    print(f"SR genes by Arabidopsis best hit: {len(sr_fam)}; by both lines of evidence: {len(sr_fam & sr_seq)}; "
          f"rule only: {len(sr_seq - sr_fam)}; hit only: {len(sr_fam - sr_seq)}")
    fam_count = collections.Counter(r["family"] for r in fam_rows)
    das_count = collections.Counter(r["family"] for r in fam_rows if r["das_events"])
    rbh_count = collections.Counter(r["family"] for r in fam_rows if r["reciprocal_best"] == "yes")
    print("family\tgenes\treciprocal_best\twith_DAS_event")
    for f in sorted(fam_count):
        print(f"{f}\t{fam_count[f]}\t{rbh_count[f]}\t{das_count[f]}")

    if a.events and a.suppa_dir and a.map:
        bg = testable_genes(a.suppa_dir, a.map, a.comparisons)
        das_bg = {g for g in das if g in bg}
        print(f"\nenrichment background: {len(bg)} reference genes with a testable event; {len(das_bg)} with a high-confidence event "
              f"({100 * len(das_bg) / len(bg):.1f} percent)")
        print("family\tmembers_in_background\twith_event\texpected\tp_hypergeometric")
        groups = collections.defaultdict(set)
        for r in fam_rows:
            groups[r["family"]].add(r["gene"])
        groups["all families"] = set().union(*groups.values()) if groups else set()
        enr = []
        for f in sorted(groups):
            members = groups[f] & bg
            k = len(members & das_bg)
            if not members:
                continue
            p = hypergeom_upper(k, len(members), len(das_bg), len(bg))
            exp = len(members) * len(das_bg) / len(bg)
            print(f"{f}\t{len(members)}\t{k}\t{exp:.1f}\t{p:.3g}")
            enr.append(dict(family=f, members_in_background=len(members), with_event=k, expected=round(exp, 2),
                            p_hypergeometric=p, background=len(bg), background_with_event=len(das_bg)))
        with open(a.out_prefix + "_enrichment.tsv", "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=list(enr[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(enr)
    print(f"tables: {a.out_prefix}_proteins.tsv  {a.out_prefix}_families.tsv")


if __name__ == "__main__":
    main()
