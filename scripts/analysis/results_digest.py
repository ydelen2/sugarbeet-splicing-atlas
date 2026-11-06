"""Digest of the summary tables: the numbers quoted in the manuscript, printed from the result files.

Reads the tables written by 09_summarize.sh and prints, in order: per comparison and per condition group
counts of the high-confidence events by type and direction; recurrence of events across comparisons;
consequences (location, frame, stop codons in retained introns); splicing-regulator families and their
enrichment; GO terms per group (FDR below --fdr, at most --top per group, cellular component excluded);
the merge and filter statistics from the log files when present.

Usage
  python results_digest.py --summary results_full/summary --domains results_full/domains \
      --comparisons metadata/comparisons_all.tsv --fdr 0.05 --top 10
"""
import argparse
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import binom_two_sided  # noqa: E402


def read_tsv(path, skip_comment=True):
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if not (skip_comment and l.startswith("#"))]
    return list(csv.DictReader(lines, delimiter="\t"))


def group_of(cid):
    return re.sub(r"\d+C$", "", cid.split("_")[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--expression", help="directory of 12_expression.sbatch (de_summary.tsv), optional")
    a = ap.parse_args()
    S = a.summary

    comps = read_tsv(a.comparisons)
    print(f"comparisons: {len(comps)}")
    by_group = collections.Counter(group_of(c["comparison_id"]) for c in comps)
    print("comparisons per group: " + ", ".join(f"{g} {n}" for g, n in sorted(by_group.items())))

    print("\n== event concordance per comparison (all types pooled) ==")
    conc = read_tsv(os.path.join(S, "event_concordance_summary.tsv"))
    per_comp = collections.defaultdict(lambda: collections.Counter())
    per_type = collections.defaultdict(lambda: collections.Counter())
    per_group_type = collections.defaultdict(lambda: collections.Counter())
    rho = collections.defaultdict(list)
    for r in conc:
        c, t = r["comparison"], r["event_type"]
        for k in ("suppa_testable", "matched_testable", "suppa_sig_nominal", "suppa_sig_bh", "rmats_sig", "sig_both",
                  "sig_both_same_direction", "same_direction_higher_in_test", "same_direction_lower_in_test"):
            v = int(r[k])
            per_comp[c][k] += v
            per_type[t][k] += v
            per_group_type[(group_of(c), t)][k] += v
        if r["dpsi_pearson"]:
            rho[t].append(float(r["dpsi_pearson"]))
    print("comparison\ttestable\tmatched\tsuppa_nominal\tsuppa_bh\trmats\tboth\tboth_same_dir\thigher\tlower")
    for c in [x["comparison_id"] for x in comps]:
        p = per_comp[c]
        print(f"{c}\t{p['suppa_testable']}\t{p['matched_testable']}\t{p['suppa_sig_nominal']}\t{p['suppa_sig_bh']}\t{p['rmats_sig']}"
              f"\t{p['sig_both']}\t{p['sig_both_same_direction']}\t{p['same_direction_higher_in_test']}\t{p['same_direction_lower_in_test']}")
    print("\n== per event type (sum over comparisons; Pearson range over comparisons) ==")
    print("type\ttestable\tmatched\tsuppa_nominal\tsuppa_bh\trmats\tboth\tboth_same_dir\thigher\tlower\tpearson_min\tpearson_max\tpearson_negative")
    for t in ("SE", "A3", "A5", "RI", "MX"):
        p = per_type[t]
        rs = rho[t]
        print(f"{t}\t{p['suppa_testable']}\t{p['matched_testable']}\t{p['suppa_sig_nominal']}\t{p['suppa_sig_bh']}\t{p['rmats_sig']}"
              f"\t{p['sig_both']}\t{p['sig_both_same_direction']}\t{p['same_direction_higher_in_test']}\t{p['same_direction_lower_in_test']}"
              f"\t{min(rs) if rs else ''}\t{max(rs) if rs else ''}\t{sum(1 for x in rs if x < 0)}")
    both = sum(p["sig_both"] for p in per_type.values())
    same = sum(p["sig_both_same_direction"] for p in per_type.values())
    print(f"events significant in both methods: {both}; same direction: {same} ({100 * same / both:.1f} percent)")
    print("\n== per group and type: high-confidence events (higher/lower in test) ==")
    for (g, t), p in sorted(per_group_type.items()):
        print(f"{g}\t{t}\t{p['sig_both_same_direction']}\t{p['same_direction_higher_in_test']}/{p['same_direction_lower_in_test']}")

    print("\n== distinct high-confidence events ==")
    ann = read_tsv(os.path.join(S, "events_high_confidence_annotated.tsv"))
    print(f"distinct events: {len(ann)}; by type: " + ", ".join(f"{t} {n}" for t, n in sorted(collections.Counter(r['event_type'] for r in ann).items())))
    multi = sum(1 for r in ann if "," in r["ref_gene"]); none = sum(1 for r in ann if not r["ref_gene"])
    print(f"single reference gene {len(ann) - multi - none}, several {multi}, none {none}; without eggNOG description {sum(1 for r in ann if not r['description'])}")
    genes = {g for r in ann for g in r["ref_gene"].split(",") if g}
    print(f"reference genes with at least one high-confidence event: {len(genes)}")
    rec = collections.Counter(int(r["n_comparisons"]) for r in ann)
    print("events by number of comparisons: " + ", ".join(f"{k}: {rec[k]}" for k in sorted(rec)))
    print(f"events in at least 4 comparisons: {sum(v for k, v in rec.items() if k >= 4)}; mixed direction: {sum(1 for r in ann if r['direction'] == 'mixed')}")
    for g in sorted({x for r in ann for x in r["groups"].split(",")}):
        rows = [r for r in ann if g in r["groups"].split(",")]
        gg = {x for r in rows for x in r["ref_gene"].split(",") if x}
        print(f"  {g}: {len(rows)} distinct events, {len(gg)} reference genes")
    cold = [r for r in ann if "cold" in r["groups"].split(",")]
    ncold = collections.Counter(sum(1 for c in r["comparisons"].split(",") if c.startswith("cold")) for r in cold)
    print("cold: events by number of cold comparisons: " + ", ".join(f"{k}: {ncold[k]}" for k in sorted(ncold) if k))
    core = [r for r in cold if sum(1 for c in r["comparisons"].split(",") if c.startswith("cold")) == 6]
    print(f"cold core (all six cold comparisons): {len(core)} events, {len({r['ref_gene'] for r in core})} reference genes, "
          f"mixed direction {sum(1 for r in core if r['direction'] == 'mixed')}")
    for r in sorted(core, key=lambda x: (x["direction"], x["ref_gene"])):
        print(f"  {r['event_type']}\t{r['direction']}\t{r['ref_gene']}\t{r['description'][:70]}\t{r['pfams']}")

    print("\n== consequences ==")
    cons = read_tsv(os.path.join(S, "events_high_confidence_consequences.tsv"))
    ctx = collections.Counter(r["cds_context"] for r in cons)
    print("location: " + ", ".join(f"{k} {v}" for k, v in sorted(ctx.items())))
    for t in ("SE", "A3", "A5", "RI", "MX"):
        rows = [r for r in cons if r["event_type"] == t]
        cds = [r for r in rows if r["cds_context"] == "CDS"]
        f0 = sum(1 for r in cds if r["frame"] == "0")
        print(f"  {t}: {len(rows)} events, CDS {len(cds)}, UTR {sum(1 for r in rows if r['cds_context'] == 'UTR')}, "
              f"frame preserved in CDS {f0} ({100 * f0 / len(cds):.1f} percent)" if cds else f"  {t}: {len(rows)} events")
    cds_no_mx = [r for r in cons if r["cds_context"] == "CDS" and r["event_type"] != "MX"]
    f0 = sum(1 for r in cds_no_mx if r["frame"] == "0")
    print(f"frame preserved in CDS, SE+A3+A5+RI: {f0} of {len(cds_no_mx)} ({100 * f0 / len(cds_no_mx):.1f} percent), "
          f"two-sided binomial test against one third p = {binom_two_sided(f0, len(cds_no_mx), 1 / 3):.2g}")
    ri = collections.Counter(r["ri_ptc"] for r in cons if r["event_type"] == "RI" and r["cds_context"] == "CDS")
    print(f"retained introns in CDS: stop codon in frame yes {ri['yes']}, no {ri['no']}, undetermined {ri['']}")

    nmd_path = os.path.join(S, "events_high_confidence_nmd.tsv")
    if os.path.isfile(nmd_path):
        print("\n== premature stop codons and NMD features (both forms on the annotated gene models) ==")
        nmd = read_tsv(nmd_path)
        cds_events = {r["suppa_event"] for r in cons if r["cds_context"] == "CDS"}
        for t in ("SE", "A3", "A5", "RI", "MX"):
            rows = [r for r in nmd if r["event_type"] == t and r["suppa_event"] in cds_events]
            if not rows:
                continue
            c = collections.Counter(r["ptc_form"] for r in rows)
            one = [r for r in rows if r["ptc_form"] in ("1", "2")]
            feat = sum(1 for r in one if r["nmd_form"])
            d = collections.Counter(r["ptc_form_in_test"] for r in one)
            print(f"  {t}: CDS events {len(rows)}, PTC in form 1 {c['1']}, form 2 {c['2']}, both {c['both']}, none {c['none']}, "
                  f"unresolved {c['unresolved']}; PTC form with NMD features {feat} of {len(one)}; "
                  f"PTC form higher in test {d['higher']}, lower {d['lower']}")
        one = [r for r in nmd if r["ptc_form"] in ("1", "2") and r["suppa_event"] in cds_events]
        print(f"  all CDS events with a PTC in one form: {len(one)}, with NMD features {sum(1 for r in one if r['nmd_form'])}; "
              f"junction rule {sum(1 for r in one if 'junction' in r['nmd_form'])}, long 3' UTR {sum(1 for r in one if 'long_utr3' in r['nmd_form'])}")
        ri = [r for r in nmd if r["event_type"] == "RI" and r["ri_ptc_consequences"] and r["ptc_form"] != "unresolved"]
        agree = sum(1 for r in ri if (r["form1_ptc"] == "yes") == (r["ri_ptc_consequences"] == "yes"))
        print(f"  RI cross-check with ri_ptc: {agree} of {len(ri)} agree")

    dom_path = os.path.join(S, "events_high_confidence_domains.tsv")
    if os.path.isfile(dom_path):
        print("\n== Pfam domains touched by the coding events ==")
        dom = read_tsv(dom_path)
        for label, keep in (("frame-preserving", lambda r: r["frame"] == "0"), ("frame-shifting", lambda r: r["frame"] != "0")):
            rows = [r for r in dom if keep(r)]
            placed = [r for r in rows if r["mrna"]]
            hit = [r for r in placed if r["n_domains"] != "0"]
            by_t = collections.Counter(r["event_type"] for r in hit)
            print(f"  {label}: {len(rows)} CDS events, placed {len(placed)}, overlapping a domain {len(hit)} ("
                  + ", ".join(f"{t} {by_t[t]}" for t in ("SE", "A3", "A5", "RI", "MX") if by_t[t]) + ")")
        names_of = lambda r: {x.split("(")[0] for x in r["domains_affected"].split(";") if x}   # one count per event
        nmd_path = os.path.join(S, "events_high_confidence_nmd.tsv")
        if os.path.isfile(nmd_path):
            ptc = {r["suppa_event"]: r["ptc_form"] for r in read_tsv(nmd_path)}
            f0 = [r for r in dom if r["frame"] == "0" and r["mrna"] and ptc.get(r["suppa_event"], "unresolved") != "unresolved"]
            groups = (("region with an in-frame stop codon (form 1 truncated)", [r for r in f0 if ptc[r["suppa_event"]] == "1"]),
                      ("stop codon only in the excluding form", [r for r in f0 if ptc[r["suppa_event"]] == "2"]),
                      ("complete reading frame in both forms (in-frame insertion or deletion)", [r for r in f0 if ptc[r["suppa_event"]] == "none"]))
            print(f"  frame-preserving events placed and reconstructed: {len(f0)}")
            for label, rows in groups:
                hit = [r for r in rows if r["n_domains"] != "0"]
                aa = sorted(float(r["region_aa"]) for r in rows)
                names = collections.Counter(n for r in hit for n in names_of(r))
                print(f"    {label}: {len(rows)} events (median region {aa[len(aa) // 2] if aa else 0:.0f} residues), "
                      f"at or within a domain {len(hit)}; events per domain: " + ", ".join(f"{k} {v}" for k, v in names.most_common(10)))
            print("  frame-preserving events at or within a domain: event type, gene, stop codon form, aa interval, domains, groups, description")
            for r in sorted((r for r in f0 if r["n_domains"] != "0"), key=lambda x: (ptc[x["suppa_event"]], x["event_type"], x["gene"])):
                print(f"    {r['event_type']}\t{r['gene']}\t{ptc[r['suppa_event']]}\t{r['aa_start']}-{r['aa_end']}\t{r['domains_affected']}\t{r['groups']}\t{r['description'][:60]}")

    aba_path = os.path.join(S, "aba_pathway_enrichment.tsv")
    if os.path.isfile(aba_path):
        print("\n== abscisic acid pathway genes ==")
        for r in read_tsv(aba_path):
            print(f"  {r['set']}: {r['pathway_genes_with_events']} of {r['genes_with_events']} genes with events "
                  f"(pathway genes in background {r['pathway_genes_in_background']} of {r['background_genes']}, "
                  f"expected {r['expected']}, p {float(r['p_hypergeometric']):.3g})")
        for r in read_tsv(os.path.join(S, "aba_pathway_events.tsv")):
            print(f"    {r['gene']}\t{r['preferred_name']}\t{r['component']}\t{r['event_type']}\t{r['direction']}\t{r['groups']}\t{r['description'][:60]}")

    if a.expression and os.path.isfile(os.path.join(a.expression, "de_summary.tsv")):
        print("\n== differential expression and its overlap with differential splicing ==")
        de = read_tsv(os.path.join(a.expression, "de_summary.tsv"))
        print("  comparison, tested, DE (up, down), DAS, population, DE and DAS, expected, p")
        for r in de:
            print(f"    {r['comparison']}\t{r['genes_tested']}\t{r['de_genes']} ({r['de_up']}, {r['de_down']})\t{r['das_genes']}\t"
                  f"{r['population']}\t{r['de_and_das']}\t{r['expected_overlap']}\t{float(r['p_hypergeometric']):.3g}")
        tot_das = sum(int(r["das_in_population"]) for r in de)
        tot_both = sum(int(r["de_and_das"]) for r in de)
        tot_exp = sum(float(r["expected_overlap"]) for r in de)
        sig = sum(1 for r in de if float(r["p_hypergeometric"]) < 0.05)
        print(f"  over all comparisons: DAS genes in the population {tot_das}, of them DE {tot_both} "
              f"({100 * tot_both / tot_das:.1f} percent), expected {tot_exp:.1f}; comparisons with p < 0.05: {sig} of {len(de)}")
        reg_path = os.path.join(a.expression, "de_regulators.tsv")
        if os.path.isfile(reg_path):
            reg = read_tsv(reg_path)
            de_reg = [r for r in reg if r["de"] == "yes"]
            print(f"  splicing-regulator genes: {len({r['gene'] for r in reg})} genes tested, DE in at least one comparison "
                  f"{len({r['gene'] for r in de_reg})}; DE rows by family: "
                  + ", ".join(f"{k} {v}" for k, v in sorted(collections.Counter(r["family"] for r in de_reg).items())))
            for r in sorted(de_reg, key=lambda x: (x["family"], x["gene"], x["comparison"])):
                print(f"    {r['family']}\t{r['gene']}\t{r['ath_gene']}\t{r['comparison']}\tlog2FC {r['log2fc']}\tpadj {r['padj']}\tDAS {r['das']}")

    cross_path = os.path.join(S, "cross_species_cold_summary.tsv")
    if os.path.isfile(cross_path):
        print("\n== cold splicing genes against the Arabidopsis cold DAS list (Calixto et al. 2018) ==")
        for r in read_tsv(cross_path):
            print(f"  {r['set']}: {r['counterpart_in_list']} of {r['genes']} linked genes have a counterpart in the list "
                  f"(share {r['share']}, expected {r['expected']}, enrichment {r['enrichment']}, p {r['p_hypergeometric']})")

    print("\n== splicing-regulator families ==")
    D = a.domains
    fam = read_tsv(os.path.join(D, "splicing_regulators_families.tsv"))
    for f in sorted({r["family"] for r in fam}):
        rows = [r for r in fam if r["family"] == f]
        print(f"  {f}: {len(rows)} genes, reciprocal best {sum(1 for r in rows if r['reciprocal_best'] == 'yes')}, "
              f"with event {sum(1 for r in rows if r['das_events'])}")
    sr = [r for r in fam if r["family"] == "SR"]
    print("  SR subfamilies: " + ", ".join(f"{k} {v}" for k, v in sorted(collections.Counter(r["subfamily"] for r in sr).items())))
    print("  SR genes: gene, subfamily, ath gene, rule, pident, reciprocal, events")
    for r in sorted(sr, key=lambda x: (x["subfamily"], x["gene"])):
        print(f"    {r['gene']}\t{r['subfamily']}\t{r['ath_gene']}\t{r['sr_by_sequence']}\t{r['best_pident']}\t{r['reciprocal_best']}\t{r['das_events']}")
    enr = read_tsv(os.path.join(D, "splicing_regulators_enrichment.tsv"))
    print("  enrichment: family, members, with event, expected, p")
    for r in enr:
        print(f"    {r['family']}\t{r['members_in_background']}\t{r['with_event']}\t{r['expected']}\t{float(r['p_hypergeometric']):.3g}")
    if enr:
        print(f"  background {enr[0]['background']} genes, {enr[0]['background_with_event']} with event")

    print("\n== GO enrichment ==")
    go = read_tsv(os.path.join(S, "go_enrichment.tsv"))
    for s in sorted({r["set"] for r in go}):
        rows = [r for r in go if r["set"] == s]
        sig = [r for r in rows if float(r["fdr"]) < a.fdr]
        print(f"  {s}: foreground {rows[0]['foreground_genes']} genes, terms tested {len(rows)}, FDR<{a.fdr}: {len(sig)}")
        for r in [x for x in sig if x["namespace"] != "cellular_component"][:a.top]:
            print(f"    {r['go']}\t{r['name']}\t{r['term_foreground']}/{r['term_background']}\tfold {r['fold']}\tFDR {float(r['fdr']):.2g}")


if __name__ == "__main__":
    main()
