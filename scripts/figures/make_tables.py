"""Tables of the manuscript, written as Markdown from the metadata and the result tables.

Table 1  data sets: BioProject, condition category, runs analyzed, comparisons built
Table 2  high-confidence events per condition group and event type, distinct events and genes
Table 3  the cold core: events found in all six cold comparisons
Table 4  splicing-regulator families and the enrichment of high-confidence events

Usage (from the repository root):  python scripts/figures/make_tables.py
"""
import collections
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T = os.path.join(ROOT, "supplementary", "tables")
META = os.path.join(ROOT, "metadata")
OUT = os.path.join(ROOT, "manuscript", "tables")
os.makedirs(OUT, exist_ok=True)
TYPES = ["SE", "A3", "A5", "RI", "MX"]
GROUP_ORDER = ["cold", "imbibition", "tablebeet", "aphano", "transgenic", "bctv", "sclerotinia"]
GROUP_NAMES = {"cold": "Cold (4 °C and 0 °C)", "imbibition": "Cold-induced seed dormancy", "tablebeet": "Taproot development (table beet)",
               "aphano": "Aphanomyces cochlioides infection", "transgenic": "Rhizomania-resistant transgenic lines", "bctv": "Beet curly top virus infection",
               "sclerotinia": "Sclerotinia sclerotiorum infection"}
CATEGORY = {"biotic_oomycete": "Aphanomyces cochlioides infection", "postharvest": "postharvest storage", "abiotic_cold": "cold",
            "biotic_virus": "beet curly top virus", "development": "development", "genotype_panel": "genotype panel",
            "abiotic_cold_dormancy": "cold-induced seed dormancy", "abiotic_multi": "salt, heat and high light",
            "abiotic_germination": "seed germination and root development", "biotic_virus_transgenic": "transgenic lines",
            "biotic_fungal": "Sclerotinia sclerotiorum infection", "misc": "single runs"}


def read(path, comment=True):
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if not (comment and l.startswith("#"))]
    return list(csv.DictReader(lines, delimiter="\t"))


def group_of(cid):
    return re.sub(r"\d+C$", "", cid.split("_")[0])


def md(rows, header, path, caption):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(caption + "\n\n| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(str(x) for x in r) + " |\n")
    print("wrote", os.path.relpath(path, ROOT))


def table1():
    samples = [r for r in read(os.path.join(META, "samples_all.tsv")) if r["use_das"] != "excluded"]
    comps = read(os.path.join(META, "comparisons_all.tsv"))
    eligible = [r for r in samples if r["use_das"] == "yes"]
    used = set()
    for c in comps:
        sel = [r for r in eligible if r["study"] == c["study"] and r["condition"] in (c["reference"], c["test"])]
        if c["subset"] not in ("", "-"):
            for kv in c["subset"].split(";"):
                k, v = kv.split("=", 1)
                sel = [r for r in sel if r[k] == v]
        used |= {r["run"] for r in sel}
    per_study = collections.defaultdict(lambda: {"runs": 0, "das": 0, "used": 0, "cat": collections.Counter(), "tissues": set(), "comps": 0})
    for r in samples:
        d = per_study[r["study"]]
        d["runs"] += 1; d["das"] += r["use_das"] == "yes"; d["used"] += r["run"] in used; d["cat"][CATEGORY.get(r["category"], r["category"])] += 1; d["tissues"].add(r["tissue"].replace("+", " and "))
    for c in comps:
        per_study[c["study"]]["comps"] += 1
    rows = []
    for st, d in sorted(per_study.items(), key=lambda x: -x[1]["runs"]):
        # a study whose runs fall into several categories lists each with its run count
        cat = next(iter(d["cat"])) if len(d["cat"]) == 1 else "; ".join(f"{k} ({v})" for k, v in d["cat"].most_common())
        rows.append([st, cat, ", ".join(sorted(t for t in d["tissues"] if t)), d["runs"], d["das"], d["used"], d["comps"]])
    rows.append(["Total", "", "", sum(d["runs"] for d in per_study.values()), sum(d["das"] for d in per_study.values()),
                 sum(d["used"] for d in per_study.values()), len(comps)])
    md(rows, ["BioProject", "Condition category", "Tissues", "Runs analyzed", "Runs eligible for comparisons", "Runs in comparisons", "Comparisons"],
       os.path.join(OUT, "table1_datasets.md"),
       "**Table 1** Public RNA-seq data sets used in this study. Runs eligible for comparisons are those of studies with replicated designs "
       "whose library strandedness is consistent within the study (flagged in the sample table of the repository); the others, "
       "including studies without replicated groups that differ in one factor, contributed to junction and transcript "
       "discovery only. Eligible runs that entered no comparison lacked an isolate-matched control (35 Aphanomyces runs), "
       "a second replicate in their group, or belonged to the GT2 cold leaf group used only in the organ comparisons. The 96 runs of the postharvest storage study have no condition labels in the archive and were used "
       "for transcript discovery only.")


def table2():
    conc = read(os.path.join(T, "event_concordance_summary.tsv"))
    ann = read(os.path.join(T, "events_high_confidence_annotated.tsv"))
    comps = read(os.path.join(META, "comparisons_all.tsv"))
    samples = read(os.path.join(META, "samples_all.tsv"))
    eligible = [r for r in samples if r["use_das"] == "yes"]

    def runs_of(c):
        sel = [r for r in eligible if r["study"] == c["study"] and r["condition"] in (c["reference"], c["test"])]
        if c["subset"] not in ("", "-"):
            for kv in c["subset"].split(";"):
                k, v = kv.split("=", 1)
                sel = [r for r in sel if r[k] == v]
        return {r["run"] for r in sel}

    rows = []
    tot = collections.Counter()
    all_runs = set()
    for g in GROUP_ORDER:
        cs = [c for c in comps if group_of(c["comparison_id"]) == g]
        runs = set()
        for c in cs:
            runs |= runs_of(c)
        all_runs |= runs
        by_type = collections.Counter()
        for r in conc:
            if group_of(r["comparison"]) == g:
                by_type[r["event_type"]] += int(r["sig_both_same_direction"])
        events = [a for a in ann if g in a["groups"].split(",")]
        genes = {x for a in events for x in a["ref_gene"].split(",") if x}
        rows.append([GROUP_NAMES[g], len(cs), len(runs)] + [by_type[t] for t in TYPES] + [sum(by_type.values()), len(events), len(genes)])
        tot["comps"] += len(cs); tot["runs"] += len(runs)
        for t in TYPES:
            tot[t] += by_type[t]
    all_genes = {x for a in ann for x in a["ref_gene"].split(",") if x}
    rows.append(["Total", tot["comps"], len(all_runs)] + [tot[t] for t in TYPES] + [sum(tot[t] for t in TYPES), len(ann), len(all_genes)])
    md(rows, ["Condition group", "Comparisons", "Runs", "SE", "A3", "A5", "RI", "MX", "Event-comparison pairs", "Distinct events", "Genes"],
       os.path.join(OUT, "table2_events_per_group.md"),
       "**Table 2** High-confidence differential splicing events (called by SUPPA2 and rMATS with the same direction) per condition group. "
       "SE skipped exon, A3 and A5 alternative 3' and 5' splice site, RI retained intron, MX mutually exclusive exons. Distinct events and "
       "genes are counted once per group; the totals count each event and gene once.")


def table3():
    hc = read(os.path.join(T, "event_concordance_high_confidence.tsv"))
    ann = {r["suppa_event"]: r for r in read(os.path.join(T, "events_high_confidence_annotated.tsv"))}
    cold = ["cold4C_GT1_leaf", "cold4C_GT3_leaf", "cold0C_GT1_leaf", "cold4C_GT1_taproot", "cold4C_GT2_taproot", "cold4C_GT3_taproot"]
    per = collections.defaultdict(dict)
    for r in hc:
        if r["comparison"] in cold:
            per[r["suppa_event"]][r["comparison"]] = ((float(r["suppa_dpsi"]) + float(r["rmats_dpsi"])) / 2)
    rows = []
    for e, d in per.items():
        if len(d) < 6:
            continue
        a = ann[e]
        vals = [d[c] for c in cold]
        chrom = re.sub(r"^0_Bv(chr\d+)_bp\d+$", r"\1", e.split(";")[1].split(":")[1])
        coords = ":".join(e.split(";")[1].split(":")[2:-1])
        rows.append([a["event_type"], a["ref_gene"] or "none", (a["description"].split(" | ")[0] or "no description"), a["pfams"].split(" | ")[0],
                     f"{chrom}:{coords}", f"{min(vals):+.2f} to {max(vals):+.2f}", a["direction"], a["n_comparisons"]])
    rows.sort(key=lambda r: (r[6], r[1]))
    md(rows, ["Type", "Gene", "eggNOG description", "Pfam", "Coordinates (RefBeet-3.0)", "Change in PSI, cold minus control", "Direction", "Comparisons in the atlas"],
       os.path.join(OUT, "table3_cold_core.md"),
       "**Table 3** Events with high-confidence differential splicing in all six cold comparisons (two studies, three genotypes, leaf and taproot). "
       "The change in PSI is the mean of the SUPPA2 and rMATS estimates; the range covers the six comparisons. Direction is given over all "
       "comparisons of the atlas in which the event was called, and the last column counts those comparisons.")


def table4():
    fam = read(os.path.join(T, "splicing_regulators_families.tsv"))
    enr = {r["family"]: r for r in read(os.path.join(T, "splicing_regulators_enrichment.tsv"))}
    order = ["SR", "SR-like", "U2AF", "SF1", "hnRNP-like", "PRP", "snRNP", "other splicing factor"]
    rows = []
    for f in order:
        rs = [r for r in fam if r["family"] == f]
        subs = collections.Counter(r["subfamily"] for r in rs)
        e = enr.get(f, {})
        rows.append([f, len(rs), ", ".join(f"{k} {v}" for k, v in sorted(subs.items()) if k != "-"), sum(1 for r in rs if r["reciprocal_best"] == "yes"),
                     sum(1 for r in rs if r["das_events"]), e.get("members_in_background", ""), e.get("with_event", ""), e.get("expected", ""),
                     f"{float(e['p_hypergeometric']):.2g}" if e else ""])
    e = enr["all families"]
    rows.append(["All families", sum(1 for r in fam), "", sum(1 for r in fam if r["reciprocal_best"] == "yes"), sum(1 for r in fam if r["das_events"]),
                 e["members_in_background"], e["with_event"], e["expected"], f"{float(e['p_hypergeometric']):.2g}"])
    md(rows, ["Family", "Genes", "Subfamilies", "Reciprocal best hits", "Genes with events", "Members in test background", "With events", "Expected", "p (hypergeometric)"],
       os.path.join(OUT, "table4_families.md"),
       "**Table 4** Splicing-regulator gene families of sugar beet defined by best hits to a curated list of Arabidopsis thaliana proteins "
       "(given in the repository), and the enrichment of high-confidence differential splicing events among their members "
       f"(the background is the {int(e['background']):,} genes with a testable event, {int(e['background_with_event']):,} of them with a high-confidence event). The 14 SR genes with events comprise 12 of the 14 genes with complete domain architecture "
       "(hypergeometric p = 3.6e-9) and the two fragmentary RS31-like models.")


if __name__ == "__main__":
    table1(); table2(); table3(); table4()
