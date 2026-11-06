"""Figures of the manuscript from the result tables in supplementary/tables and the metadata.

Figure 1  data set: runs per condition category, STAR unique mapping, transcript filtering, event types
Figure 2  high-confidence events per comparison by event type, and their direction
Figure 3  the cold core: change in PSI of the events found in all six cold comparisons
Figure 4  consequences (location and frame), family enrichment and GO terms
Figure 5 (the SR protein tree) is drawn by make_tree_figure.py.

Usage (from the repository root):  python scripts/figures/make_figures.py
"""
import collections
import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T = os.path.join(ROOT, "supplementary", "tables")
META = os.path.join(ROOT, "metadata")
OUT = os.path.join(ROOT, "manuscript", "figures")
os.makedirs(OUT, exist_ok=True)
TYPES = ["SE", "A3", "A5", "RI", "MX"]
TYPE_NAMES = {"SE": "skipped exon", "A3": "alternative 3' site", "A5": "alternative 5' site", "RI": "retained intron", "MX": "mutually exclusive"}
GROUP_NAMES = {"cold": "cold", "imbibition": "seed dormancy", "tablebeet": "taproot development", "aphano": "Aphanomyces",
               "transgenic": "transgenic", "bctv": "BCTV", "sclerotinia": "Sclerotinia"}
FIRST = ["cold4C_GT1_leaf", "cold4C_GT3_leaf", "cold0C_GT1_leaf", "cold4C_GT1_taproot", "cold4C_GT2_taproot", "cold4C_GT3_taproot",
         "imbibition_10C_60h_vs_dry", "imbibition_10C_82h_vs_60h", "imbibition_5C_vs_10C_60h", "imbibition_5C_160h_vs_60h",
         "imbibition_5C_272h_vs_60h", "imbibition_5C_432h_vs_60h", "imbibition_5C_272h_vs_10C_82h"]
GROUP_ORDER = ["cold", "imbibition", "tablebeet", "aphano", "transgenic", "bctv", "sclerotinia"]
COLORS = {"SE": "#4477AA", "A3": "#66CCEE", "A5": "#228833", "RI": "#CCBB44", "MX": "#EE6677"}
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})



def flatten_tiff(path):
    """RGBA TIFF written by Matplotlib to RGB on white, LZW, same resolution."""
    from PIL import Image
    im = Image.open(path)
    dpi = im.info.get("dpi", (600, 600))
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im = bg
    im.save(path, compression="tiff_lzw", dpi=dpi)

def read(path, comment=True):
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if not (comment and l.startswith("#"))]
    return list(csv.DictReader(lines, delimiter="\t"))


def group_of(cid):
    return re.sub(r"\d+C$", "", cid.split("_")[0])


def save(fig, name):
    for ext in ("png", "svg", "tif"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600 if ext == "tif" else 300 if ext == "png" else None, bbox_inches="tight", metadata=({"Software": None} if ext == "png" else {"Creator": None, "Date": None} if ext == "svg" else None), **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
        if ext == "tif":
            flatten_tiff(os.path.join(OUT, f"{name}.tif"))
    plt.close(fig)
    print("wrote", name)


def figure1():
    samples = read(os.path.join(META, "samples_all.tsv"))
    star = read(os.path.join(T, "star_summary.tsv"))
    struct = read(os.path.join(T, "structural_report.tsv"))
    events = read(os.path.join(T, "suppa_event_counts.tsv"))
    fig, ax = plt.subplots(2, 2, figsize=(6.8, 5.3))

    cat = collections.Counter(r["category"] for r in samples if r["use_das"] != "excluded")
    labels = {"biotic_oomycete": "Aphanomyces", "postharvest": "postharvest storage", "abiotic_cold": "cold", "biotic_virus": "beet curly top virus",
              "development": "development", "genotype_panel": "genotype panel", "abiotic_cold_dormancy": "seed dormancy", "abiotic_multi": "salt, heat, light",
              "abiotic_germination": "germination", "biotic_virus_transgenic": "transgenic lines", "biotic_fungal": "Sclerotinia", "misc": "other"}
    items = sorted(cat.items(), key=lambda x: -x[1])
    a = ax[0, 0]
    a.barh([labels.get(k, k) for k, _ in items][::-1], [v for _, v in items][::-1], color="#777777")
    a.set_xlabel("RNA-seq runs")
    a.set_title("a", loc="left", fontweight="bold")

    a = ax[0, 1]
    uniq = [float(r["uniq_pct"]) for r in star]
    a.hist(uniq, bins=np.arange(0, 101, 5), color="#4477AA")
    a.axvline(60, color="black", linestyle="--", linewidth=0.8)
    a.set_xlabel("uniquely mapped reads (%)"); a.set_ylabel("runs")
    a.set_title("b", loc="left", fontweight="bold")

    a = ax[1, 0]
    kept = {r["class_code"]: int(r["kept"]) for r in struct}
    dropped = sum(int(r["dropped"]) for r in struct)
    stages = ["merged" + chr(10) + "assembly", "structural" + chr(10) + "filter", "expression" + chr(10) + "filter"]
    values = [sum(kept.values()) + dropped, sum(kept.values()), 73743]
    a.bar(stages, values, color=["#BBBBBB", "#999999", "#4477AA"])
    for i, v in enumerate(values):
        a.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=7)
    a.set_ylabel("transcripts"); a.set_ylim(0, max(values) * 1.15)
    a.set_title("c", loc="left", fontweight="bold")
    a.tick_params(axis="x", labelsize=7)

    a = ax[1, 1]
    ev = sorted(events, key=lambda r: -int(r["events"]))
    names = {"A3": "alt. 3' site", "A5": "alt. 5' site", "SE": "skipped exon", "RI": "retained intron", "AF": "alt. first exon", "AL": "alt. last exon", "MX": "mutually excl."}
    a.barh([names[r["event_type"]] for r in ev][::-1], [int(r["events"]) for r in ev][::-1], color="#4477AA")
    a.set_xlabel("events in the expressed transcriptome")
    a.set_title("d", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "figure1_dataset")


def pretty_label(c):
    """Readable comparison label: temperatures with the degree sign, units separated, isolate and line names restored."""
    s = re.sub(r"^(aphano|bctv|imbibition|tablebeet|transgenic)_", "", c).replace("_vs_", " vs ").replace("_", " ")
    s = re.sub(r"cold(\d+)C", lambda m: m.group(1) + " " + chr(176) + "C", s)
    s = re.sub(r"\b(\d+)C\b", lambda m: m.group(1) + " " + chr(176) + "C", s)
    s = re.sub(r"\b(\d+)(h|hpi|dpi)\b", r"\1 \2", s)
    s = s.replace("01 7 6", "01.7.6").replace("EL 10", "EL10").replace("FC709 2", "FC709-2").replace("sclerotinia", "Sclerotinia")
    return s


def figure2():
    conc = read(os.path.join(T, "event_concordance_summary.tsv"))
    comps = [r["comparison_id"] for r in read(os.path.join(META, "comparisons_all.tsv"))]
    same = collections.defaultdict(dict); up = collections.defaultdict(int); down = collections.defaultdict(int)
    for r in conc:
        same[r["comparison"]][r["event_type"]] = int(r["sig_both_same_direction"])
        up[r["comparison"]] += int(r["same_direction_higher_in_test"]); down[r["comparison"]] += int(r["same_direction_lower_in_test"])
    order = sorted(comps, key=lambda c: (GROUP_ORDER.index(group_of(c)), FIRST.index(c) if c in FIRST else 99, c))
    short = [pretty_label(c) for c in order]
    fig, ax = plt.subplots(2, 1, figsize=(6.8, 6.1), sharex=True, gridspec_kw={"height_ratios": [3, 1.6]})
    x = np.arange(len(order)); bottom = np.zeros(len(order))
    for t in TYPES:
        vals = np.array([same[c].get(t, 0) for c in order])
        ax[0].bar(x, vals, bottom=bottom, color=COLORS[t], label=TYPE_NAMES[t], width=0.8)
        bottom += vals
    ax[0].set_ylabel("high-confidence events"); ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a", loc="left", fontweight="bold")
    # group brackets
    pos = 0
    for g in GROUP_ORDER:
        n = sum(1 for c in order if group_of(c) == g)
        if n:
            ax[0].annotate(GROUP_NAMES[g].replace("taproot development", "taproot" + chr(10) + "development").replace("seed dormancy", "seed" + chr(10) + "dormancy"), xy=((pos + pos + n - 1) / 2, bottom.max() * 1.02), ha="center", va="bottom", fontsize=7, annotation_clip=False)
            pos += n
    ax[1].bar(x, [up[c] for c in order], color="#EE6677", width=0.8, label="higher inclusion in test")
    ax[1].bar(x, [-down[c] for c in order], color="#4477AA", width=0.8, label="lower inclusion in test")
    ax[1].axhline(0, color="black", linewidth=0.6); ax[1].legend(frameon=False, fontsize=7); ax[1].set_ylabel("events")
    ax[1].set_title("b", loc="left", fontweight="bold")
    ax[1].set_xticks(x); ax[1].set_xticklabels(short, rotation=90, fontsize=7)
    fig.tight_layout()
    save(fig, "figure2_events_per_comparison")


def figure3():
    hc = read(os.path.join(T, "event_concordance_high_confidence.tsv"))
    ann = {r["suppa_event"]: r for r in read(os.path.join(T, "events_high_confidence_annotated.tsv"))}
    cold = ["cold4C_GT1_leaf", "cold4C_GT3_leaf", "cold0C_GT1_leaf", "cold4C_GT1_taproot", "cold4C_GT2_taproot", "cold4C_GT3_taproot"]
    per = collections.defaultdict(dict)
    for r in hc:
        if r["comparison"] in cold:
            per[r["suppa_event"]][r["comparison"]] = (float(r["suppa_dpsi"]) + float(r["rmats_dpsi"])) / 2
    core = [e for e, d in per.items() if len(d) == 6]
    rows = []
    for e in core:
        a = ann[e]
        desc = a["description"].split(" | ")[0][:80] or a["ref_gene"].split(",")[0] or a["merged_gene"]
        rows.append((np.mean([per[e][c] for c in cold]), f"{a['event_type']}  {desc}", [per[e][c] for c in cold]))
    rows.sort(key=lambda r: r[0])
    mat = np.array([r[2] for r in rows])
    fig, ax = plt.subplots(figsize=(6.8, 0.22 * len(rows) + 1.2))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[1] for r in rows], fontsize=7)
    ax.set_xticks(range(6)); ax.set_xticklabels(["4 °C GT1 leaf", "4 °C GT3 leaf", "0 °C GT1 leaf", "4 °C GT1 taproot", "4 °C GT2 taproot", "4 °C GT3 taproot"], rotation=45, ha="right", fontsize=7)
    for i in range(len(rows)):
        for j in range(6):
            ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=6.5, color="white" if abs(mat[i, j]) > 0.4 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("change in PSI, cold minus control (mean of both methods)", fontsize=7)
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    save(fig, "figure3_cold_core")


def figure4():
    cons = read(os.path.join(T, "events_high_confidence_consequences.tsv"))
    enr = read(os.path.join(T, "splicing_regulators_enrichment.tsv"))
    go = read(os.path.join(T, "go_enrichment.tsv"))
    fig = plt.figure(figsize=(6.8, 7.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.5])

    a = fig.add_subplot(gs[0, 0])
    loc = {t: collections.Counter(r["cds_context"] for r in cons if r["event_type"] == t) for t in TYPES}
    bottom = np.zeros(len(TYPES))
    for ctx, col in (("CDS", "#4477AA"), ("UTR", "#CCBB44"), ("none", "#BBBBBB"), ("noncoding", "#999999")):
        vals = np.array([loc[t].get(ctx, 0) for t in TYPES])
        if vals.sum():
            a.bar(TYPES, vals, bottom=bottom, color=col, label="outside gene models" if ctx == "none" else ctx); bottom += vals
    a.set_ylabel("distinct events"); a.set_ylim(0, 820); a.legend(frameon=False, fontsize=7, loc="upper right")
    a.set_title("a", loc="left", fontweight="bold")

    a = fig.add_subplot(gs[0, 1])
    fr = []
    for t in TYPES[:4]:
        cds = [r for r in cons if r["event_type"] == t and r["cds_context"] == "CDS"]
        fr.append(100 * sum(1 for r in cds if r["frame"] == "0") / len(cds))
    ri = collections.Counter(r["ri_ptc"] for r in cons if r["event_type"] == "RI" and r["cds_context"] == "CDS")
    a.bar(TYPES[:4], fr, color="#4477AA"); a.axhline(100 / 3, color="black", linestyle="--", linewidth=0.8)
    a.set_ylabel("coding events with length multiple of 3 (%)"); a.set_ylim(0, 60)
    a.text(3.6, 52, f"retained introns in CDS:\n{ri['yes']} with in-frame stop,\n{ri['no']} without,\n{ri['']} undetermined", fontsize=7, ha="right", va="top")
    a.set_title("b", loc="left", fontweight="bold")

    a = fig.add_subplot(gs[1, 0])
    fams = [r for r in enr if r["family"] not in ("all families",)]
    fams.sort(key=lambda r: float(r["p_hypergeometric"]))
    y = np.arange(len(fams))
    a.barh(y, [float(r["expected"]) for r in fams], color="#BBBBBB", height=0.4, label="expected")
    a.barh(y + 0.4, [int(r["with_event"]) for r in fams], color="#EE6677", height=0.4, label="observed")
    for i, r in enumerate(fams):
        a.text(max(float(r["expected"]), int(r["with_event"])) + 0.3, i + 0.2, f"{int(r['with_event'])}/{r['members_in_background']}  p={float(r['p_hypergeometric']):.2g}", fontsize=6, va="center")
    a.set_yticks(y + 0.2); a.set_yticklabels([r["family"] for r in fams], fontsize=7); a.invert_yaxis()
    a.set_xlabel("genes with high-confidence events"); a.legend(frameon=False, fontsize=7, loc="lower right"); a.set_xlim(0, 22)
    a.set_title("c", loc="left", fontweight="bold")

    a = fig.add_subplot(gs[1, 1])
    sets = [("all", "all"), ("cold", "cold"), ("imbibition", "seed dormancy"), ("tablebeet", "taproot development")]
    terms = ["GO:0008380", "GO:0000398", "GO:0043484", "GO:0006397", "GO:0003723", "GO:0009081", "GO:0006551"]
    names = {}
    grid = {}
    for r in go:
        if r["go"] in terms and r["set"] in dict(sets):
            grid[(r["set"], r["go"])] = (float(r["fold"]), float(r["fdr"])); names[r["go"]] = r["name"]
    for j, (s, lab) in enumerate(sets):
        for i, t in enumerate(terms):
            if (s, t) in grid:
                fold, fdr = grid[(s, t)]
                a.scatter(j, i, s=12 + 40 * fold, color="#4477AA" if fdr < 0.05 else "#CCCCCC", edgecolor="black", linewidth=0.3)
    a.set_xticks(range(len(sets))); a.set_xticklabels([l for _, l in sets], fontsize=7, rotation=30, ha="right")
    a.set_yticks(range(len(terms))); a.set_yticklabels([names.get(t, t)[:48] for t in terms], fontsize=7); a.invert_yaxis()
    a.set_xlim(-0.6, len(sets) - 0.4)
    a.set_title("d", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "figure4_consequences_regulators")


if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4()
