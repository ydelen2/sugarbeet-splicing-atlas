"""Supplementary figures from the packaged result files in supplementary/pkg.

Figure S1  length and GC content of retained introns and skipped exons: annotated regions of BeetSet-3,
           all events of the event universe, and the high-confidence (consensus) events

Usage (from the repository root):  python scripts/figures/make_supplementary_figures.py
"""
import collections
import csv
import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PKG = os.path.join(ROOT, "supplementary", "pkg")
OUT = os.path.join(ROOT, "manuscript", "figures")
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
SETS = [("annotated", "annotated (BeetSet-3)", "#BBBBBB"), ("tested", "all events", "#88CCEE"), ("consensus", "high-confidence events", "#CC3311")]



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

def figure_s1():
    path = os.path.join(PKG, "features", "region_features_regions.tsv.gz")
    if not os.path.isfile(path):
        print("missing", path); return
    vals = collections.defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vals[(r["type"], r["set"])].append((int(r["length"]), float(r["gc"])))
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.4))
    for j, (etype, title) in enumerate((("RI", "retained introns"), ("SE", "skipped exons"))):
        for i, (feature, label, idx, log) in enumerate((("length", "length (nt)", 0, True), ("GC", "GC content", 1, False))):
            a = ax[i, j]
            data = [[v[idx] for v in vals[(etype, s)]] for s, _, _ in SETS]
            if log:
                data = [np.log10(d) for d in data]
            parts = a.violinplot(data, showmedians=True, widths=0.8)
            for body, (_, _, col) in zip(parts["bodies"], SETS):
                body.set_facecolor(col); body.set_alpha(0.8)
            a.set_xticks([1, 2, 3]); a.set_xticklabels([f"{lab}\n(n = {len(d):,})" for (_, lab, _), d in zip(SETS, data)], fontsize=7)
            a.set_ylabel("log10 length (nt)" if log else "GC fraction")
            a.set_title("abcd"[i * 2 + j], loc="left", fontweight="bold")
    fig.tight_layout()
    for ext in ("png", "svg", "tif"):
        fig.savefig(os.path.join(OUT, f"figure_s1_region_features.{ext}"), dpi=600 if ext == "tif" else 300 if ext == "png" else None, bbox_inches="tight", metadata=({"Software": None} if ext == "png" else {"Creator": None, "Date": None} if ext == "svg" else None), **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
        if ext == "tif":
            flatten_tiff(os.path.join(OUT, "figure_s1_region_features.tif"))
    print("wrote figure_s1_region_features")


if __name__ == "__main__":
    figure_s1()
