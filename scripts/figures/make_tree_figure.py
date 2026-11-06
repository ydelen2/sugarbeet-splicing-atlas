"""Figure 5: the SR protein tree, drawn from the IQ-TREE Newick file with matplotlib.

The tree is rooted at the midpoint of the longest leaf-to-leaf path and drawn as a rectangular
phylogram. Leaves are colored by species, labeled with the species prefix and the gene, and the
subfamily of each clade (from the sequence table) is written at the right margin. Ultrafast bootstrap
values of 95 or more are marked with a filled circle at the node.

Usage (from the repository root):  python scripts/figures/make_tree_figure.py
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts", "analysis"))
from tree_summary import parse_newick, leaves, Node  # noqa: E402

T = os.path.join(ROOT, "supplementary", "tables")
OUT = os.path.join(ROOT, "manuscript", "figures")
SPECIES_COLOR = {"arabidopsis": "#000000", "beet": "#CC3311", "quinoa": "#0077BB", "spinach": "#009988", "amaranth": "#EE7733"}
SPECIES_NAME = {"arabidopsis": "Arabidopsis thaliana", "beet": "Beta vulgaris", "quinoa": "Chenopodium quinoa",
                "spinach": "Spinacia oleracea", "amaranth": "Amaranthus hypochondriacus"}
ABBREV = {"arabidopsis": "At", "beet": "Bv", "quinoa": "Cq", "spinach": "So", "amaranth": "Ah"}



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

def reroot_midpoint(root):
    """Return a new rooted tree whose root lies on the midpoint of the longest leaf-to-leaf path."""
    # undirected graph of the tree
    nodes, adj = [], {}
    def collect(n):
        nodes.append(n); adj[id(n)] = []
        for c in n.children:
            collect(c)
            adj[id(n)].append((c, c.length, c.support)); adj[id(c)].append((n, c.length, c.support))
    collect(root)
    byid = {id(n): n for n in nodes}
    lv = [n for n in nodes if not n.children]

    def distances(src):
        d, stack = {id(src): 0.0}, [src]
        while stack:
            u = stack.pop()
            for v, w, _ in adj[id(u)]:
                if id(v) not in d:
                    d[id(v)] = d[id(u)] + w; stack.append(v)
        return d
    far = max(((distances(a)[id(b)], a, b) for a in lv for b in lv), key=lambda x: x[0])
    total, a, b = far
    # path from a to b
    prev = {id(a): None}; stack = [a]
    while stack:
        u = stack.pop()
        for v, w, _ in adj[id(u)]:
            if id(v) not in prev:
                prev[id(v)] = u; stack.append(v)
    path = []; n = b
    while n is not None:
        path.append(n); n = prev[id(n)]
    path.reverse()                                     # a ... b
    half, acc = total / 2, 0.0
    for u, v in zip(path, path[1:]):
        w = next(w for x, w, _ in adj[id(u)] if x is v)
        if acc + w >= half:
            frac = half - acc                          # distance from u along the edge u-v
            sup = next(sp for x, _, sp in adj[id(u)] if x is v)
            break
        acc += w
    # rebuild a rooted tree from the new root placed on edge u-v
    def build(n, parent_node, length, support, seen):
        m = Node(); m.name, m.length, m.support = n.name, length, support
        m.parent = parent_node; seen.add(id(n))
        for v2, w2, sp2 in adj[id(n)]:
            if id(v2) not in seen:
                m.children.append(build(v2, m, w2, sp2, seen))
        return m
    new_root = Node()
    seen = {id(u), id(v)}
    cu = build(u, new_root, frac, sup, set(seen) - {id(u)} | {id(v)})
    cv = build(v, new_root, w - frac, sup, set(seen) - {id(v)} | {id(u)})
    new_root.children = [cu, cv]
    return new_root


def layout(root):
    """y coordinate per leaf in tree order, x = distance from the root."""
    y = {}
    counter = [0]
    def walk(n, x):
        n.x = x
        if not n.children:
            n.y = counter[0]; counter[0] += 1
        else:
            for c in n.children:
                walk(c, x + c.length)
            n.y = (n.children[0].y + n.children[-1].y) / 2
    walk(root, 0.0)
    return root


def main():
    tree = os.path.join(T, "sr_proteins.treefile")
    table = os.path.join(T, "sr_proteins.tsv")
    if not (os.path.isfile(tree) and os.path.isfile(table)):
        print("tree files not found in supplementary/tables"); return
    info = {r["label"]: r for r in csv.DictReader(open(table, encoding="utf-8"), delimiter="\t")}
    root = parse_newick(open(tree, encoding="utf-8").read())
    root = reroot_midpoint(root)
    layout(root)
    lv = leaves(root)
    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    fig, ax = plt.subplots(figsize=(6.5, 0.085 * len(lv) + 1.0))

    def draw(n):
        for c in n.children:
            ax.plot([n.x, n.x], [n.y, c.y], color="#444444", linewidth=0.7)
            ax.plot([n.x, c.x], [c.y, c.y], color="#444444", linewidth=0.7)
            if c.children and c.support is not None and c.support >= 95:
                ax.plot(c.x, c.y, "o", color="black", markersize=2.2)
            draw(c)
    draw(root)
    xmax = max(l.x for l in lv)
    for l in lv:
        r = info.get(l.name, {})
        sp = r.get("species", "")
        gene = r.get("gene", l.name)
        label = f"{ABBREV.get(sp, '')} {gene}".strip() if sp else l.name
        ax.text(l.x + xmax * 0.01, l.y, label, va="center", fontsize=6.5, color=SPECIES_COLOR.get(sp, "black"))
    # subfamily brackets at the right margin
    # one bracket per contiguous block of leaves with the same subfamily label (in tree order)
    ordered = sorted(lv, key=lambda l: l.y)
    blocks, subs = [], {}
    for l in ordered:
        s = info.get(l.name, {}).get("subfamily", "")
        subs.setdefault(s, []).append(l.y)
        if blocks and blocks[-1][0] == s:
            blocks[-1][2] = l.y
        else:
            blocks.append([s, l.y, l.y])
    for s, y0, y1 in blocks:
        ax.plot([xmax * 1.42, xmax * 1.42], [y0 - 0.3, y1 + 0.3], color="#888888", linewidth=1.2)
        ax.text(xmax * 1.45, (y0 + y1) / 2, s if s != "SR45" else "SR45-like", va="center", fontsize=7 if y1 - y0 > 2 else 6.5)
    ax.plot([0, 0.2], [-2, -2], color="black", linewidth=1); ax.text(0.25, -2, "0.2 substitutions per site", ha="left", va="center", fontsize=7)
    for sp, col in SPECIES_COLOR.items():
        ax.plot([], [], "s", color=col, label=SPECIES_NAME[sp])
    fig.legend(frameon=False, fontsize=7, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.005))
    ax.set_xlim(-0.02 * xmax, xmax * 1.6); ax.set_ylim(-5, len(lv))
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "svg", "tif"):
        fig.savefig(os.path.join(OUT, f"figure5_sr_tree.{ext}"), dpi=600 if ext == "tif" else 300 if ext == "png" else None, bbox_inches="tight", metadata=({"Software": None} if ext == "png" else {"Creator": None, "Date": None} if ext == "svg" else None), **({"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}))
        if ext == "tif":
            flatten_tiff(os.path.join(OUT, "figure5_sr_tree.tif"))
    print("wrote figure5_sr_tree;", len(lv), "leaves; subfamilies:", ", ".join(f"{s} {len(v)}" for s, v in sorted(subs.items())))


if __name__ == "__main__":
    main()
