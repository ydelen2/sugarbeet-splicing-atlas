"""Summary of the SR protein tree without a figure.

Reads the IQ-TREE Newick file (branch lengths and ultrafast bootstrap values on the internal nodes) and
the sequence table written by 11_sr_phylogeny.sbatch (label, species, gene, subfamily, length), and
reports
  1. for every non-Arabidopsis leaf, the nearest Arabidopsis SR protein by patristic distance and its
     subfamily, next to the subfamily given by the best DIAMOND hit;
  2. for every subfamily, whether all leaves labeled with it (Arabidopsis members and candidates of the
     other species) form one split of the unrooted tree, with the bootstrap support of that split; when
     they do not, the smallest split that contains them all and the number of other leaves it contains;
  3. counts of SR proteins per species and subfamily.
The tree is treated as unrooted, so the results do not depend on a root.

Usage
  python tree_summary.py --tree results_full/phylogeny/sr_proteins.treefile \
      --table results_full/phylogeny/sr_proteins.tsv --out results_full/phylogeny/sr_tree_summary.tsv
"""
import argparse
import collections
import csv
import itertools


class Node:
    __slots__ = ("name", "length", "support", "children", "parent", "x", "y")

    def __init__(self):
        self.name, self.length, self.support, self.children, self.parent = "", 0.0, None, [], None


def parse_newick(text):
    """Recursive descent parser; internal node labels are read as support values."""
    text = text.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        # label
        start = pos
        while pos < len(text) and text[pos] not in ":,)":
            pos += 1
        label = text[start:pos].strip()
        if node.children:
            if label:
                try:
                    node.support = float(label.split("/")[0])
                except ValueError:
                    node.support = None
        else:
            node.name = label
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",)":
                pos += 1
            node.length = float(text[start:pos])
        return node

    return parse_node()


def leaves(node):
    return [node] if not node.children else [l for c in node.children for l in leaves(c)]


def patristic(root):
    """Distances between all leaves through the tree (Floyd style over the tree paths)."""
    lv = leaves(root)
    index = {l.name: i for i, l in enumerate(lv)}
    # depth from root and path to root for each leaf
    def path(n):
        p = []
        while n is not None:
            p.append(n)
            n = n.parent
        return p
    paths = {l.name: path(l) for l in lv}
    dist = {}
    for a, b in itertools.combinations(lv, 2):
        pa, pb = paths[a.name], paths[b.name]
        sb = set(id(n) for n in pb)
        lca = next(n for n in pa if id(n) in sb)
        d = 0.0
        for n in pa:
            if n is lca:
                break
            d += n.length
        for n in pb:
            if n is lca:
                break
            d += n.length
        dist[(a.name, b.name)] = dist[(b.name, a.name)] = d
    return dist, [l.name for l in lv]


def splits(root):
    """Leaf sets below every internal edge (child side), with the support of the child node."""
    out = []
    def walk(n):
        for c in n.children:
            if c.children:
                out.append((frozenset(l.name for l in leaves(c)), c.support))
            walk(c)
    walk(root)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    root = parse_newick(open(a.tree, encoding="utf-8").read())
    info = {r["label"]: r for r in csv.DictReader(open(a.table, encoding="utf-8"), delimiter="\t")}
    dist, names = patristic(root)
    all_leaves = set(names)
    missing = [n for n in names if n not in info]
    if missing:
        print(f"leaves without a table entry: {missing[:5]}")
    ath = [n for n in names if info.get(n, {}).get("species") == "arabidopsis"]

    rows = []
    print("== nearest Arabidopsis SR protein per candidate ==")
    print("label\tspecies\tsubfamily_by_hit\tnearest_ath\tnearest_subfamily\tdistance\tagree")
    for n in names:
        if n in ath or n not in info:
            continue
        near = min(ath, key=lambda x: dist[(n, x)])
        sub_hit, sub_tree = info[n]["subfamily"], info[near]["subfamily"]
        rows.append(dict(label=n, species=info[n]["species"], subfamily_by_hit=sub_hit, nearest_ath=near,
                         nearest_subfamily=sub_tree, distance=round(dist[(n, near)], 4), agree="yes" if sub_hit == sub_tree else "no"))
        print(f"{n}\t{info[n]['species']}\t{sub_hit}\t{near}\t{sub_tree}\t{dist[(n, near)]:.3f}\t{'yes' if sub_hit == sub_tree else 'no'}")
    agree = sum(1 for r in rows if r["agree"] == "yes")
    print(f"agreement between best-hit subfamily and nearest Arabidopsis protein: {agree} of {len(rows)}")

    print("\n== subfamilies as splits of the unrooted tree ==")
    sp = splits(root)
    print("subfamily\tleaves\tmonophyletic\tsupport\tsmallest_split_size\tother_leaves_inside")
    for sub in sorted({info[n]["subfamily"] for n in names if n in info}):
        members = frozenset(n for n in names if info.get(n, {}).get("subfamily") == sub)
        if len(members) < 2:
            print(f"{sub}\t{len(members)}\tsingle leaf\t\t\t")
            continue
        exact = [s for s, sup in sp if s == members or s == all_leaves - members]
        if exact:
            sup = next(sup for s, sup in sp if s == members or s == all_leaves - members)
            print(f"{sub}\t{len(members)}\tyes\t{sup}\t{len(members)}\t0")
            continue
        containing = [(s, sup) for s, sup in sp if members <= s] + [(all_leaves - s, sup) for s, sup in sp if members <= all_leaves - s]
        s, sup = min(containing, key=lambda x: len(x[0]))
        others = sorted(s - members)
        print(f"{sub}\t{len(members)}\tno\t{sup}\t{len(s)}\t{len(others)}: {', '.join(others[:8])}")

    print("\n== SR proteins per species and subfamily ==")
    tab = collections.Counter((info[n]["species"], info[n]["subfamily"]) for n in names if n in info)
    species = sorted({s for s, _ in tab})
    subs = sorted({t for _, t in tab})
    print("species\t" + "\t".join(subs) + "\ttotal")
    for s in species:
        print(s + "\t" + "\t".join(str(tab[(s, t)]) for t in subs) + f"\t{sum(tab[(s, t)] for t in subs)}")

    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"table: {a.out}")


if __name__ == "__main__":
    main()
