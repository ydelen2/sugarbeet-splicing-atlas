"""Recompute the numbers quoted in the manuscript from the local tables, for a final check of the text.

Prints every quantity that the text states, computed from supplementary/tables, supplementary/pkg and the
metadata, so the manuscript can be compared line by line with its sources.

Usage (from the repository root):  python scripts/figures/verify_numbers.py
"""
import collections
import csv
import gzip
import os
import re
import statistics

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
T = os.path.join(ROOT, "supplementary", "tables")
P = os.path.join(ROOT, "supplementary", "pkg")
META = os.path.join(ROOT, "metadata")


def read(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def main():
    print("== data set ==")
    samples = read(os.path.join(META, "samples_all.tsv"))
    used = [r for r in samples if r["use_das"] != "excluded"]
    print(f"runs in sheet {len(samples)}, analyzed {len(used)}, BioProjects analyzed {len({r['study'] for r in used})}, "
          f"eligible for comparisons {sum(1 for r in used if r['use_das'] == 'yes')}, total Gb {sum(float(r['gb']) for r in used):.1f}")
    star = read(os.path.join(T, "star_summary.tsv"))
    u = [float(r["uniq_pct"]) for r in star]
    print(f"STAR runs {len(u)}, unique mean {statistics.mean(u):.1f}, min {min(u)}, max {max(u)}, below 80: {sum(1 for x in u if x < 80)}, "
          f"annotated junction read pct min {min(float(r['annotated_junction_read_pct']) for r in star)}")
    sal = read(os.path.join(T, "salmon_qc.tsv"))
    m = [float(r["percent_mapped"]) for r in sal]
    low = [r for r in sal if float(r["percent_mapped"]) < 50]
    print(f"Salmon runs {len(m)}, median {statistics.median(m):.1f}, min {min(m)}, below 50: {len(low)} ({', '.join(r['run'] + ' ' + r['condition'] for r in low)})")
    print(f"libtype agreement: {sum(1 for r in sal if r['agree'] == 'yes')} of {len(sal)}")

    print("\n== merge and junctions ==")
    mr = read(os.path.join(P, "junctions", "merge_runs.tsv"))
    print(f"runs merged {sum(1 for r in mr if r['used'] == 'yes')} of {len(mr)}")
    ann = nov = 0; support = collections.Counter()
    with gzip.open(os.path.join(P, "junctions", "pooled_junctions.tsv.gz"), "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[-1] == "annotated":
                ann += 1
            else:
                nov += 1; support[int(f[5])] += 1
    print(f"pooled junctions: annotated {ann}, novel {nov}, novel in >= 50 runs {sum(v for k, v in support.items() if k >= 50)}, "
          f"reference junctions {sum(1 for _ in open(os.path.join(P, 'junctions', 'ref_junctions.tsv'), encoding='utf-8'))}")
    print(f"pooled header: {header}")

    print("\n== transcript sets ==")
    tmap = read(os.path.join(T, "transcript_map.tsv"))
    kept = [r for r in tmap if r["kept_structural"] == "yes"]
    print(f"merged transcripts {len(tmap)}, structural kept {len(kept)}, with reference gene {sum(1 for r in kept if r['ref_gene'])}, "
          f"loci {len({r['gene'] for r in tmap})}, novel loci {len({r['gene'] for r in tmap if not r['ref_gene']})}")
    expressed = set()
    with gzip.open(os.path.join(P, "annotation", "bvas_expressed_73743.gtf.gz"), "rt", encoding="utf-8") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "transcript":
                expressed.add(re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
    cls = collections.Counter(r["class_code"] for r in tmap if r["transcript"] in expressed)
    print(f"expressed transcripts {len(expressed)}: " + ", ".join(f"{k} {v}" for k, v in sorted(cls.items(), key=lambda x: -x[1])))
    non_ref = len(expressed) - cls["="]
    at_genes = len(expressed) - cls["u"]
    print(f"additional isoforms {non_ref}; transcripts at annotated genes {at_genes}; per annotated gene (28,111) {at_genes / 28111:.2f}")
    ev = collections.Counter()
    with gzip.open(os.path.join(P, "events", "suppa_all_events.ioe.gz"), "rt", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            ev[line.split("\t")[2].split(";")[1].split(":")[0]] += 1
    print(f"events {sum(ev.values())}: " + ", ".join(f"{k} {v}" for k, v in sorted(ev.items(), key=lambda x: -x[1])))

    print("\n== tree ==")
    iq = open(os.path.join(T, "sr_proteins.iqtree"), encoding="utf-8").read()
    print(re.search(r"Best-fit model according to BIC: (\S+)", iq).group(1),
          re.search(r"Input data: (\d+) sequences with (\d+) amino-acid sites", iq).groups())
    print("species:", dict(collections.Counter(r["species"] for r in read(os.path.join(T, "sr_proteins.tsv")))))

    print("\n== see results_digest.txt for the concordance, consequence, family and GO numbers ==")


if __name__ == "__main__":
    main()
