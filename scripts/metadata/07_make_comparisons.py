"""Build the comparison table for differential splicing from the labeled sample sheet.

Every comparison is a control-versus-treatment contrast inside one BioProject, within one genotype,
tissue and time point. Only runs with use_das = yes take part. A comparison is written only when both
groups have at least MIN_REPS runs.

Design notes
  PRJNA1072272  Aphanomyces cochlioides, roots. Non-inoculated controls were collected separately for
                each isolate (Rossi et al. 2024, https://doi.org/10.1186/s12870-024-05910-y), so an
                inoculated group is only compared with the control of the same isolate. For genotypes G12
                and G17 the archive holds controls for isolate Arhill2012 only.
  PRJNA1111959  Beet curly top virus, leaf, two genotypes, 1, 7 and 14 dpi.
  PRJNA602804   4 C versus 20 C, leaf and taproot, three genotypes.
  PRJNA748559   0 C versus 20 C, leaf; genotype GT1 only (two cold GT3 libraries are unstranded).
  PRJNA675466   Cold-induced secondary dormancy of seeds (Hourston et al. 2022, doi:10.1111/pce.14264):
                5 C versus 10 C at 60 h, the 5 C time course versus 60 h at 5 C, the physiological time
                points (1 percent germination: 272 h at 5 C, 82 h at 10 C), the 10 C germination course
                and 60 h imbibition at 10 C versus dry seed.
  PRJNA1473699  Table beet root, tuber enlargement versus rosette stage, per color type.
  PRJNA1394332  Rhizomania-resistant transgenic lines versus wild type (two replicates per group).
  PRJNA574280   Sclerotinia sclerotiorum lesion periphery versus uninfected leaf.

Usage:
  python scripts/metadata/07_make_comparisons.py
Output:
  metadata/comparisons_all.tsv
"""
import collections
import csv
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
META = os.path.join(REPO, "metadata")
MIN_REPS = 2


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def main():
    with open(os.path.join(META, "samples_all.tsv"), encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["use_das"] == "yes"]
    by = collections.defaultdict(list)
    for r in rows:
        by[r["study"]].append(r)

    comps = []

    def add(cid, study, subset, ref, test, note):
        sel = by[study]
        for col, val in subset:
            sel = [r for r in sel if r[col] == val]
        n_ref = sum(r["condition"] == ref for r in sel)
        n_test = sum(r["condition"] == test for r in sel)
        if n_ref < MIN_REPS or n_test < MIN_REPS:
            print(f"skipped {cid}: reference {n_ref}, test {n_test} runs")
            return
        comps.append(dict(comparison_id=cid, study=study,
                          subset=";".join(f"{c}={v}" for c, v in subset) or "-",
                          group_column="condition", reference=ref, test=test,
                          n_reference=n_ref, n_test=n_test, note=note))

    # Aphanomyces: genotype x time point x isolate, matched control only
    st = "PRJNA1072272"
    cells = sorted(set((r["genotype"], r["timepoint"], r["condition"]) for r in by[st]))
    for g, tp, cond in cells:
        if cond.startswith("aphanomyces_"):
            iso = cond[len("aphanomyces_"):]
            add(f"aphano_{g}_{tp}_{slug(iso)}", st, [("genotype", g), ("timepoint", tp)], f"mock_{iso}", cond,
                "root; control collected with the same isolate")

    st = "PRJNA1111959"
    for g, tp in sorted(set((r["genotype"], r["timepoint"]) for r in by[st])):
        add(f"bctv_{slug(g)}_{tp}", st, [("genotype", g), ("timepoint", tp)], "mock", "BCTV", "leaf")

    st = "PRJNA602804"
    for g, tis in sorted(set((r["genotype"], r["tissue"]) for r in by[st])):
        add(f"cold4C_{g}_{tis}", st, [("genotype", g), ("tissue", tis)], "control_20C", "cold_4C", "10 weeks")

    add("cold0C_GT1_leaf", "PRJNA748559", [("genotype", "GT1")], "control", "cold", "15 weeks, all libraries reverse-stranded")

    st = "PRJNA675466"
    add("imbibition_5C_vs_10C_60h", st, [], "imbibed_for_60_hours_10C", "imbibed_for_60_hours_5C", "seed")
    for h in ("160", "272", "432"):
        add(f"imbibition_5C_{h}h_vs_60h", st, [], "imbibed_for_60_hours_5C", f"imbibed_for_{h}_hours_5C", "seed")
    add("imbibition_5C_272h_vs_10C_82h", st, [], "imbibed_for_82_hours_10C", "imbibed_for_272_hours_5C",
        "seed, physiological time point (1 percent germination) at each temperature")
    add("imbibition_10C_82h_vs_60h", st, [], "imbibed_for_60_hours_10C", "imbibed_for_82_hours_10C",
        "seed, germination course at 10 C")
    add("imbibition_10C_60h_vs_dry", st, [], "none_20C", "imbibed_for_60_hours_10C", "seed, dry versus imbibed")

    st = "PRJNA1473699"
    for color in ("red", "yellow"):
        add(f"tablebeet_{color}_enlargement_vs_rosette", st, [], f"{color}_rosette_stage",
            f"{color}_tuber_enlargement_stage", "table beet root")

    st = "PRJNA1394332"
    for line in ("S3", "S6"):
        add(f"transgenic_{line}_vs_WT", st, [], "WT", f"transgenic_{line}", "root; two replicates per group")

    add("sclerotinia", "PRJNA574280", [], "control", "sclerotinia",
        "leaf; 59 to 69 percent of reads in infected samples are fungal, low power")

    cols = ["comparison_id", "study", "subset", "group_column", "reference", "test", "n_reference", "n_test", "note"]
    ids = [c["comparison_id"] for c in comps]
    assert len(ids) == len(set(ids)), "duplicate comparison ids"
    with open(os.path.join(META, "comparisons_all.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(comps)
    print(f"comparisons written: {len(comps)}")
    for study, n in collections.Counter(c["study"] for c in comps).items():
        print(f"  {study}\t{n}")
    used = set()
    for c in comps:
        sel = by[c["study"]]
        if c["subset"] != "-":
            for cond in c["subset"].split(";"):
                k, v = cond.split("=", 1)
                sel = [r for r in sel if r[k] == v]
        used |= {r["run"] for r in sel if r["condition"] in (c["reference"], c["test"])}
    print(f"runs used in at least one comparison: {len(used)} of {len(rows)} eligible")


if __name__ == "__main__":
    main()
