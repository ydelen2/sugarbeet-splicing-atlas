"""Build the labeled sample sheets used by the Swan workflow.

Inputs (metadata/raw/, produced by scripts 01-04):
  samples_candidate.tsv        ENA runs passing the splice-analysis filter (PE, >=90 bp, >=2 Gb)
  samples_candidate_attrs.tsv  BioSample attributes per run
  run_aliases.tsv              run, experiment, library and sample aliases per run
Outputs (metadata/):
  samples_all.tsv    one row per run with condition, tissue, genotype, time point and replicate labels
  samples_pilot.tsv  pilot subset (3 BioProjects)
Labels are assigned by study-specific rules below. Rows labeled "unlabeled" still need
labels from the source publication.

Usage:
  python scripts/metadata/05_make_sample_sheets.py
"""
import csv, re, collections, os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(REPO, "metadata", "raw")
OUT = os.path.join(REPO, "metadata")


def load(name):
    with open(os.path.join(RAW, name), encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

cand = load("samples_candidate.tsv")
attrs = {r["run_accession"]: r["attributes"] for r in load("samples_candidate_attrs.tsv")}
alias = {r["run_accession"]: r for r in load("run_aliases.tsv")}

def parse_attrs(s):
    d = {}
    for kv in s.split(" | "):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip().lower()] = v.strip()
    return d

def reported_tissue(a):
    """Tissue exactly as written by the submitter (BioSample attribute), or empty."""
    return a.get("tissue") or a.get("tissue type") or a.get("tissue_type") or ""


def normalize_tissue(t):
    """Map the submitter's wording to a small controlled vocabulary."""
    t = t.lower()
    has_leaf = "leaf" in t or "leaves" in t
    if has_leaf and "root" in t:
        return "leaf+root"
    if has_leaf:
        return "leaf"
    if "taproot" in t or "tap root" in t:
        return "taproot"
    if "root" in t:
        return "root"
    if "seed" in t:
        return "seed"
    if "inflorescence" in t:
        return "inflorescence"
    return "not reported" if not t else t


def label(r, a, al):
    """Return dict(condition, tissue, genotype, rep, timepoint, category, use_das).

    Tissue always comes from the submitter's BioSample attribute; a study rule only
    fills it in when the attribute is missing (PRJEB35822)."""
    st = r["study_accession"]
    lib = al.get("library_name", "") or ""
    sa = al.get("sample_alias", "") or ""
    out = dict(condition="unlabeled", tissue=normalize_tissue(reported_tissue(a)),
               genotype=a.get("cultivar") or a.get("isolate") or a.get("genotype") or "",
               rep=a.get("replicate", ""), timepoint=a.get("dev_stage", ""), category="misc", use_das="no")
    if st == "PRJNA748559":  # cold vs control, source leaf, 15 weeks
        # Genotype code from the submitter's sample alias, e.g. 0_C_15w_Leaf_GT3_Rep2.
        # The library names of the GT3 samples say GT2 (Bv_RNA_0C_GT2-2); the sample alias is used.
        m = re.search(r"_(GT\d)_", sa)
        out.update(condition="cold" if a.get("treatment", "").startswith("cold") else "control",
                   genotype=m.group(1) if m else "", timepoint="15w",
                   category="abiotic_cold", use_das="yes", rep=a.get("replicate", ""))
    elif st == "PRJEB35822":  # salt / heat / high light, KWS2320DH
        m = lib
        tissue = "root" if "root" in m.lower() else ("leaf" if "leaf" in m.lower() else "not reported")
        cond = m.lower()
        if cond.startswith("control"):
            c = "control"
        elif "salt" in cond:
            c = "salt_300mM"
        elif "heat" in cond:
            c = "heat"
        elif "light" in cond:
            c = "high_light"
        else:
            c = "unlabeled"
        out.update(condition=c, tissue=tissue, genotype="KWS2320DH", category="abiotic_multi",
                   use_das="junction_only")  # no biological replicates -> discovery only
        out["timepoint"] = m
    elif st == "PRJNA574280":  # Sclerotinia
        out.update(condition="sclerotinia" if "infected" in a.get("infection", "") and "control" not in a.get("infection", "") else "control",
                   category="biotic_fungal", use_das="yes")
        out["rep"] = re.sub(r".*rep", "", al.get("sample_alias", ""))
    elif st == "PRJNA602804":  # 4_C_10w_Tap_GT1_Rep1
        m = re.match(r"(\d+)_C_(\d+w)_(Tap|Leaf)_(GT\d)_Rep(\d)", sa)
        if m:
            out.update(condition=("cold_4C" if m.group(1) == "4" else "control_20C"), timepoint=m.group(2),
                       genotype=m.group(4), rep=m.group(5), category="abiotic_cold", use_das="yes")
    elif st == "PRJNA1111959":  # 1V_FC_1  rep V/NV cultivar dpi
        m = re.match(r"(\d)(N?V)_(FC|El)_(\d+)", sa)
        if m:
            out.update(rep=m.group(1), condition=("BCTV" if m.group(2) == "V" else "mock"),
                       genotype=("FC709-2" if m.group(3) == "FC" else "EL-10"), timepoint=m.group(4) + "dpi",
                       category="biotic_virus", use_das="yes")
    elif st == "PRJNA1072272":  # G1_30hpi_Hokkaido_C
        m = re.match(r"(G\d+)_(\d+ ?hpi)_([A-Za-z0-9\.]+?)(NC)?_([ABC])$", sa)
        if m:
            treat = a.get("treatment", "").lower()
            inoc = treat.startswith("inoc")
            out.update(genotype=m.group(1), timepoint=m.group(2).replace(" ", ""),
                       condition=(("aphanomyces_" if inoc else "mock_") + m.group(3)), rep=m.group(5),
                       category="biotic_oomycete", use_das="yes")
    elif st == "PRJNA450098":
        t = a.get("treatment", "")
        if t:
            out.update(condition=t, category="abiotic_germination", use_das="junction_only")
        else:
            out.update(condition="dev_root_" + re.sub(r" ", "_", a.get("age", "") or a.get("dev_stage", "")),
                       category="development", use_das="junction_only")
    elif st == "PRJNA675466":
        out.update(condition=a.get("treatment", "").replace(" ", "_") + "_" + a.get("temp", ""),
                   category="abiotic_cold_dormancy", use_das="yes", rep=a.get("replicate", ""))
    elif st == "PRJNA1276275":
        out.update(condition="line_" + a.get("cultivar", ""), category="genotype_panel", use_das="junction_only")
    elif st == "PRJNA1473699":  # table beet root, YG/RG x stage
        out.update(condition=("yellow" if sa.startswith("Y") else "red") + "_" + a.get("dev_stage", "").replace(" ", "_"),
                   category="development", use_das="yes", rep=sa.split("-")[-1])
    elif st == "PRJNA343692":
        out.update(condition="dev_" + a.get("dev stage", a.get("dev_stage", "")).replace(" ", "_"),
                   category="development", use_das="junction_only")
    elif st == "PRJNA1394332":
        out.update(condition=("WT" if sa.startswith("WT") else "transgenic_" + sa.split("-")[0]),
                   category="biotic_virus_transgenic", use_das="yes", rep=sa.split("-")[-1])
    elif st == "PRJNA254489":
        c = sa.replace("Bvu_leaf_", "")
        out.update(condition=c, genotype="KWS2320", category="abiotic_multi", use_das="junction_only")
    elif st == "PRJNA610534":
        out.update(condition="unlabeled_storability", category="postharvest", use_das="pending_curation")
    else:
        out.update(condition="misc_" + st, use_das="junction_only")
    if not out["genotype"] or out["genotype"].lower() == "missing":
        out["genotype"] = "not reported"
    return out

rows = []
for r in cand:
    a = parse_attrs(attrs.get(r["run_accession"], ""))
    al = alias.get(r["run_accession"], {})
    lab = label(r, a, al)
    ftp = r["fastq_ftp"].split(";"); md5 = r["fastq_md5"].split(";")
    if len(ftp) < 2:
        continue  # PE but only one file listed -> skip
    # ENA sometimes lists an unpaired file first; take the two ending in _1/_2
    pairs = [(u, m) for u, m in zip(ftp, md5) if re.search(r"_[12]\.fastq\.gz$", u)]
    if len(pairs) != 2:
        continue
    rows.append(dict(run=r["run_accession"], study=r["study_accession"], sample=r["sample_accession"],
                     category=lab["category"], condition=lab["condition"], tissue=lab["tissue"],
                     genotype=lab["genotype"], timepoint=lab["timepoint"], rep=lab["rep"], use_das=lab["use_das"],
                     read_len=r["read_len"], gb=r["gb"], instrument=r["instrument_model"],
                     fq1="https://" + pairs[0][0], fq2="https://" + pairs[1][0], md5_1=pairs[0][1], md5_2=pairs[1][1],
                     sample_alias=al.get("sample_alias", ""), library_name=al.get("library_name", ""),
                     tissue_reported=reported_tissue(a)))

# Manual curation: documented overrides (metadata/curation/run_overrides.tsv: run, column, value, reason).
# The reason is copied to the "curation_note" column so every manual decision stays traceable.
for r in rows:
    r["curation_note"] = ""
ov_path = os.path.join(OUT, "curation", "run_overrides.tsv")
if os.path.isfile(ov_path):
    by_run = {r["run"]: r for r in rows}
    with open(ov_path, encoding="utf-8") as f:
        for o in csv.DictReader(f, delimiter="\t"):
            if o["run"] not in by_run:
                raise SystemExit(f"override for unknown run {o['run']}")
            if o["column"] not in by_run[o["run"]] or o["column"] in ("run", "curation_note"):
                raise SystemExit(f"override for unknown column {o['column']}")
            by_run[o["run"]][o["column"]] = o["value"]
            note = f"{o['column']} set to {o['value']}: {o['reason']}"
            prev = by_run[o["run"]]["curation_note"]
            by_run[o["run"]]["curation_note"] = (prev + " | " if prev else "") + note
    print("overrides applied:", sum(1 for r in rows if r["curation_note"]))

cols = list(rows[0].keys())
rows.sort(key=lambda x: (x["study"], x["condition"], x["run"]))
with open(os.path.join(OUT, "samples_all.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)

PILOT = {"PRJNA748559", "PRJNA574280", "PRJEB35822"}
with open(os.path.join(OUT, "samples_pilot.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader()
    w.writerows([r for r in rows if r["study"] in PILOT])

print("all:", len(rows), " pilot:", sum(r["study"] in PILOT for r in rows))
c = collections.Counter((r["study"], r["condition"]) for r in rows)
for (s, cnd), n in sorted(c.items()):
    print(f"{s}\t{cnd}\t{n}")
