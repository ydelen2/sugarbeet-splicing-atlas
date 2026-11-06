#!/bin/bash
#SBATCH --job-name=bv_sr_tree
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
# SR protein family across Amaranthaceae and Arabidopsis: reference proteomes of quinoa, spinach and
# amaranth from NCBI Datasets (resolved from the species name, logged), the same SR rule in every
# proteome (RRM_1 HMM from Pfam-A plus the RS domain test), best hits to the Arabidopsis UniProt
# proteome, and a maximum likelihood tree of the SR proteins (MAFFT L-INS-i, ClipKIT, IQ-TREE with
# ModelFinder and 1000 ultrafast bootstrap replicates). Run after 10_domains.sbatch.
# Usage (from $NRDSTOR/bvas):  ANALYSIS=full sbatch workflow/11_sr_phylogeny.sbatch
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-.}/workflow/config.sh"
# The tool modules depend on the compiler modules of the default environment, which activate_env purges;
# they are loaded first and the conda environment is activated afterwards without a purge.
module load compiler/gcc/11 openmpi/4.1          # toolchain the tool modules are built against (module spider hmmer/3.4)
module load hmmer/3.4 diamond/2.1 mafft/7.526 iqtree/2.3 clipkit/2.7
set +u; module load anaconda; conda activate "$ENV_DIR"; set -u
for t in hmmsearch hmmfetch diamond mafft clipkit; do command -v "$t" > /dev/null || { echo "$t not found after module load" >&2; module list 2>&1; exit 1; }; done
command -v iqtree2 > /dev/null || command -v iqtree > /dev/null || { echo "iqtree not found after module load" >&2; exit 1; }
T="${SLURM_CPUS_ON_NODE:-16}"
A="${SCRIPTS_DIR}/analysis"
PROT_DIR="${REF_DIR}/proteomes"; DOM="${RES_DIR}/domains"; OUT="${RES_DIR}/phylogeny"
mkdir -p "$OUT" "$PROT_DIR"
ATH_AA="${PROT_DIR}/athaliana_UP000006548.fa"
FAM="${META_DIR}/curation/splicing_factor_families.tsv"
API="https://api.ncbi.nlm.nih.gov/datasets/v2"
for f in "$ATH_AA" "$DOM/beet_proteins.fa" "${REF_DIR}/pfam/Pfam-A.hmm"; do
  [ -s "$f" ] || { echo "missing $f (run 10_domains.sbatch)" >&2; exit 1; }
done

echo "== versions =="; hmmsearch -h | grep -m1 HMMER; diamond version; mafft --version 2>&1 | head -n 1; iqtree2 --version 2>/dev/null | head -n 1 || iqtree --version | head -n 1; clipkit --version

echo; echo "== reference proteomes from NCBI Datasets =="
declare -A TAXA=( [quinoa]="Chenopodium quinoa" [spinach]="Spinacia oleracea" [amaranth]="Amaranthus hypochondriacus" )
SPECIES=()
for sp in quinoa spinach amaranth; do
  D="${PROT_DIR}/${sp}"; mkdir -p "$D"
  if [ ! -s "$D/protein.faa" ]; then
    taxon=$(printf '%s' "${TAXA[$sp]}" | sed 's/ /%20/g')
    curl -sS "${API}/genome/taxon/${taxon}/dataset_report?filters.reference_only=true" > "$D/report.json"
    ACC=$(python - "$D/report.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
reps = d.get("reports", [])
for r in reps:
    ai = r.get("annotation_info", {})
    print(r["accession"], r["organism"]["organism_name"], r.get("assembly_info", {}).get("assembly_name", ""),
          "annotation:", ai.get("name", "none"), ai.get("release_date", ""), file=sys.stderr)
ok = [r["accession"] for r in reps if r.get("annotation_info")]
print(ok[0] if ok else "")
PY
)
    if [ -z "$ACC" ]; then echo "$sp: no annotated reference assembly found, skipped"; continue; fi
    echo "$sp: downloading $ACC (protein FASTA and GFF3)"
    curl -sS -o "$D/dataset.zip" "${API}/genome/accession/${ACC}/download?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF"
    unzip -qo "$D/dataset.zip" -d "$D/unzipped"
    cp "$D/unzipped/ncbi_dataset/data/${ACC}/protein.faa" "$D/protein.faa"
    cp "$D/unzipped/ncbi_dataset/data/${ACC}/genomic.gff" "$D/genomic.gff"
    echo "$ACC" > "$D/accession.txt"; date -u +"downloaded %Y-%m-%d" > "$D/download_date.txt"
    rm -rf "$D/unzipped" "$D/dataset.zip"
  fi
  echo "$sp: $(cat "$D/accession.txt") $(cat "$D/download_date.txt"), proteins: $(grep -c '>' "$D/protein.faa")"
  SPECIES+=("$sp")
done

echo; echo "== RRM_1 profile and SR candidates per proteome =="
# hmmfetch indexes the pressed binary when it exists (Pfam-A.hmm.h3m.ssi), otherwise the text file (Pfam-A.hmm.ssi)
[ -s "${REF_DIR}/pfam/Pfam-A.hmm.h3m.ssi" ] || [ -s "${REF_DIR}/pfam/Pfam-A.hmm.ssi" ] || hmmfetch --index "${REF_DIR}/pfam/Pfam-A.hmm" > /dev/null
hmmfetch "${REF_DIR}/pfam/Pfam-A.hmm" RRM_1 > "$OUT/RRM_1.hmm"
grep -c "^NAME" "$OUT/RRM_1.hmm" | sed "s/^/profiles fetched: /"
[ -s "$DOM/ath.dmnd" ] || diamond makedb --in "$ATH_AA" -d "$DOM/ath" --quiet
FMT="qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore"
scan() { # name fasta [gff]
  local sp="$1" fa="$2" gff="${3:-}"
  hmmsearch --cpu "$T" --cut_ga --noali --domtblout "$OUT/${sp}.rrm.domtblout" "$OUT/RRM_1.hmm" "$fa" > /dev/null
  diamond blastp -q "$fa" -d "$DOM/ath" -o "$OUT/${sp}.vs_ath.tsv" --very-sensitive --evalue 1e-5 \
    --max-target-seqs 5 --threads "$T" --outfmt 6 $FMT --quiet
  python "$A/sr_candidates.py" --species "$sp" --proteins "$fa" --domtbl "$OUT/${sp}.rrm.domtblout" \
    --diamond "$OUT/${sp}.vs_ath.tsv" --families "$FAM" ${gff:+--gff "$gff"} --out "$OUT/${sp}.sr_candidates.tsv"
}
scan beet "$DOM/beet_proteins.fa"
for sp in "${SPECIES[@]}"; do scan "$sp" "${PROT_DIR}/${sp}/protein.faa" "${PROT_DIR}/${sp}/genomic.gff"; done

echo; echo "== tree of the SR proteins (Arabidopsis entries of the family table plus the candidates with an Arabidopsis SR best hit) =="
python - "$ATH_AA" "$FAM" "$OUT" beet "${SPECIES[@]}" <<'PY'
import csv, sys
ath, fam, out, species = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
entries = {r["uniprot_entry"]: r for r in csv.DictReader(open(fam, encoding="utf-8"), delimiter="\t") if r["family"] in ("SR", "SR-like")}
seqs, name = {}, None
for line in open(ath, encoding="utf-8"):
    if line.startswith(">"):
        parts = line[1:].split("|")
        name = parts[2].split()[0] if len(parts) == 3 else None
        if name in entries:
            seqs[name] = []
        else:
            name = None
    elif name:
        seqs[name].append(line.strip())
n = 0
with open(f"{out}/sr_proteins.fa", "w", newline="\n") as fh, open(f"{out}/sr_proteins.tsv", "w", newline="\n") as tab:
    tab.write("label\tspecies\tgene\tsubfamily\tlength\n")
    for e, s in seqs.items():
        seq = "".join(s); label = f"Ath_{entries[e]['gene']}"
        fh.write(f">{label}\n{seq}\n"); tab.write(f"{label}\tarabidopsis\t{entries[e]['gene']}\t{entries[e]['subfamily']}\t{len(seq)}\n"); n += 1
    for sp in species:
        for r in csv.DictReader(open(f"{out}/{sp}.sr_candidates.tsv", encoding="utf-8"), delimiter="\t"):
            if r["subfamily"] == "SR-related":
                continue
            label = f"{sp[:3].capitalize()}_{r['gene']}".replace("|", "_")
            fh.write(f">{label}\n{r['sequence']}\n"); tab.write(f"{label}\t{sp}\t{r['gene']}\t{r['subfamily']}\t{r['length']}\n"); n += 1
print(f"sequences in the tree: {n} (Arabidopsis {len(seqs)})")
PY
cd "$OUT"
mafft --localpair --maxiterate 1000 --thread "$T" --quiet sr_proteins.fa > sr_proteins.aln.fa
clipkit sr_proteins.aln.fa -m smart-gap -o sr_proteins.trim.fa -l
IQ=$(command -v iqtree2 || command -v iqtree)
"$IQ" -s sr_proteins.trim.fa -m MFP -B 1000 -T AUTO --threads-max "$T" --prefix sr_proteins -quiet -redo   # -redo: ignore the checkpoint of an earlier run
grep -m1 'Best-fit model' sr_proteins.iqtree || true
alen() { awk '/^>/{if(s){print length(s); exit}} !/^>/{s=s $0}' "$1"; }
echo "alignment columns: $(alen sr_proteins.aln.fa); after trimming: $(alen sr_proteins.trim.fa)"
echo "outputs: $OUT/sr_proteins.treefile  $OUT/sr_proteins.tsv  $OUT/*.sr_candidates.tsv"
