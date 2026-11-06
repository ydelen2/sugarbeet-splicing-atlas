# Running the workflow on Swan (University of Nebraska-Lincoln HCC)

## Storage layout

The project lives on NRDStor (`BULK_ON_NRDSTOR=1` in `workflow/config.sh`). NRDStor has no purge policy but is not backed up; the code and metadata are therefore kept in this repository.

| Location | Content |
|---|---|
| `$NRDSTOR/bvas` | code (`workflow`, `scripts`, `environment`, `metadata`), conda environment, reference genome and STAR index, logs, QC, results |
| `$NRDSTOR/bvas/bulk` | FASTQ and BAM files |
| node-local `$TMPDIR`, or `$NRDSTOR/bvas/bulk/tmp` when it is too small | temporary files of fastp and STAR |

All commands are run from `$NRDSTOR/bvas`. The sbatch scripts read `workflow/config.sh` from the directory the job was submitted from. The project path is `$NRDSTOR/bvas` (`/mnt/nrdstor/<group>/<user>/bvas`).

## 0. Copying the files to Swan

From a local terminal:
```bash
scp -r workflow metadata environment scripts <user>@swan.unl.edu:/mnt/nrdstor/<group>/<user>/bvas/
```
On Swan:
```bash
cd $NRDSTOR/bvas && chmod +x workflow/*.sh && ls workflow metadata environment scripts
```
Do not overwrite a script while a job is running it; bash reads the file as it runs.

## 1. Storage check
```bash
bash workflow/00_storage_check.sh
```
Prints the group quota and the write speed of the two storage areas.

## 2. Conda environment (once)
```bash
bash workflow/00_env_setup.sh
```
The last lines print the tool versions. The exact versions used in the study are in `environment/bvas.exported.yml`.

## 3. Reference genome, GTF and STAR index
```bash
bash workflow/01a_get_reference.sh
sbatch workflow/01b_star_index.sbatch
```
For the full analysis the organellar genomes are appended with `workflow/01c_add_organelles.sh` before the index is built (see `config.sh`, `ANALYSIS=full`).

## 4. FASTQ download

The files are written to `$NRDSTOR/bvas/bulk/fastq`. Check first whether the compute nodes can reach the archive:
```bash
srun --time=00:05:00 --mem=1G --pty bash -c 'curl -sI https://ftp.sra.ebi.ac.uk | head -1'
```
If the response code is 200:
```bash
N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%8 workflow/02_download.sbatch
```
Otherwise, on the transfer node in the background:
```bash
nohup bash workflow/02_download_xfer.sh metadata/samples_pilot.tsv > logs/download_xfer.log 2>&1 &
```
Number of completed runs:
```bash
ls $NRDSTOR/bvas/bulk/fastq/*/.download_ok | wc -l
```

## 5. Alignment
```bash
N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%20 workflow/03_align.sbatch
```
When the jobs have finished:
```bash
bash workflow/05_qc_summary.sh
```
The share of junction reads on annotated junctions is computed from the first-pass junction file (`SJ.pass1.out.tab`); the "Annotated (sjdb)" value of `Log.final.out` is not used because in two-pass mode it counts the junctions found in the first pass as annotated. The summary is `qc/star_summary.tsv`.

## 6. Transcript discovery
```bash
N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%30 workflow/04a_stringtie.sbatch
```
After these jobs:
```bash
sbatch workflow/04b_merge.sbatch
```
The class-code summary is `results/merged/class_codes.txt`. Then the structural filter and the Salmon index:
```bash
sbatch workflow/04c_filter_index.sbatch
```
The filter keeps reference transcripts (gffcompare classes `=`, `k`), multi-exon new isoforms of known genes (`j`, `m`, `n`, `o`) and multi-exon new loci (`u`); single-exon new transcripts, antisense (`x`) and the other classes are removed. The summary is `results/filtered/structural_report.tsv`.

## 7. Isoform quantification (Salmon)
```bash
N=$(($(wc -l < metadata/samples_pilot.tsv)-1)); sbatch --array=1-${N}%20 workflow/06_salmon_quant.sbatch
```
The mapping rate of each run is the `percent_mapped` line of its log.

## 8. Differential splicing

Comparisons are defined in `metadata/comparisons_pilot.tsv` (pilot) and `metadata/comparisons_all.tsv` (full analysis).

SUPPA2 (one job, after Salmon). The job first applies the expression filter (`MIN_TPM`, `MIN_RUNS` and `EXPR_FILTER` in `config.sh`) and writes `results/filtered/final.gtf`:
```bash
sbatch workflow/07_suppa.sbatch
```
rMATS (uses the same `final.gtf`, so after 07; one task per comparison):
```bash
N=$(($(wc -l < metadata/comparisons_pilot.tsv)-1)); sbatch --array=1-${N} workflow/08_rmats.sbatch
```
Summary tables, event annotation, consequences, families, GO enrichment and digest:
```bash
bash workflow/09_summarize.sh
```
rMATS events count only when every replicate has at least 20 supporting reads. The direction of the SUPPA2 dPSI is checked against the group PSI means.

## 9. Full analysis

Set `ANALYSIS=full` (465 runs, reference with the organellar genomes, results in `results_full/`) and repeat steps 3 to 8 with `metadata/samples_all.tsv` and `metadata/comparisons_all.tsv`, for example `ANALYSIS=full sbatch workflow/03_align.sbatch`. FASTQ files are deleted after Salmon (`DELETE_FASTQ_AFTER_QUANT=1`), because STAR and Salmon both read them. The full set is 465 runs and 2,339 GB of compressed FASTQ.

Protein domains and splicing-regulator families (Pfam, HMMER, DIAMOND against the Arabidopsis proteome):
```bash
ANALYSIS=full sbatch workflow/10_domains.sbatch
```
SR protein phylogeny with the proteomes of quinoa, spinach and amaranth:
```bash
ANALYSIS=full sbatch workflow/11_sr_phylogeny.sbatch
```
Gene expression (PyDESeq2) and its overlap with differential splicing; the package is installed once:
```bash
conda install -p "$NRDSTOR/bvas/envs/bvas" -c conda-forge -c bioconda pydeseq2
```
```bash
ANALYSIS=full sbatch workflow/12_expression.sbatch
```
Leaf and taproot comparisons (`metadata/comparisons_tissue.tsv`), written to a separate folder without touching the main tables:
```bash
COMPARISONS_FILE=metadata/comparisons_tissue.tsv ANALYSIS=full sbatch workflow/07b_suppa_diff.sbatch
```
```bash
COMPARISONS_FILE=metadata/comparisons_tissue.tsv ANALYSIS=full sbatch --array=1-5 workflow/08_rmats.sbatch
```
```bash
COMPARISONS_FILE=metadata/comparisons_tissue.tsv SUMMARY_DIR=$NRDSTOR/bvas/results_full/summary_tissue ANALYSIS=full bash workflow/09_summarize.sh
```
Junction evidence (two public PacBio full-length cDNA runs and the Illumina support of the junctions outside the pooled table); minimap2 is installed once:
```bash
conda install -p "$NRDSTOR/bvas/envs/bvas" -c conda-forge -c bioconda minimap2
```
```bash
ANALYSIS=full sbatch workflow/13_junction_evidence.sbatch
```
The results are written to `results_full/junction_evidence/` (`longread_summary.txt`, `longread_events.tsv`, `support_summary.txt`).

Comparison with the cold response of Arabidopsis. Two external files are needed: the cold DAS gene list of Calixto et al. 2018 (Plant Cell, https://doi.org/10.1105/tpc.18.00177), Supplemental Data Set 1 (text or CSV; the script reads the AGI codes in the file) and the UniProt id mapping of the Arabidopsis reference proteome (`UP000006548_3702.idmapping.gz`, UniProt reference_proteomes). The lists used in the study are in `metadata/external/`.
```bash
ANALYSIS=full bash -c 'source workflow/config.sh; activate_env; python scripts/analysis/cross_species_cold.py --events results_full/summary/events_high_confidence_annotated.tsv --suppa-dir results_full/suppa --comparisons metadata/comparisons_all.tsv --map results_full/filtered/transcript_map.tsv --beet-ath results_full/domains/beet_vs_ath.tsv --idmapping ref/proteomes/UP000006548_3702.idmapping.gz --ath-list metadata/external/calixto2018_cold_das_genes.txt --ath-expressed metadata/external/calixto2018_expressed_genes.txt --out-prefix results_full/summary/cross_species_cold'
```

Two further checks run inside `09_summarize.sh` when their inputs exist: the abundance of the premature-stop form of every coding event over all runs (`ptc_form_abundance.py`, reads `results_full/suppa/psi/all_events.psi`, writes `summary/ptc_form_abundance_summary.tsv`) and the expression of the NMD factors UPF1, UPF2, UPF3 and SMG7 (`nmd_factor_expression.py`, needs `results_full/expression/de_results` and the id mapping, writes `summary/nmd_factors_summary.txt`). To run only these two steps:
```bash
ANALYSIS=full bash -c 'source workflow/config.sh; activate_env; S=results_full/summary; python scripts/analysis/ptc_form_abundance.py --nmd $S/events_high_confidence_nmd.tsv --consequences $S/events_high_confidence_consequences.tsv --psi results_full/suppa/psi/all_events.psi --out-prefix $S/ptc_form_abundance; python scripts/analysis/nmd_factor_expression.py --beet-ath results_full/domains/beet_vs_ath.tsv --idmapping ref/proteomes/UP000006548_3702.idmapping.gz --de-dir results_full/expression/de_results --comparisons metadata/comparisons_all.tsv --out-prefix $S/nmd_factors'
```

Two matched controls of the enrichment and long-read results:
```bash
ANALYSIS=full bash -c 'source workflow/config.sh; activate_env; S=results_full/summary; python scripts/analysis/matched_permutation.py --suppa-dir results_full/suppa --comparisons metadata/comparisons_all.tsv --map results_full/filtered/transcript_map.tsv --events $S/events_high_confidence_annotated.tsv --families results_full/domains/splicing_regulators_families.tsv --functional ref/BeetSet-3.genes.2604.functional --go-term GO:0008380 --go-term GO:0000398 --out $S/matched_permutation.tsv; python scripts/analysis/longread_event_depth.py --junctions results_full/junction_evidence/longread_junctions.tsv --ioe results_full/suppa/events/all_events.ioe --events $S/events_high_confidence_annotated.tsv --psi results_full/suppa/psi/all_events.psi --psi-min 0.1 --psi-max 0.9 --out results_full/junction_evidence/longread_event_depth_psi.tsv'
```
The first repeats the SR, family and GO enrichment tests with random gene sets matched on the number of testable events per gene (`matched_permutation.tsv`); the second restricts the long-read comparison to events with a mean PSI between 0.1 and 0.9 (`longread_event_depth_psi.tsv`).

## Result tables

The tables that the manuscript tables and figures are built from are copied from `results_full/summary/`, `results_full/domains/`, `results_full/expression/` and `results_full/junction_evidence/` into `supplementary/tables/` (and `results_full/summary_tissue/` into `supplementary/tables_tissue/`), as listed in `supplementary/README.md`. `scripts/figures/make_tables.py`, `make_figures.py`, `make_supplementary_figures.py` and `make_tree_figure.py` then regenerate the tables and figures, and `scripts/figures/verify_numbers.py` recomputes the quantities quoted in the text.
