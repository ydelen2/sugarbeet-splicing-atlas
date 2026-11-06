# Alternative splicing atlas of sugar beet (Beta vulgaris)

Code and metadata for a genome-wide analysis of alternative splicing in sugar beet built from public RNA-seq data mapped to the RefBeet-3.0 reference genome. The repository is organized so that every table and result can be regenerated from public sources.

## Repository layout

```
environment/        conda specification (bvas.yml); exported exact versions after installation
metadata/
  raw/              ENA run metadata, BioSample attributes and submitter aliases
  samples_all.tsv   465 labeled runs selected for the study (464 analyzed, one excluded)
  samples_pilot.tsv 33 runs from three BioProjects used for the pilot
scripts/metadata/   scripts that rebuild everything in metadata/ from ENA
scripts/analysis/   analysis scripts called by the workflow
scripts/figures/    scripts that build the tables and figures of the manuscript
workflow/           SLURM scripts for the Swan cluster (University of Nebraska-Lincoln)
docs/               running instructions for the Swan cluster
supplementary/      result tables, annotation files and the event list (supplementary/README.md)
manuscript/         figures and tables of the manuscript
```

## Data

All sequencing data are public and are downloaded from the European Nucleotide Archive. Run selection:

| Step | Criterion | Runs |
|---|---|---|
| Transcriptomic runs for Beta vulgaris (taxon 3555) | ENA query `tax_tree(3555) AND library_source="TRANSCRIPTOMIC"` | 756 |
| Suitable for splice junction analysis | RNA-Seq or ssRNA-seq, paired-end, mean read length at least 90 bp, at least 2 Gb per run | 465 |

The reference genome RefBeet-3.0 and the gene set BeetSet-3 are downloaded from https://bvseq.boku.ac.at (Dohm et al. 2025, https://doi.org/10.1093/nargab/lqaf142).

## Rebuilding the metadata

Requires Python 3.8 or later and internet access. ENA is a live archive; the records in `metadata/raw/` were retrieved on 17 September 2026, and a rebuild writes its own retrieval date to `metadata/raw/ena_runs_raw.retrieved.txt`.

```bash
python scripts/metadata/01_fetch_ena_runs.py
python scripts/metadata/02_filter_runs.py
python scripts/metadata/03_fetch_biosample_attributes.py
python scripts/metadata/04_fetch_run_aliases.py
python scripts/metadata/05_make_sample_sheets.py
python scripts/metadata/06_summarize_dataset.py
python scripts/metadata/07_make_comparisons.py
```

## Running the workflow

The workflow was written for the Swan cluster at the Holland Computing Center, University of Nebraska-Lincoln (SLURM, Lmod, Anaconda). Storage paths are set in `workflow/config.sh`. Steps:

| Script | Purpose |
|---|---|
| `00_storage_check.sh` | quota and read/write speed of the storage areas |
| `00_env_setup.sh` | create the conda environment from `environment/bvas.yml` |
| `01a_get_reference.sh` | download RefBeet-3.0 and BeetSet-3, convert GFF3 to GTF |
| `01b_star_index.sbatch` | STAR genome index |
| `01c_add_organelles.sh` | full analysis only: add the chloroplast (KR230391.1) and mitochondrial (NC_002511.2) genomes to the reference |
| `02_download.sbatch`, `02_download_xfer.sh` | FASTQ download from ENA with MD5 check |
| `03_align.sbatch` | fastp trimming and two-pass STAR alignment |
| `04a_stringtie.sbatch` | per-run reference-guided transcript assembly |
| `04b_merge.sbatch` | merged transcriptome, gffcompare class codes, pooled junctions |
| `04c_filter_index.sbatch` | structural transcript filter and Salmon index |
| `05_qc_summary.sh` | alignment summary and MultiQC report |
| `06_salmon_quant.sbatch` | transcript quantification against the filtered transcriptome |
| `06b_salmon_summary.sbatch` | Salmon QC, transcript TPM matrix and a scan of expression-filter thresholds |
| `07_suppa.sbatch` | expression filter, SUPPA2 events, PSI values and differential splicing per comparison |
| `07b_suppa_diff.sbatch` | SUPPA2 diffSplice for an additional comparison table (for example the leaf and taproot comparisons), reusing the existing events and PSI matrices; reports the effect of the additional groups on the expression filter |
| `08_rmats.sbatch` | rMATS junction-count based differential splicing per comparison |
| `09_summarize.sh` | summary table of both tools with a read-support filter for rMATS events |
| (same script) | `event_concordance.py`: SE, A3, A5, RI and MX events of both tools matched on splice-site coordinates, BH-adjusted SUPPA2 p-values, high-confidence event list |
| (same script) | `annotate_events.py`: reference gene by coordinate, eggNOG description and recurrence across comparisons for the high-confidence events |
| (same script) | `event_consequences.py`: CDS or UTR location, reading frame and stop codons in retained introns for the high-confidence events |
| (same script) | `splicing_regulators.py`: splicing-regulator families from the Pfam scan and the Arabidopsis searches (SR proteins by sequence rule and by best hit, U2AF, hnRNP-like, PRP, snRNP), with the high-confidence events of the members |
| (same script) | `go_enrichment.py`: GO term enrichment (hypergeometric, BH) of the genes with high-confidence events per condition group, background = genes with a testable event |
| (same script) | `nmd_features.py`: both transcript forms of every event built on the annotated gene models, premature stop codons, distance to the last exon junction and 3' UTR length (50-nucleotide rule, long 3' UTR) |
| (same script) | `ptc_form_abundance.py`: share of the non-annotated form of every coding event in the SUPPA2 PSI over all runs, premature-stop forms against forms without one (Mann-Whitney, retained introns apart); `nmd_factor_expression.py`: DESeq2 results of the UPF1, UPF2, UPF3 and SMG7 homologs in every comparison |
| (same script) | `domain_impact.py`: protein coordinates of the alternative region on the annotated protein and the Pfam domains it overlaps, for the coding events |
| (same script) | `pathway_genes.py`: abscisic acid pathway genes (annotation terms in `metadata/curation/aba_pathway_terms.tsv`) among the genes with events, with a hypergeometric test against the testable genes |
| (standalone) | `intron_features.py`: length and GC content of the retained introns and skipped exons of the high-confidence set against the tested events and the annotation (Mann-Whitney U) |
| `10_domains.sbatch` | Pfam-A domain scan of the BeetSet-3 proteome (HMMER, gathering thresholds) and DIAMOND searches against the Arabidopsis UniProt reference proteome in both directions |
| `11_sr_phylogeny.sbatch` | SR protein family across quinoa, spinach, amaranth (NCBI reference proteomes), sugar beet and Arabidopsis: same sequence rule in every proteome (`sr_candidates.py`), MAFFT, ClipKIT and IQ-TREE |
| (after 11) | `tree_summary.py`: nearest Arabidopsis SR protein of every candidate by patristic distance and the bootstrap support of each subfamily as a split of the unrooted tree |
| `12_expression.sbatch` | differential gene expression per comparison from the STAR gene counts (PyDESeq2) and its overlap with the high-confidence splicing genes (`differential_expression.py`) |
| `13_junction_evidence.sbatch` | independent junction evidence: the two public PacBio full-length cDNA runs of KWS2320 (`metadata/longread_runs.tsv`) aligned with minimap2 and compared with the annotated, novel and event junctions (`longread_junctions.py`), and the Illumina support of the expressed-set junctions outside the pooled table (`junction_support.py`) |
| (standalone) | `junction_motifs.py`: motif classes and intron lengths of annotated and novel junctions; `longread_event_depth.py`: both-forms observation of events in long reads within depth strata (Mantel-Haenszel; `--psi` restricts the events to an intermediate mean PSI); `matched_permutation.py`: enrichment of events in the SR, family and GO gene sets against random gene sets matched on the number of testable events per gene |
| (standalone) | `cross_species_cold.py`: cold splicing genes against the Arabidopsis cold DAS gene list of Calixto et al. 2018 through the DIAMOND best hits and the UniProt id mapping |

The scripts work on one of two analyses, selected with the `ANALYSIS` variable: `pilot` (33 runs, default) or `full` (465 runs, reference with organelle genomes, results in `results_full/`), for example `ANALYSIS=full sbatch workflow/03_align.sbatch`.

Comparisons are defined in `metadata/comparisons_pilot.tsv` and `metadata/comparisons_all.tsv`; the latter is built by `scripts/metadata/07_make_comparisons.py`. Input files for SUPPA2 and rMATS are built by `scripts/analysis/prepare_splicing_inputs.py`.

Detailed instructions: `docs/swan_guide.md`.

## License and citation

The code (scripts/, workflow/, environment/) is released under the MIT License (LICENSE). The curated metadata, the result tables and the generated figures and tables are released under the Creative Commons Attribution 4.0 International license (LICENSE-DATA). A manuscript describing this work is under review; the citation will be added here when it is published.

## Tables and figures

`supplementary/tables/` holds the result tables copied from the cluster (event concordance, high-confidence events with annotation and
consequences, GO enrichment, splicing-regulator families, quality summaries). `scripts/figures/make_tables.py` and
`scripts/figures/make_figures.py` build the manuscript tables (`manuscript/tables/`) and figures (`manuscript/figures/`) from these files
and the metadata, so every number in the text can be traced to a table. `results_digest.txt` in the same folder is the printed digest of
`scripts/analysis/results_digest.py`. `scripts/figures/make_tree_figure.py` draws the SR protein tree (Fig. 5). `scripts/figures/verify_numbers.py` recomputes the quantities quoted in the text from the same tables for a final check.
