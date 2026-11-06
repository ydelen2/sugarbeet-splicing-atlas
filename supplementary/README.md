# Supplementary data

Result files copied from the cluster (`results_full/` on the analysis system) that the manuscript tables, figures and
supplementary files are built from. Every file was produced by the scripts of this repository; the digest
`tables/results_digest.txt` is the printed summary from which the numbers in the text are taken.

## tables/

| File | Content | Produced by |
|---|---|---|
| `event_concordance_summary.tsv` | SUPPA2 and rMATS agreement per comparison and event type, with the overlap expected by chance and its hypergeometric p-value | `scripts/analysis/event_concordance.py` |
| `event_concordance_high_confidence.tsv` | events called by both methods with the same direction, one row per event and comparison | same |
| `events_high_confidence_annotated.tsv` | one row per distinct high-confidence event: reference gene, eggNOG description, Pfam domains, recurrence, direction | `scripts/analysis/annotate_events.py` |
| `events_high_confidence_consequences.tsv` | location of the alternative region (CDS, UTR, outside the models), reading frame, stop codons in retained introns | `scripts/analysis/event_consequences.py` |
| `events_high_confidence_nmd.tsv` | both transcript forms of every event on the annotated gene models: premature stop codons, distance to the last exon junction, 3' UTR length, NMD features, direction of the PTC form | `scripts/analysis/nmd_features.py` |
| `ptc_form_abundance_summary.tsv`, `ptc_form_abundance_events.tsv` | abundance of the premature-stop form among the non-annotated forms of the coding events, over all quantified runs | `scripts/analysis/ptc_form_abundance.py` |
| `matched_permutation.tsv` | permutation test of the SR, family and GO enrichments with random gene sets matched on testable events per gene | `scripts/analysis/matched_permutation.py` |
| `longread_event_depth_psi.tsv` | the depth-stratified long-read comparison restricted to events with a mean PSI between 0.1 and 0.9 | `scripts/analysis/longread_event_depth.py --psi` |
| `nmd_factors_summary.txt`, `nmd_factors_table.tsv`, `nmd_factors_genes.tsv` | expression of the NMD factor homologs (UPF1, UPF2, UPF3, SMG7) in every comparison | `scripts/analysis/nmd_factor_expression.py` |
| `events_high_confidence_domains.tsv` | protein coordinates of the alternative region of the coding events and the Pfam domains it overlaps | `scripts/analysis/domain_impact.py` |
| `aba_pathway_events.tsv`, `aba_pathway_enrichment.tsv` | abscisic acid pathway genes (annotation terms) with high-confidence events and their enrichment per condition group | `scripts/analysis/pathway_genes.py` |
| `de_summary.tsv`, `de_das_genes.tsv`, `de_regulators.tsv` | differential expression per comparison (PyDESeq2 on the STAR gene counts) and its overlap with the differentially spliced genes; expression of the differentially spliced genes and of the splicing-regulator genes | `scripts/analysis/differential_expression.py` |
| `cross_species_cold_summary.tsv`, `cross_species_cold_genes.tsv` | sugar beet cold splicing genes against the Arabidopsis cold DAS genes of Calixto et al. 2018 (background, cold genes, cold core) and the linked genes | `scripts/analysis/cross_species_cold.py` |
| `junction_motifs.txt` | splice-site motif classes and intron lengths of annotated and novel junctions, and junction support of the high-confidence events | `scripts/analysis/junction_motifs.py` |
| `longread_summary.txt`, `longread_events.tsv`, `longread_junctions.tsv.gz`, `longread_event_depth.tsv` | junctions and events observed in the two public PacBio cDNA runs, and the depth-stratified comparison of high-confidence and other events | `workflow/13_junction_evidence.sbatch`, `scripts/analysis/longread_junctions.py`, `longread_event_depth.py` |
| `support_summary.txt`, `support_junctions.tsv.gz` | Illumina support and motifs of the expressed-set junctions outside the pooled table | `scripts/analysis/junction_support.py` |
| `splicing_summary.tsv`, `splicing_summary_gene_overlap.tsv` | per-method counts and gene-level overlap per comparison | `scripts/analysis/summarize_splicing.py` |
| `go_enrichment.tsv` | GO term enrichment per condition group | `scripts/analysis/go_enrichment.py` |
| `splicing_regulators_proteins.tsv`, `splicing_regulators_families.tsv`, `splicing_regulators_enrichment.tsv` | splicing-regulator families of sugar beet and their enrichment of high-confidence events | `scripts/analysis/splicing_regulators.py` |
| `*.sr_candidates.tsv` | SR protein candidates per species (beet, quinoa, spinach, amaranth) | `scripts/analysis/sr_candidates.py` |
| `sr_proteins.tsv`, `sr_proteins.trim.fa`, `sr_proteins.treefile`, `sr_proteins.iqtree`, `sr_tree_summary.tsv` | SR protein tree: sequence table, trimmed alignment, IQ-TREE tree and report, tree summary | `workflow/11_sr_phylogeny.sbatch`, `scripts/analysis/tree_summary.py` |
| `star_summary.tsv`, `salmon_qc.tsv` | alignment and quantification quality per run | `workflow/05_qc_summary.sh`, `scripts/analysis/salmon_qc.py` |
| `structural_report.tsv`, `transcript_map.tsv` | structural filter counts per class code; merged transcript to reference gene map | `scripts/analysis/filter_transcripts.py` |
| `suppa_event_counts.tsv` | number of SUPPA2 events per type in the expressed transcriptome (from the SUPPA2 job log) | `workflow/07_suppa.sbatch` |
| `results_digest.txt` | printed digest of all tables | `scripts/analysis/results_digest.py` |

## tables_tissue/

The same tables for the five leaf against taproot comparisons of PRJNA602804 (`metadata/comparisons_tissue.tsv`), written by `workflow/09_summarize.sh` with `COMPARISONS_FILE` and `SUMMARY_DIR`. The testable-gene background of these tables is built from the five tissue comparisons only; they are separate from the 40 condition comparisons.

## pkg/ (annotation, event list and source files of the supplementary tables of the manuscript)

| File | Content |
|---|---|
| `annotation/bvas_structural_96221.gtf.gz` | transcript models after the structural filter (96,221 transcripts), RefBeet-3.0 coordinates |
| `annotation/bvas_expressed_73743.gtf.gz` | transcript models after the expression filter (73,743 transcripts), the annotation used for all splicing analyses |
| `annotation/structural_report.tsv` | class code counts of the structural filter |
| `junctions/pooled_junctions.tsv.gz` | splice junctions pooled from 444 runs with their run support and annotation status |
| `junctions/ref_junctions.tsv`, `junctions/merge_runs.tsv` | BeetSet-3 junctions; runs used for merging with their unique mapping rate |
| `events/suppa_all_events.ioe.gz` | the 34,363 local splicing events with their transcripts |
| `features/region_features_summary.tsv`, `features/region_features_regions.tsv.gz` | length and GC content of retained introns and skipped exons (annotated, tested, high-confidence) |
